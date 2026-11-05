"""Minimal S3 client: SigV4 signed with the standard library only, so the
project gains no new dependency (CLAUDE.md section 7). Path-style addressing,
which is what s3.colgate.edu (Ceph RGW) serves.

Credentials come from the environment first -- that is how cluster pods get
them, from the s3-credentials secret, under AWS_* names -- and fall back to
Implementation/.env for local use. Credential values are never printed.

  python s3.py ls [prefix]
  python s3.py put <local_file> <key>
  python s3.py get <key> <local_file>
  python s3.py putdir <local_dir> <key_prefix>
  python s3.py getdir <key_prefix> <local_dir>
"""
import datetime, hashlib, hmac, os, sys, urllib.error, urllib.parse, urllib.request

REGION, SERVICE = "us-east-1", "s3"
CONTENT_TYPE = "application/octet-stream"
_UNRESERVED_SAFE = "/~_-."


def _load_cfg():
    alt = {"s3_endpoint": ("S3_ENDPOINT", "AWS_ENDPOINT_URL"),
           "s3_access":   ("S3_ACCESS", "AWS_ACCESS_KEY_ID"),
           "s3_key":      ("S3_KEY", "AWS_SECRET_ACCESS_KEY"),
           "s3_bucket":   ("S3_BUCKET",)}
    cfg = {k: next((os.environ[n] for n in names if os.environ.get(n)), "")
           for k, names in alt.items()}
    if all(cfg.values()):
        if not cfg["s3_endpoint"].startswith("http"):
            cfg["s3_endpoint"] = "https://" + cfg["s3_endpoint"]
        return cfg
    here = os.path.dirname(os.path.abspath(__file__))
    env = os.environ.get("S3_ENV_FILE") or os.path.join(here, "..", "..", ".env")
    if os.path.exists(env):
        with open(env, encoding="utf-8-sig") as fh:
            for ln in fh:
                ln = ln.strip()
                if ln and not ln.startswith("#") and "=" in ln:
                    k, v = ln.split("=", 1)
                    k = k.strip().lower()
                    if k in cfg and not cfg[k]:
                        cfg[k] = v.strip().strip('"').strip("'")
    missing = [k for k, v in cfg.items() if not v]
    if missing:
        sys.exit("missing S3 settings: %s" % ", ".join(sorted(missing)))
    if not cfg["s3_endpoint"].startswith("http"):
        cfg["s3_endpoint"] = "https://" + cfg["s3_endpoint"]
    return cfg


CFG = _load_cfg()
HOST = urllib.parse.urlparse(CFG["s3_endpoint"]).netloc
SCHEME = urllib.parse.urlparse(CFG["s3_endpoint"]).scheme


def _sign(key, msg):
    return hmac.new(key, msg.encode(), hashlib.sha256).digest()


def request(method, key="", query="", body=b""):
    """Signed request. `body` is bytes; SigV4 needs its hash before sending."""
    now = datetime.datetime.now(datetime.timezone.utc)
    amzdate = now.strftime("%Y%m%dT%H%M%SZ")
    datestamp = now.strftime("%Y%m%d")
    phash = hashlib.sha256(body).hexdigest()

    uri = "/" + urllib.parse.quote(CFG["s3_bucket"]) + "/"
    if key:
        uri += urllib.parse.quote(key, safe=_UNRESERVED_SAFE)

    # Ceph RGW refuses, with AccessDenied, any request carrying a Content-Type
    # that is not listed in SignedHeaders -- the value does not matter, only
    # that it is signed. urllib adds a Content-Type to every request with a
    # body, so leaving it unsigned makes every write fail while reads, which
    # have no body, succeed. That looks exactly like a missing s3:PutObject
    # permission and is not one.
    sig_headers = {"host": HOST, "x-amz-content-sha256": phash, "x-amz-date": amzdate}
    if body:
        sig_headers["content-type"] = CONTENT_TYPE
    names = sorted(sig_headers)
    headers = "".join("%s:%s\n" % (n, sig_headers[n]) for n in names)
    signed = ";".join(names)
    creq = "\n".join([method, uri, query, headers, signed, phash])
    scope = "%s/%s/%s/aws4_request" % (datestamp, REGION, SERVICE)
    sts = "\n".join(["AWS4-HMAC-SHA256", amzdate, scope,
                     hashlib.sha256(creq.encode()).hexdigest()])
    k = _sign(_sign(_sign(_sign(("AWS4" + CFG["s3_key"]).encode(), datestamp),
                          REGION), SERVICE), "aws4_request")
    sig = hmac.new(k, sts.encode(), hashlib.sha256).hexdigest()
    auth = ("AWS4-HMAC-SHA256 Credential=%s/%s, SignedHeaders=%s, Signature=%s"
            % (CFG["s3_access"], scope, signed, sig))

    url = "%s://%s%s" % (SCHEME, HOST, uri) + (("?" + query) if query else "")
    hdrs = {"Host": HOST, "x-amz-date": amzdate, "x-amz-content-sha256": phash,
            "Authorization": auth}
    if body:
        hdrs["Content-Type"] = CONTENT_TYPE
        hdrs["Content-Length"] = str(len(body))
    req = urllib.request.Request(url, data=body if body else None, method=method,
                                 headers=hdrs)
    with urllib.request.urlopen(req, timeout=600) as r:
        return r.status, r.read()


def ls(prefix=""):
    import re
    keys, token = [], None
    while True:
        q = "list-type=2&max-keys=1000&prefix=" + urllib.parse.quote(prefix, safe="")
        if token:
            q += "&continuation-token=" + urllib.parse.quote(token, safe="")
        q = "&".join(sorted(q.split("&")))
        _, body = request("GET", "", q)
        txt = body.decode("utf-8", "replace")
        keys += re.findall(r"<Key>([^<]+)</Key>", txt)
        m = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", txt)
        if not m:
            return keys
        token = m.group(1)


def put(local, key):
    with open(local, "rb") as fh:
        body = fh.read()
    st, _ = request("PUT", key, "", body)
    return st, len(body)


def get(key, local):
    st, body = request("GET", key)
    d = os.path.dirname(os.path.abspath(local))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(local, "wb") as fh:
        fh.write(body)
    return st, len(body)


def _walk(root):
    for dp, _, fns in os.walk(root):
        for fn in fns:
            p = os.path.join(dp, fn)
            yield p, os.path.relpath(p, root).replace(os.sep, "/")


def main(argv):
    if len(argv) < 2:
        sys.exit(__doc__)
    cmd = argv[1]
    if cmd == "ls":
        for k in ls(argv[2] if len(argv) > 2 else ""):
            print(k)
    elif cmd == "put":
        st, n = put(argv[2], argv[3])
        print("PUT %s -> %s  HTTP %s  %.2f MB" % (argv[2], argv[3], st, n / 1e6))
    elif cmd == "get":
        st, n = get(argv[2], argv[3])
        print("GET %s -> %s  HTTP %s  %.2f MB" % (argv[2], argv[3], st, n / 1e6))
    elif cmd == "putdir":
        total = 0
        for p, rel in sorted(_walk(argv[2])):
            st, n = put(p, argv[3].rstrip("/") + "/" + rel)
            total += n
            print("  PUT %-55s HTTP %s  %7.2f MB" % (rel, st, n / 1e6))
        print("uploaded %.2f MB" % (total / 1e6))
    elif cmd == "getdir":
        pre = argv[2].rstrip("/") + "/"
        for k in ls(pre):
            st, n = get(k, os.path.join(argv[3], k[len(pre):]))
            print("  GET %-55s HTTP %s  %7.2f MB" % (k[len(pre):], st, n / 1e6))
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main(sys.argv)
