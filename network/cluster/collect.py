#!/usr/bin/env python3
"""Pull Arm 3 cluster results out of S3 and report what landed.

    python network/cluster/collect.py            # status only, downloads nothing
    python network/cluster/collect.py --fetch    # also unpack results into network/runs/
    python network/cluster/collect.py --put FILE KEY   # send an input the other way

A pod's own files and logs are deleted a minute after its job ends, so S3 is
the only record. This reconciles what is there against what the queue says was
asked for, and prints the headline number from every row that finished.
"""
import argparse, io, os, re, subprocess, sys, tarfile, csv
from s3ls import client, BUCKET, ROOT

RUNS = "network/runs"
QUEUE = "network/RUNLOG/queue.tsv"


def listing(s3, sub=""):
    out, tok = {}, None
    while True:
        kw = dict(Bucket=BUCKET, Prefix=ROOT + sub)
        if tok:
            kw["ContinuationToken"] = tok
        r = s3.list_objects_v2(**kw)
        for o in r.get("Contents", []):
            out[o["Key"][len(ROOT):]] = o["Size"]
        if not r.get("IsTruncated"):
            break
        tok = r["NextContinuationToken"]
    return out


def queue_rows():
    rows = {}
    if not os.path.exists(QUEUE):
        return rows
    with io.open(QUEUE, encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            rows[r["id"]] = r["status"]
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true", help="download and unpack")
    ap.add_argument("--only", help="space-separated run ids; fetch only these")
    ap.add_argument("--force", action="store_true",
                    help="allow a fetch to overwrite committed, unmodified result files")
    ap.add_argument("--put", nargs=2, metavar=("FILE", "KEY"),
                    help="upload an input under the prefix, e.g. data/bmind_lineage_expr.rds")
    a = ap.parse_args()
    s3 = client()
    if a.put:
        src, key = a.put
        size = os.path.getsize(src)
        print(f"uploading {src} ({size/1048576:.0f} MB) -> {ROOT}{key}")
        s3.upload_file(src, BUCKET, ROOT + key)
        print("  done")
        return
    keys = listing(s3)

    # A sliced row uploads runs/<id>.sliceNN.tar.gz, one per pod. They all
    # unpack into the same network/runs/<id>/ directory, which is what the
    # row's aggregation step reads, so the id is the part before ".slice".
    run_keys = [k for k in keys if k.startswith("runs/") and k.endswith(".tar.gz")]
    runs = sorted({re.sub(r"[.]slice\d+$", "",
                          k[len("runs/"):-len(".tar.gz")]) for k in run_keys})
    if a.only:
        want = set(a.only.split())
        unknown = want - set(runs)
        if unknown:
            sys.exit("--only: no uploads for " + ", ".join(sorted(unknown)))
        runs = [r for r in runs if r in want]
        run_keys = [k for k in run_keys
                    if re.sub(r"[.]slice\d+$", "",
                              k[len("runs/"):-len(".tar.gz")]) in want]
    draws = sorted(k for k in keys if k.startswith("perm031/perm_draw_"))
    logs = sorted(k for k in keys if k.startswith("logs/"))
    q = queue_rows()

    print(f"S3 under {ROOT}: {len(keys)} objects\n")
    print(f"  finished rows uploaded : {len(runs)}")
    print(f"  perm-031 draws uploaded: {len(draws)}  (target 200, 1-20 pre-existing)")
    print(f"  pod logs               : {len(logs)}")
    print(f"  reference cache        : {'present' if any('ref_control_all' in k for k in keys) else 'MISSING'}")

    pending = [i for i, s in q.items() if s == "PENDING"]
    done_here = [r for r in runs if r in q]
    todo = sorted(set(pending) - set(runs))
    print(f"\n  queue PENDING at last sync: {len(pending)}")
    print(f"  of those, results in S3   : {len(set(pending) & set(runs))}")
    if todo:
        print(f"  still outstanding         : {' '.join(todo)}")

    if draws:
        n = sorted(int(re.search(r"(\d+)", os.path.basename(k)).group(1)) for k in draws)
        miss = [i for i in range(21, 201) if i not in n]
        print(f"\n  perm-031 draws present: {len(n)} ({min(n)}-{max(n)})"
              f"{'; missing ' + str(len(miss)) if miss else '; slice complete'}")

    if not a.fetch:
        print("\n(status only — pass --fetch to download)")
        return

    os.makedirs(RUNS, exist_ok=True)

    # Which files under network/runs/ git is tracking, and are currently clean.
    # A sliced row uploads one tarball per pod and EVERY pod's tarball carries
    # that row's whole directory, including any aggregate the row wrote from
    # its own slice alone. extractall() let the last slice win, so on 14 Nov a
    # routine fetch silently rewrote cent-056's committed 20-split null with a
    # 1-split one -- K=20 became K=1, null_sd became NA and p went 0.048 ->
    # 0.5. It still looked like a result, which is the dangerous part.
    # So: never overwrite a committed, unmodified file with different content
    # unless asked. Untracked and already-dirty files are extracted normally.
    def tracked_clean():
        def git(*args):
            r = subprocess.run(["git"] + list(args), capture_output=True,
                               text=True, check=True)
            return {ln.strip() for ln in r.stdout.splitlines() if ln.strip()}
        try:
            tracked = git("ls-files", RUNS)
            return tracked - git("diff", "--name-only", "--", RUNS)
        except Exception:
            return set()          # not a checkout, or no git: guard is a no-op

    protected = set() if a.force else tracked_clean()

    def extract_guarded(tf, dest, label):
        skipped = []
        for m in tf.getmembers():
            out = os.path.normpath(os.path.join(dest, m.name)).replace(os.sep, "/")
            if m.isfile() and out in protected and os.path.exists(out):
                src = tf.extractfile(m)
                if src is not None and src.read() != open(out, "rb").read():
                    skipped.append(out)
                    continue
            tf.extract(m, dest)
        if skipped:
            print(f"    !! {label}: kept {len(skipped)} committed file(s) rather than "
                  f"overwrite with different content (--force to override):")
            for x in skipped[:6]:
                print(f"       {x}")

    for r in runs:
        mine = sorted(k for k in run_keys
                      if re.sub(r"[.]slice\d+$", "",
                                k[len("runs/"):-len(".tar.gz")]) == r)
        for key in mine:
            tmp = os.path.join(RUNS, "_" + os.path.basename(key))
            s3.download_file(BUCKET, ROOT + key, tmp)
            with tarfile.open(tmp) as tf:
                extract_guarded(tf, RUNS, os.path.basename(key))
            os.remove(tmp)
        print(f"  unpacked {r}" + (f" ({len(mine)} slices)" if len(mine) > 1 else ""))
    dd = os.path.join(RUNS, "perm-031")
    os.makedirs(dd, exist_ok=True)
    for k in draws:
        dst = os.path.join(dd, os.path.basename(k))
        if not os.path.exists(dst):
            s3.download_file(BUCKET, ROOT + k, dst)
    print(f"  {len(draws)} perm-031 draws in {dd}")

    # headline line from each row's log, so a glance says what happened
    print("\n== last line of each row's log ==")
    for r in runs:
        lg = os.path.join(RUNS, r, "log.txt")
        if os.path.exists(lg):
            tail = [l.rstrip() for l in io.open(lg, encoding="utf-8", errors="replace") if l.strip()]
            print(f"  {r:22s} {tail[-1][:100] if tail else '(empty)'}")


if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    main()
