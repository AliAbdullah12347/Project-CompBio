#!/usr/bin/env python3
"""List keys under network/arm3/<prefix> in the project bucket.

Credentials come from the namespace's own s3-credentials secret via kubectl, so
nothing is stored on disk. Used by the job-chaining watcher and for pulling
results down; pods write to S3 because their own files vanish a minute after
the job ends.
"""
import base64, json, os, subprocess, sys
import boto3, botocore

# Each shell invocation is fresh, so fall back to the kubeconfig sitting
# beside the repo rather than relying on the caller having exported it.
if not os.environ.get("KUBECONFIG"):
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    kc = os.path.join(here, "kubeconfig-user-aabdullah.yaml")
    if os.path.exists(kc):
        os.environ["KUBECONFIG"] = kc

def client():
    raw = subprocess.run(["kubectl", "get", "secret", "s3-credentials",
                          "-o", "jsonpath={.data}"],
                         capture_output=True, text=True).stdout
    d = {k: base64.b64decode(v).decode() for k, v in json.loads(raw).items()}
    return boto3.client("s3", endpoint_url=d["AWS_ENDPOINT_URL"],
                        aws_access_key_id=d["AWS_ACCESS_KEY_ID"],
                        aws_secret_access_key=d["AWS_SECRET_ACCESS_KEY"],
                        region_name=d.get("AWS_REGION", "us-east-1"),
                        config=botocore.config.Config(s3={"addressing_style": "path"}))

BUCKET = "aabdullah-lithium-bipolar-project-fall-2026"
ROOT = "network/arm3/"

if __name__ == "__main__":
    sub = sys.argv[1] if len(sys.argv) > 1 else ""
    s3 = client()
    tok, n = None, 0
    while True:
        kw = dict(Bucket=BUCKET, Prefix=ROOT + sub)
        if tok:
            kw["ContinuationToken"] = tok
        r = s3.list_objects_v2(**kw)
        for o in r.get("Contents", []):
            print(f"{o['Size']/1048576:9.2f} MB  {o['Key'][len(ROOT):]}")
            n += 1
        if not r.get("IsTruncated"):
            break
        tok = r["NextContinuationToken"]
    print(f"-- {n} objects", file=sys.stderr)
