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

    runs = sorted(k[len("runs/"):-len(".tar.gz")] for k in keys if k.startswith("runs/") and k.endswith(".tar.gz"))
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
    for r in runs:
        key = f"runs/{r}.tar.gz"
        tmp = os.path.join(RUNS, f"_{r}.tar.gz")
        s3.download_file(BUCKET, ROOT + key, tmp)
        with tarfile.open(tmp) as tf:
            tf.extractall(RUNS)
        os.remove(tmp)
        print(f"  unpacked {r}")
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
