#!/usr/bin/env python3
"""
build_index.py -- build a queryable index of everything under Implementation/.

WHY THIS EXISTS

DIRECTORY.txt is a hand-written description of the project. It is good to read
and bad to trust: it goes stale the moment a file is added, and nothing checks
it. This script keeps the two halves separate and joins them:

    the FILESYSTEM is the source of truth for what exists, how big it is and
    when it changed -- discovered by walking the tree, never typed by hand;

    descriptions.tsv is the source of truth for what a file MEANS -- curated,
    because no amount of walking tells you that a .rds is three hours of
    deconvolution or that a folder is quarantined.

A file on disk with no matching description shows up as UNDESCRIBED rather than
being silently omitted, so the curation can never drift out of date without
saying so.

OUTPUTS
    INDEX.db    SQLite, the queryable form
    INDEX.csv   the same file table, for spreadsheets and diffs

USAGE
    python index/build_index.py                 rebuild both outputs
    python index/build_index.py --check         rebuild, exit 1 if anything is
                                                undescribed (for a pre-commit hook)
    python index/build_index.py --query "SQL"   run a query against INDEX.db
    python index/build_index.py --undescribed   list files with no description
"""

import csv, fnmatch, os, sqlite3, sys, time
from datetime import datetime, timezone

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # Implementation/
DESC = os.path.join(ROOT, "index", "descriptions.tsv")
DB   = os.path.join(ROOT, "INDEX.db")
CSVF = os.path.join(ROOT, "INDEX.csv")
SKIP_DIRS = {".git", "__pycache__", ".ipynb_checkpoints"}


def human(n):
    for unit in ("B", "K", "M", "G", "T"):
        if n < 1024 or unit == "T":
            return f"{n:.0f}{unit}" if unit == "B" else f"{n:.1f}{unit}"
        n /= 1024


def load_descriptions():
    """Curated rules. Later-listed rules win, so a file-specific rule can
    override the folder-wide rule above it. Blank lines and # comments ignored."""
    rules = []
    with open(DESC, encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln.strip() or ln.lstrip().startswith("#"):
                continue
            parts = ln.split("\t")
            if len(parts) < 5:
                raise SystemExit(f"descriptions.tsv: need 5 tab-separated fields, got "
                                 f"{len(parts)} in: {ln[:80]}")
            pattern, kind, status, produced_by, description = parts[:5]
            rules.append(dict(pattern=pattern.strip(), kind=kind.strip(),
                              status=status.strip(), produced_by=produced_by.strip(),
                              description=description.strip()))
    return rules


def match(relpath, rules, is_dir=False):
    """Last matching rule wins. Patterns are fnmatch globs on the POSIX-style
    path relative to Implementation/; a trailing /** matches everything under
    a folder."""
    hit = None
    rp = relpath.replace(os.sep, "/")
    for r in rules:
        pat = r["pattern"]
        if pat.endswith("/**"):
            base = pat[:-3]
            if rp == base or rp.startswith(base + "/"):
                hit = r
        elif fnmatch.fnmatch(rp, pat):
            hit = r
    return hit


def walk(rules):
    files, folders = [], []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = sorted(d for d in dirnames if d not in SKIP_DIRS)
        rel_dir = os.path.relpath(dirpath, ROOT).replace(os.sep, "/")
        if rel_dir == ".":
            rel_dir = ""
        if rel_dir:
            m = match(rel_dir, rules, is_dir=True)
            folders.append(dict(
                path=rel_dir,
                purpose=(m or {}).get("description", "UNDESCRIBED"),
                status=(m or {}).get("status", "unknown")))
        for fn in sorted(filenames):
            full = os.path.join(dirpath, fn)
            rel = (f"{rel_dir}/{fn}" if rel_dir else fn)
            try:
                st = os.stat(full)
            except OSError:
                continue
            m = match(rel, rules)
            files.append(dict(
                path=rel, folder=rel_dir or ".", name=fn,
                ext=os.path.splitext(fn)[1].lstrip(".").lower(),
                kind=(m or {}).get("kind", "?"),
                status=(m or {}).get("status", "unknown"),
                produced_by=(m or {}).get("produced_by", ""),
                description=(m or {}).get("description", "UNDESCRIBED"),
                size_bytes=st.st_size, size_human=human(st.st_size),
                modified=datetime.fromtimestamp(st.st_mtime, timezone.utc)
                         .strftime("%Y-%m-%d %H:%M")))
    return files, folders


def build(files, folders):
    if os.path.exists(DB):
        os.remove(DB)
    con = sqlite3.connect(DB)
    c = con.cursor()
    c.executescript("""
    CREATE TABLE files (
      path TEXT PRIMARY KEY, folder TEXT, name TEXT, ext TEXT,
      kind TEXT, status TEXT, produced_by TEXT, description TEXT,
      size_bytes INTEGER, size_human TEXT, modified TEXT);
    CREATE TABLE folders (
      path TEXT PRIMARY KEY, purpose TEXT, status TEXT,
      n_files INTEGER, n_files_recursive INTEGER,
      size_bytes INTEGER, size_human TEXT);
    CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT);
    CREATE INDEX idx_files_folder ON files(folder);
    CREATE INDEX idx_files_kind   ON files(kind);
    CREATE INDEX idx_files_status ON files(status);
    -- Convenience views: the questions actually asked of this index.
    CREATE VIEW v_current_outputs AS
      SELECT path, size_human, produced_by, description FROM files
      WHERE kind='OUTPUT' AND status='current' ORDER BY path;
    CREATE VIEW v_code AS
      SELECT path, size_human, description FROM files
      WHERE kind='CODE' ORDER BY path;
    CREATE VIEW v_do_not_use AS
      SELECT path, status, description FROM files
      WHERE status IN ('superseded','archived') ORDER BY path;
    CREATE VIEW v_biggest AS
      SELECT path, size_human, kind, status, description FROM files
      ORDER BY size_bytes DESC LIMIT 25;
    CREATE VIEW v_undescribed AS
      SELECT path, size_human FROM files
      WHERE description='UNDESCRIBED' ORDER BY path;
    """)
    c.executemany("""INSERT INTO files
      (path,folder,name,ext,kind,status,produced_by,description,size_bytes,size_human,modified)
      VALUES (:path,:folder,:name,:ext,:kind,:status,:produced_by,:description,
              :size_bytes,:size_human,:modified)""", files)

    # Folder aggregates are computed from the file table, never typed.
    for f in folders:
        direct = [x for x in files if x["folder"] == f["path"]]
        rec = [x for x in files
               if x["path"] == f["path"] or x["path"].startswith(f["path"] + "/")]
        tot = sum(x["size_bytes"] for x in rec)
        c.execute("""INSERT INTO folders
          (path,purpose,status,n_files,n_files_recursive,size_bytes,size_human)
          VALUES (?,?,?,?,?,?,?)""",
                  (f["path"], f["purpose"], f["status"], len(direct), len(rec), tot, human(tot)))

    und = sum(1 for x in files if x["description"] == "UNDESCRIBED")
    for k, v in [("built_utc", datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")),
                 ("root", "Implementation/"),
                 ("n_files", str(len(files))), ("n_folders", str(len(folders))),
                 ("total_bytes", str(sum(x["size_bytes"] for x in files))),
                 ("total_human", human(sum(x["size_bytes"] for x in files))),
                 ("n_undescribed", str(und)),
                 ("builder", "index/build_index.py"),
                 ("descriptions", "index/descriptions.tsv")]:
        c.execute("INSERT INTO meta VALUES (?,?)", (k, v))
    con.commit()
    con.close()

    cols = ["path", "folder", "kind", "status", "produced_by", "size_human",
            "modified", "description"]
    with open(CSVF, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in sorted(files, key=lambda x: x["path"]):
            w.writerow(r)
    return und


def query(sql):
    con = sqlite3.connect(DB)
    con.row_factory = sqlite3.Row
    try:
        rows = con.execute(sql).fetchall()
    except sqlite3.Error as e:
        raise SystemExit(f"SQL error: {e}")
    if not rows:
        print("(no rows)")
        return
    hdr = rows[0].keys()
    w = [max(len(h), max(len(str(r[h])) for r in rows)) for h in hdr]
    w = [min(x, 78) for x in w]
    print("  ".join(h.ljust(x) for h, x in zip(hdr, w)))
    print("  ".join("-" * x for x in w))
    for r in rows:
        print("  ".join(str(r[h])[:x].ljust(x) for h, x in zip(hdr, w)))
    print(f"\n{len(rows)} row(s)")


def main():
    args = sys.argv[1:]
    if args and args[0] == "--query":
        return query(" ".join(args[1:]))
    if args and args[0] == "--undescribed":
        return query("SELECT * FROM v_undescribed")

    t0 = time.time()
    rules = load_descriptions()
    files, folders = walk(rules)
    und = build(files, folders)
    total = sum(f["size_bytes"] for f in files)
    print(f"indexed {len(files)} files in {len(folders)} folders, {human(total)}, "
          f"{time.time()-t0:.1f}s")
    print(f"  {DB}")
    print(f"  {CSVF}")
    by_kind = {}
    for f in files:
        by_kind[f["kind"]] = by_kind.get(f["kind"], 0) + 1
    print("  by kind: " + "  ".join(f"{k}={v}" for k, v in sorted(by_kind.items())))
    if und:
        print(f"\n  {und} file(s) UNDESCRIBED -- add a rule to index/descriptions.tsv")
        print("  list them with:  python index/build_index.py --undescribed")
        if "--check" in args:
            sys.exit(1)
    else:
        print("  every file has a description")


if __name__ == "__main__":
    main()
