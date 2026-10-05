# The project index

`INDEX.db` is a SQLite database of every file and folder under
`Implementation/` — where it is, how big it is, when it changed, what produced
it, whether it is still current, and what it means.

It exists because `DIRECTORY.txt` is good to read and bad to trust: it is
hand-written, so it goes stale the moment a file is added and nothing checks it.
The index keeps the two halves apart and joins them:

- the **filesystem** is the source of truth for what exists — discovered by
  walking the tree, never typed;
- **`descriptions.tsv`** is the source of truth for what a file *means* —
  curated, because no amount of walking tells you that a `.rds` is three hours
  of deconvolution or that a folder is quarantined.

A file on disk that matches no description shows up as `UNDESCRIBED` rather than
being quietly omitted, so the curation cannot drift out of date in silence.

---

## Rebuild

```bash
python index/build_index.py
```

Takes under a second. Run it after adding, moving or deleting files. It rewrites
`INDEX.db` and `INDEX.csv` from scratch, so there is no stale state to clear.

```bash
python index/build_index.py --check
```

Same, but exits non-zero if anything is undescribed — suitable for a pre-commit
hook.

---

## Ask it things

```bash
python index/build_index.py --query "SELECT path, size_human FROM v_biggest"
```

Five views cover the questions actually asked:

| view | answers |
|:--|:--|
| `v_current_outputs` | What results are current, and which script made each? |
| `v_code` | What did we write? |
| `v_do_not_use` | What is superseded or archived? |
| `v_biggest` | What is eating the disk? |
| `v_undescribed` | What has fallen out of the curation? |

Useful queries beyond those:

```sql
-- Where does a result come from?
SELECT path, produced_by, description FROM files WHERE name LIKE '%permutation%';

-- What does one script produce?
SELECT path, size_human FROM files WHERE produced_by LIKE '%01_refit%';

-- Disk by folder, biggest first
SELECT path, size_human, n_files_recursive, purpose FROM folders
ORDER BY size_bytes DESC LIMIT 15;

-- Everything current in the live analysis
SELECT path, description FROM files
WHERE folder LIKE 'de_v2%' AND status='current' ORDER BY path;

-- Anything that must not be cited
SELECT path, description FROM v_do_not_use;

-- How much is regenerable vs irreplaceable?
SELECT kind, COUNT(*) n, SUM(size_bytes)/1048576 mb FROM files
GROUP BY kind ORDER BY mb DESC;
```

Any SQLite client works too — `INDEX.db` is an ordinary database file.

---

## Schema

**`files`** — one row per file. `path` is relative to `Implementation/`.

| column | meaning |
|:--|:--|
| `path`, `folder`, `name`, `ext` | location |
| `kind` | `RAW` source data never edited · `CODE` written by us · `OUTPUT` generated, delete and re-run · `CACHE` expensive intermediate, regenerable but slow · `DOC` for humans · `EXT` from outside |
| `status` | `current` · `superseded` · `archived` · `external` · `raw` |
| `produced_by` | the script that creates it |
| `description` | one line, for someone who has never seen the project |
| `size_bytes`, `size_human`, `modified` | from the filesystem |

**`folders`** — one row per folder, with `purpose`, `status`, and file counts and
sizes computed from the file table rather than typed.

**`meta`** — build time, totals, and the undescribed count.

---

## Adding a description

Append a line to `descriptions.tsv`: five **tab**-separated fields.

```
pattern <TAB> kind <TAB> status <TAB> produced_by <TAB> description
```

`pattern` is an fnmatch glob on the path with forward slashes; `foo/**` matches
a folder and everything beneath it. **The last matching rule wins**, so put
folder-wide rules first and let file-specific rules below override them.
