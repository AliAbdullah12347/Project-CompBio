# Layout conventions — written to avoid merge conflicts

This repository is worked on from more than one machine (the iMac runs unattended;
a laptop pulls). The rules below exist so that two sides can work at once and `git pull`
never stops on a conflict.

## The rule: new work goes in NEW files, never into shared ones

| what | where | conflict risk |
|:--|:--|:--|
| a new analysis | `network/analyses/<NNN>-<slug>/` — its own folder, nothing shared | **none** |
| a queue experiment | `network/runs/<id>/` — its own folder | **none** |
| a finding | `network/findings/<date>-<slug>.md` — one file, never edited after writing | **none** |
| an experiment script | `network/experiments/<id>.R` — one file per row | **none** |

Because every file above is written once by exactly one author, two machines can both add
work and `git pull --rebase` merges cleanly.

## Files that ARE shared — treat with care

| file | nature | if it conflicts |
|:--|:--|:--|
| `network/RUNLOG/queue.tsv` | rewritten by the driver constantly | **take the iMac's version** (`git checkout --theirs`). It is the live run state; a laptop copy is stale by definition. |
| `network/RANKING.md`, `network/ranking.csv` | generated | **regenerate**: `Rscript network/build_index.R` |
| `network/RESULTS.md`, `network/results_index.csv` | generated | **regenerate** (same command) |
| `network/INDEX.md`, `network/catalog.csv` | generated | **regenerate**: `Rscript network/build_catalog.R` |
| `network/METHODS.md` | append-only log | keep **both** sides — the entries are independent and dated. |
| `network/BLOCKED.md` | append-only log | keep **both** sides. |
| `network/STATE.md` | overwritten each session | take the iMac's version. |

`RANKING.md`, `ranking.csv`, `RESULTS.md`, `results_index.csv`, `INDEX.md` and
`catalog.csv` are **derived files**. Never hand-edit them; edit the per-finding files
in `network/findings/` and regenerate.

## Finding anything

Start at `network/INDEX.md`. To search:

```bash
bash network/find.sh darkturquoise        # catalog + text + result files
bash network/find.sh --type finding composition
bash network/find.sh --files lithium      # just the file list
grep ',"finding",' network/catalog.csv
awk -F'","' '$4=="RETRACTED"' network/catalog.csv
```
That way a conflict in any of them is resolved by running one command rather than by
merging prose.

## From the laptop

```bash
git pull --rebase            # should be clean for anything under analyses/, runs/, findings/
# if queue.tsv or a generated file conflicts:
git checkout --theirs network/RUNLOG/queue.tsv && git add network/RUNLOG/queue.tsv
Rscript network/build_index.R && git add network/RANKING.md network/ranking.csv \
        network/RESULTS.md network/results_index.csv
git rebase --continue
```

## Do not

- Do not edit a file under `network/runs/` or `network/analyses/` from the laptop while
  the iMac may be writing it. Add a new folder instead.
- Do not hand-edit the generated files.
