#!/bin/bash
# Queue driver for the unattended network-analysis arm (Arm 3).
#
# CONTRACT WITH THE SESSION
#   The session writes network/experiments/<id>.R, THEN calls this script.
#   Use `run_next.sh --peek` to see which row is next without claiming it.
#   If this script claims a row whose script is missing it RELEASES the row
#   unchanged and exits 2 - it never retires work the session was about to do.
#
# NO STATUS IS ABSORBING EXCEPT BLOCKED, and BLOCKED is reachable only after
# MAX_ATTEMPTS genuine execution failures. A missing script, a cap hit and a
# crash are all recoverable; the queue can never drain itself to terminal.
#
# HARD CONSTRAINT: only paths under network/ are ever staged or committed.

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT" || exit 1
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:$HOME/.local/bin"

QUEUE="network/RUNLOG/queue.tsv"
PROGRESS="network/RUNLOG/progress.jsonl"
LOCKDIR="network/RUNLOG/.driver.lock"      # gitignored
STALE_SECS=21600          # 6 h before a RUNNING row is presumed dead
LOCK_MAX_AGE=28800        # 8 h before a lock is presumed dead regardless of PID
MAX_ATTEMPTS=3
DEFAULT_CAP=3600

ts()  { date -u +%Y-%m-%dT%H:%M:%SZ; }
now() { date +%s; }

# JSON-safe: escape backslash, quote, and strip control characters. A malformed
# line here corrupts the one file the next session reads to find its place.
beat() {
  local msg
  msg="$(printf '%s' "$4" | sed -e 's/\\/\\\\/g' -e 's/"/\\"/g' | tr -d '\000-\037')"
  printf '{"ts":"%s","run_id":"%s","phase":"%s","pct":%s,"msg":"%s"}\n' \
    "$(ts)" "$1" "$2" "$3" "$msg" >> "$PROGRESS"
}

# ---------------------------------------------------------------- queue I/O --
# Every rewrite goes through this so a crash mid-write can never truncate the
# queue: write to .tmp, verify it is non-empty and has the header, then mv.
rewrite() {
  local tmp="$QUEUE.tmp.$$"
  cat > "$tmp"
  # A short write can still carry the header, so compare row counts: the queue may
  # only ever grow (rows are appended), never shrink.
  local oldn newn
  oldn=$(wc -l < "$QUEUE" 2>/dev/null || echo 0); newn=$(wc -l < "$tmp" 2>/dev/null || echo 0)
  if [ ! -s "$tmp" ] || ! head -1 "$tmp" | grep -q '^id' || [ "$newn" -lt "$oldn" ]; then
    echo "FATAL: refusing to install a malformed queue" >&2
    rm -f "$tmp"; return 1
  fi
  mv "$tmp" "$QUEUE"
}

# PARTIAL and stale RUNNING both return to PENDING while attempts remain.
requeue() {
  local cutoff=$(( $(now) - STALE_SECS ))
  awk -F'\t' -v OFS='\t' -v cutoff="$cutoff" -v maxa="$MAX_ATTEMPTS" '
    NR==1 { print; next }
    $3=="RUNNING" {
      # -u is required: $5 is UTC and BSD `date -j` would otherwise read it as
      # local time, silently stretching the staleness window by the UTC offset.
      started = ""
      cmd = "date -j -u -f %Y-%m-%dT%H:%M:%SZ \"" $5 "\" +%s 2>/dev/null"
      cmd | getline started; close(cmd)
      # An unparseable stamp must be treated as stale, not skipped, or the row
      # is stuck RUNNING forever. getline leaves `started` untouched on failure,
      # which is why it is reset to "" above.
      if (started == "" || started+0 < cutoff) {
        # Judged on attempts alone a long, correctly-progressing experiment that is
        # killed mid-session three times would be retired. The driver refunds attempts
        # when a run produced real checkpoints, so by the time $4 reaches maxa here the
        # row has genuinely failed that many times without producing anything.
        if ($4+0 >= maxa) { $3="BLOCKED"; $8="dead RUNNING after " $4 " attempts, no progress" }
        else              { $3="PENDING"; $8="requeued from stale RUNNING" }
      }
    }
    $3=="PARTIAL" && $4+0 < maxa { $3="PENDING"; $8="requeued from PARTIAL (attempt " $4 ")" }
    { print }
  ' "$QUEUE" | rewrite
}

next_pending() { awk -F'\t' 'NR>1 && $3=="PENDING" { print $1; exit }' "$QUEUE"; }
field() { awk -F'\t' -v id="$1" -v n="$2" '$1==id { print $n; exit }' "$QUEUE"; }

# attempts is incremented HERE AND NOWHERE ELSE.
claim() {
  local id="$1"
  awk -F'\t' -v OFS='\t' -v id="$id" -v t="$(ts)" '
    NR==1 { print; next }
    $1==id && $3=="PENDING" { $3="RUNNING"; $5=t; $4=$4+1 }
    { print }
  ' "$QUEUE" | rewrite || return 1
  [ "$(field "$id" 3)" = "RUNNING" ] || return 1   # verify, do not assume
}

release() {   # undo a claim without burning an attempt
  awk -F'\t' -v OFS='\t' -v id="$1" -v n="$2" '
    NR==1 { print; next }
    $1==id { $3="PENDING"; $5=""; $4=($4+0>0 ? $4-1 : 0); $8=n }
    { print }
  ' "$QUEUE" | rewrite
}

mark() {      # mark <id> <status> <wall> <note>
  awk -F'\t' -v OFS='\t' -v id="$1" -v st="$2" -v w="$3" -v n="$4" -v t="$(ts)" '
    NR==1 { print; next }
    $1==id { $3=st; $6=t; $7=w; $8=n }
    { print }
  ' "$QUEUE" | rewrite
}

# ------------------------------------------------------------------- git -----
sync_up() {
  local msg="$1"

  # Commit by PATHSPEC rather than staging and unstaging. `git commit -- network/`
  # ignores the index entirely and commits only those paths, which removes three
  # separate failure modes at once:
  #   - rename detection could hide a deletion outside network/
  #   - `git restore --staged` left the outside file as a DIRTY WORKTREE change,
  #     which then made `git pull --rebase` refuse and stopped all pushes
  #   - anything a session happened to stage could ride along
  git add -- network/ 2>/dev/null

  # Public repo: refuse live credentials. Check the files that would be committed.
  local SECRET_RE="sk-ant-[a-z0-9]+-[A-Za-z0-9_-]{30,}|ghp_[A-Za-z0-9]{36}|github_pat_[A-Za-z0-9_]{50,}|AKIA[0-9A-Z]{16}|BEGIN [A-Z ]*PRIVATE KEY"
  local dirty=""
  while IFS= read -r f; do
    [ -n "$f" ] && [ -f "$f" ] || continue
    grep -qE "$SECRET_RE" "$f" 2>/dev/null && dirty="$dirty $f"
  done < <(git diff --cached --name-only --no-renames -- network/ 2>/dev/null; \
           git ls-files -m -- network/ 2>/dev/null)
  if [ -n "$dirty" ]; then
    echo "SECRET DETECTED in:$dirty"
    for f in $dirty; do git restore --staged -- "$f" 2>/dev/null; done
    { echo "## $(ts)  SECRET BLOCKED FROM A PUBLIC REPO"
      echo "Refused to commit these files; they match a live-credential pattern:"
      for f in $dirty; do echo "  - $f"; done
      echo "They remain on disk. Remove the secret before re-running."; echo
    } >> network/BLOCKED.md
  fi

  if git diff --quiet HEAD -- network/ 2>/dev/null && git diff --cached --quiet -- network/ 2>/dev/null; then
    echo "nothing to commit"; return 0
  fi
  git commit -q -m "$msg" -- network/ || { echo "commit failed"; return 1; }

  # Verify, do not assume: nothing outside network/ may be in that commit.
  local leaked
  leaked="$(git show --name-only --format= HEAD | grep -v '^network/' | grep . || true)"
  if [ -n "$leaked" ]; then
    echo "FATAL: commit contains paths outside network/:"; echo "$leaked"
    { echo "## $(ts)  SCOPE VIOLATION IN COMMIT $(git rev-parse --short HEAD)"
      echo "$leaked"; echo "Run halted. A human must reconcile."; echo; } >> network/BLOCKED.md
    return 1
  fi

  local n=0
  while [ $n -lt 3 ]; do
    # --autostash: tracked files outside network/ may be dirty (a session or a human
    # left them so), and without it rebase refuses and every push stops for the run.
    local err rc
    err="$(git -c user.name="network-arm" -c user.email="aabdullah@colgate.edu" \
             pull --rebase --autostash origin main 2>&1)"; rc=$?
    if [ $rc -ne 0 ]; then
      echo "PULL --REBASE FAILED rc=$rc"
      git rebase --abort 2>/dev/null
      { echo "## $(ts)  PULL --REBASE FAILED (rc=$rc), ABORTED"
        echo "Local commits are intact on main and unpushed. A human must reconcile."
        echo '```'; printf '%s\n' "$err" | head -12; echo '```'; echo
      } >> network/BLOCKED.md
      return 1
    fi
    if git push -q origin main 2>/dev/null; then echo "pushed"; return 0; fi
    n=$((n+1)); echo "push failed ($n/3)"; sleep $((n*20))
  done
  echo "push failed 3x; work is committed locally"
  return 1
}

# =============================== main ========================================
if [ "${1:-}" = "--peek" ]; then
  # READ-ONLY by design. It must never call requeue(): that rewrites queue.tsv,
  # and --peek holds no lock, so a peek racing a running driver would silently
  # discard the driver's status update. The driver requeues when it runs.
  id="$(next_pending)"
  if [ -z "$id" ]; then
    echo "QUEUE_EMPTY_PENDING"
    awk -F'\t' 'NR>1{c[$3]++} END{for(s in c) printf "%s\t%d\n", s, c[s]}' "$QUEUE"
    echo "note: PARTIAL and stale RUNNING rows are returned to PENDING by the next run"
    exit 3
  fi
  printf 'id\t%s\nfamily\t%s\ncap_sec\t%s\nspec\t%s\nscript\tnetwork/experiments/%s.R\nexists\t%s\n' \
    "$id" "$(field "$id" 2)" "$(field "$id" 9)" "$(field "$id" 10)" "$id" \
    "$([ -f "network/experiments/$id.R" ] && echo yes || echo NO)"
  awk -F'\t' 'NR>1{c[$3]++} END{printf "counts"; for(s in c) printf "\t%s=%d", s, c[s]; printf "\n"}' "$QUEUE"
  exit 0
fi

if ! mkdir "$LOCKDIR" 2>/dev/null; then
  stale=0
  if [ -f "$LOCKDIR/pid" ]; then
    kill -0 "$(cat "$LOCKDIR/pid" 2>/dev/null)" 2>/dev/null || stale=1
    # PIDs are recycled, so liveness alone is not enough: an old lock whose PID
    # now belongs to an unrelated process would wedge the run forever.
    age=$(( $(now) - $(stat -f %m "$LOCKDIR" 2>/dev/null || echo 0) ))
    [ "$age" -gt "$LOCK_MAX_AGE" ] && stale=1
  else
    stale=1
  fi
  if [ $stale -eq 1 ]; then rm -rf "$LOCKDIR"; mkdir "$LOCKDIR" 2>/dev/null || { echo "lock race"; exit 0; }
  else echo "another driver is running"; exit 0; fi
fi
echo $$ > "$LOCKDIR/pid"
trap 'rm -rf "$LOCKDIR"' EXIT

[ -f "$QUEUE" ] || { echo "no queue at $QUEUE"; exit 1; }

# A wedged repository must stop the run, not be worked around. Committing onto a
# detached HEAD produces results that `git rebase --abort` would silently destroy.
if [ -d .git/rebase-merge ] || [ -d .git/rebase-apply ]; then
  echo "FATAL: repository is mid-rebase; refusing to run"
  { echo "## $(ts)  REPOSITORY MID-REBASE"
    echo "An earlier pull --rebase did not finish. Every commit made now would land on a"
    echo "detached HEAD and be lost by \`git rebase --abort\`. A human must reconcile:"
    echo "  git -C \"$ROOT\" rebase --abort   # or --continue, after inspecting"
    echo; } >> network/BLOCKED.md
  exit 1
fi
if ! git symbolic-ref -q HEAD >/dev/null 2>&1; then
  echo "FATAL: detached HEAD; refusing to run"
  { echo "## $(ts)  DETACHED HEAD"
    echo "Commits made here are unreachable from main. A human must reconcile before"
    echo "the run continues. Do NOT \`rebase --abort\` without first checking"
    echo "\`git log HEAD --not main\` for orphaned results."
    echo; } >> network/BLOCKED.md
  exit 1
fi

requeue

ID="$(next_pending)"
[ -z "$ID" ] && { echo "QUEUE_EMPTY - no PENDING rows"; exit 3; }

SCRIPT="network/experiments/${ID}.R"
# Check BEFORE claiming. A missing script is the session's cue to write one,
# never a reason to retire the row.
if [ ! -f "$SCRIPT" ]; then
  echo "NEEDS_SCRIPT $ID -> $SCRIPT"
  echo "spec: $(field "$ID" 10)"
  beat "$ID" "needs_script" 0 "no script at $SCRIPT"
  exit 2
fi

claim "$ID" || { echo "FATAL: could not claim $ID"; beat "$ID" "error" 0 "claim failed"; exit 1; }

FAMILY="$(field "$ID" 2)"
CAP="$(field "$ID" 9)"
case "$CAP" in ''|*[!0-9]*) CAP=$DEFAULT_CAP ;; esac
[ "$CAP" -lt 60 ] && CAP=$DEFAULT_CAP

OUTDIR="network/runs/${ID}"; mkdir -p "$OUTDIR"
echo "=== $ID (family=$FAMILY, cap=${CAP}s) at $(ts) ==="
beat "$ID" "run" 5 "starting $SCRIPT cap=${CAP}s"

# Snapshot the tree outside network/ so anything a library writes there can be
# detected and cleaned: WGCNA's modulePreservation does exactly this by default.
# Once network/runs/<id>/ files are tracked, ' M network/...' lines appear too;
# filtering only '?? network/' would report them as strays outside the scope.
BEFORE="$(git status --porcelain | grep -vE '^(\?\?|.[ MADRCU]) network/' | sort)"

# Count checkpoints before the run. A long experiment (base-001 is 100 separate
# preservation runs and cannot fit one cap) is expected to span several sessions,
# marked PARTIAL each time. Burning an attempt on a run that MADE PROGRESS would
# retire it at MAX_ATTEMPTS while it was still working, so progress is measured
# rather than assumed.
CKPT_MARK="network/RUNLOG/.ckpt_mark_${ID}"; : > "$CKPT_MARK"   # mtime reference

START=$(now)
Rscript "$SCRIPT" "$ID" > "$OUTDIR/log.txt" 2>&1 &
RPID=$!
( # Poll in short sleeps. One long `sleep $CAP` survives `kill $WPID` as an
  # orphan holding this script's stdout, which can hang the caller. Output is
  # redirected so it can never hold that fd in the first place.
  while kill -0 $RPID 2>/dev/null; do
    [ $(( $(now) - START )) -ge "$CAP" ] || { sleep 5; continue; }
    echo "CAP EXCEEDED after ${CAP}s" >> "$OUTDIR/log.txt"
    pkill -TERM -P $RPID 2>/dev/null; kill -TERM $RPID 2>/dev/null
    sleep 10
    pkill -KILL -P $RPID 2>/dev/null; kill -KILL $RPID 2>/dev/null
    break
  done ) >/dev/null 2>&1 &
WPID=$!
wait $RPID; RC=$?
kill $WPID 2>/dev/null; wait $WPID 2>/dev/null
WALL=$(( $(now) - START ))

AFTER="$(git status --porcelain | grep -vE '^(\?\?|.[ MADRCU]) network/' | sort)"
if [ "$BEFORE" != "$AFTER" ]; then
  STRAY="$(comm -13 <(printf '%s\n' "$BEFORE") <(printf '%s\n' "$AFTER") | sed 's/^?? //')"
  echo "STRAY FILES OUTSIDE network/:"; printf '%s\n' "$STRAY"
  while IFS= read -r f; do
    [ -n "$f" ] || continue
    [ -f "$f" ] || continue
    if git check-ignore -q "$f" 2>/dev/null; then continue; fi   # already ignored: leave it
    mv "$f" "$OUTDIR/stray-$(basename "$f")" 2>/dev/null
  done <<< "$STRAY"
  { echo "## $(ts)  $ID wrote outside network/"; printf '%s\n' "$STRAY"
    echo "Moved into $OUTDIR/. Set explicit output paths (e.g. permutedStatisticsFile)."
    echo; } >> network/BLOCKED.md
fi

if   [ $RC -eq 0 ]; then STATUS=DONE;    NOTE="completed in ${WALL}s"
elif grep -q "CAP EXCEEDED" "$OUTDIR/log.txt" 2>/dev/null; then
     STATUS=PARTIAL; NOTE="hit ${CAP}s cap; partial results kept"
else STATUS=PARTIAL; NOTE="exit $RC after ${WALL}s; see log.txt"
fi
# Real progress = a file this run actually produced: newer than the marker, NON-EMPTY,
# and not a stray the driver itself quarantined here.
NEW_CKPT=$(find "$OUTDIR" -type f \( -name '*.rds' -o -name '*.csv' \) \
             ! -name 'stray-*' -size +0 -newer "$CKPT_MARK" 2>/dev/null | wc -l | tr -d ' ')
rm -f "$CKPT_MARK"
PROGRESS_MADE=0
[ "${NEW_CKPT:-0}" -gt 0 ] && PROGRESS_MADE=1

# Cap the refunds. An experiment that writes one junk .csv per run would otherwise
# refund forever and could never reach BLOCKED. 20 refunds x 3 h exceeds the whole run.
RF="network/RUNLOG/.refunds_${ID}"
RFN=$(cat "$RF" 2>/dev/null || echo 0)
if [ "$STATUS" = "PARTIAL" ] && [ "$PROGRESS_MADE" = "1" ] && [ "$RFN" -ge 20 ]; then
  PROGRESS_MADE=0
  NOTE="$NOTE; refund cap reached after $RFN refunds"
fi
if [ "$STATUS" = "PARTIAL" ] && [ "$PROGRESS_MADE" = "1" ]; then
  echo $((RFN + 1)) > "$RF"
  # Refund the attempt: this run advanced the experiment, it did not fail.
  awk -F'\t' -v OFS='\t' -v id="$ID" '
    NR==1 { print; next }
    $1==id { $4 = ($4+0 > 0 ? $4-1 : 0) }
    { print }
  ' "$QUEUE" | rewrite
  NOTE="$NOTE; +$NEW_CKPT new checkpoints (attempt refunded)"
  beat "$ID" "progress" 50 "$NEW_CKPT new non-empty checkpoints"
fi

ATT="$(field "$ID" 4)"
if [ "$STATUS" = "PARTIAL" ] && [ "$PROGRESS_MADE" = "0" ] && [ "${ATT:-0}" -ge "$MAX_ATTEMPTS" ]; then
  STATUS=BLOCKED; NOTE="$NOTE (no progress in $ATT attempts)"
fi

echo "=== $ID -> $STATUS (${WALL}s) ==="
beat "$ID" "done" 100 "$ID $STATUS ${WALL}s"
mark "$ID" "$STATUS" "$WALL" "$NOTE"
[ "$STATUS" != "DONE" ] && { echo "## $(ts)  $ID  $NOTE"; echo '```'; tail -20 "$OUTDIR/log.txt"; echo '```'; echo; } >> network/BLOCKED.md

sync_up "network: $ID $STATUS - $NOTE"
