#!/bin/bash
# Launched by launchd every 5 h. Runs ONE headless Claude session that works the
# queue, then exits. Bounded, single-instance, and loud about auth failure.

export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:$HOME/.local/bin"
ROOT="$HOME/Project-CompBio"
cd "$ROOT" || exit 1

LOG="network/RUNLOG/cron.log"            # gitignored
SESSION_LOG="network/RUNLOG/last_session.out"   # gitignored
LOCK="network/RUNLOG/.wake.lock"         # gitignored
CLAUDE="$HOME/.local/bin/claude"
TOKEN_FILE="$HOME/.arm3_token"           # deliberately OUTSIDE this public repo
SESSION_CAP=14400                        # 4 h < the 5 h interval, so wakes never overlap
LOCK_MAX_AGE=21600                       # 6 h: a lock older than this is dead, whatever its PID says

mkdir -p network/RUNLOG
stamp() { date -u +%Y-%m-%dT%H:%M:%SZ; }

# --check runs every precondition and reports, but never commits, pushes, locks,
# or launches a session. Without this the only way to exercise the guards was to
# run the script for real - and its failure paths commit and push by design.
DRY=0
[ "${1:-}" = "--check" ] && DRY=1

# --- single instance, with an age fallback ----------------------------------
# PID liveness alone is not enough: PIDs are recycled, so after a reboot a stale
# lock whose number now belongs to an unrelated process would wedge every
# remaining wake-up for the rest of the run, logging only "skipped".
if [ $DRY -eq 0 ] && ! mkdir "$LOCK" 2>/dev/null; then
  stale=0
  [ -f "$LOCK/pid" ] || stale=1
  if [ -f "$LOCK/pid" ]; then
    kill -0 "$(cat "$LOCK/pid" 2>/dev/null)" 2>/dev/null || stale=1
    age=$(( $(date +%s) - $(stat -f %m "$LOCK" 2>/dev/null || echo 0) ))
    [ "$age" -gt "$LOCK_MAX_AGE" ] && stale=1
  fi
  if [ $stale -eq 1 ]; then
    echo "=== $(stamp) clearing stale lock ===" >> "$LOG"
    rm -rf "$LOCK"; mkdir "$LOCK" 2>/dev/null || { echo "=== $(stamp) lock race, skipping ===" >> "$LOG"; exit 0; }
  else
    echo "=== $(stamp) skipped: a session is still running ===" >> "$LOG"; exit 0
  fi
fi
if [ $DRY -eq 0 ]; then echo $$ > "$LOCK/pid"; trap 'rm -rf "$LOCK"' EXIT; fi

note_blocked() {
  # --check must be side-effect free: appending here would dirty a tracked file that
  # the next real run would then commit.
  [ $DRY -eq 1 ] && { echo "  [check] would record blocker: $1"; return 0; }
  { echo "## $(stamp)  $1"; shift; printf '%s\n' "$@"; echo; } >> network/BLOCKED.md
}

# Commit ONLY network/, never -A. Mirrors run_next.sh's guard so the auth-failure
# path cannot smuggle in unrelated staged work.
commit_network() {
  [ $DRY -eq 1 ] && { echo "  [check] would commit: $1"; return 0; }
  # Never operate on a wedged repository: committing onto a detached HEAD produces
  # work that `rebase --abort` silently destroys.
  if [ -d .git/rebase-merge ] || [ -d .git/rebase-apply ] || ! git symbolic-ref -q HEAD >/dev/null 2>&1; then
    echo "=== $(stamp) repo is mid-rebase or detached; refusing to commit ===" >> "$LOG"
    return 1
  fi
  git add -- network/ 2>/dev/null
  if git diff --quiet HEAD -- network/ 2>/dev/null && git diff --cached --quiet -- network/ 2>/dev/null; then
    return 0
  fi
  git commit -q -m "$1" -- network/ || return 1
  if ! git -c user.name="network-arm" -c user.email="aabdullah@colgate.edu" \
         pull --rebase --autostash -q origin main 2>/dev/null; then
    git rebase --abort 2>/dev/null
    echo "=== $(stamp) pull --rebase failed; aborted, commit kept locally ===" >> "$LOG"
    return 1
  fi
  git push -q origin main 2>/dev/null
}

# --- preconditions ----------------------------------------------------------
[ -x "$CLAUDE" ] || { echo "=== $(stamp) FATAL: no claude at $CLAUDE ===" >> "$LOG"; exit 1; }

# A truncated paste is a real failure mode: a 50-byte token passes -s and -r but
# fails every API call, which would look exactly like an expired token.
if [ ! -s "$TOKEN_FILE" ]; then
  echo "=== $(stamp) FATAL: no token at $TOKEN_FILE ===" >> "$LOG"
  note_blocked "NO AUTH TOKEN" "The auth token file is missing or empty." \
    "A human must run \`claude setup-token\` and store the FULL token in the configured file (chmod 600)."
  commit_network "network: blocked - no auth token"; exit 1
fi
TOKLEN=$(wc -c < "$TOKEN_FILE" | tr -d ' ')
if [ "$TOKLEN" -lt 80 ] || ! grep -q '^sk-ant-' "$TOKEN_FILE"; then
  echo "=== $(stamp) FATAL: token looks malformed ($TOKLEN bytes) ===" >> "$LOG"
  note_blocked "MALFORMED AUTH TOKEN" \
    "The auth token file is $TOKLEN bytes; a real token is ~108 characters." \
    "Most likely a paste truncated at a line break. Re-store the FULL token."
  commit_network "network: blocked - malformed auth token"; exit 1
fi
CLAUDE_CODE_OAUTH_TOKEN="$(cat "$TOKEN_FILE")"; export CLAUDE_CODE_OAUTH_TOKEN

if [ $DRY -eq 1 ]; then
  echo "  [check] all preconditions PASSED"
  echo "  [check] claude:  $CLAUDE"
  echo "  [check] token:   $TOKLEN bytes, prefix $(cut -c1-12 "$TOKEN_FILE")…"
  echo "  [check] cap:     ${SESSION_CAP}s"
  echo "  [check] would now launch a session; stopping here"
  exit 0
fi

{ echo; echo "============================================================"
  echo "=== wake $(stamp) ==="; echo "============================================================"
} >> "$LOG"

# --- run one bounded session ------------------------------------------------
# Without a cap, a single hung session holds the lock and every later wake logs
# only "skipped" - the run dies silently with four days still on the clock.
"$CLAUDE" -p "$(cat network/RESUME_PROMPT.md)" \
  --dangerously-skip-permissions > "$SESSION_LOG" 2>&1 < /dev/null &
CPID=$!
( sleep "$SESSION_CAP"
  if kill -0 $CPID 2>/dev/null; then
    echo "[wake.sh] SESSION CAP ${SESSION_CAP}s EXCEEDED - terminating" >> "$SESSION_LOG"
    pkill -TERM -P $CPID 2>/dev/null; kill -TERM $CPID 2>/dev/null
    sleep 15; pkill -KILL -P $CPID 2>/dev/null; kill -KILL $CPID 2>/dev/null
  fi ) &
WPID=$!
wait $CPID; RC=$?
kill $WPID 2>/dev/null; wait $WPID 2>/dev/null

cat "$SESSION_LOG" >> "$LOG"
echo "=== session exited rc=$RC at $(stamp) ===" >> "$LOG"

# --- auth failure must never be silent --------------------------------------
# Matched against the session's OWN last 40 lines only. Matching the whole log,
# or BLOCKED.md, would self-trigger on text this script itself wrote earlier.
if tail -40 "$SESSION_LOG" | grep -qE "OAuth access token has expired|Failed to authenticate|API Error: 401"; then
  note_blocked "AUTHENTICATION FAILURE" \
    "The headless CLI could not authenticate, so this wake-up produced nothing." \
    "Every later wake-up fails the same way until a human runs \`claude setup-token\`" \
    "and stores the full token in the configured file." \
    "" "The queue, the data and the scheduler are all intact." \
    "(Log lines are deliberately not copied here - this repository is public.)"
  echo "!!! AUTH FAILURE - see network/BLOCKED.md !!!" >> "$LOG"
  commit_network "network: AUTH FAILURE - headless token rejected, run is stalled"
  exit 1
fi

# --- zero-throughput watchdog ----------------------------------------------
# Every silent failure in this design presents identically: logs look normal and
# nothing completes. This is the one check that can tell the difference.
DONE_N=$(awk -F'\t' 'NR>1 && ($3=="DONE") {n++} END{print n+0}' network/RUNLOG/queue.tsv 2>/dev/null)
case "$DONE_N" in ''|*[!0-9]*)
  note_blocked "QUEUE UNREADABLE" "Could not count DONE rows in network/RUNLOG/queue.tsv." \
    "The zero-throughput watchdog cannot run until this is fixed."
  DONE_N="" ;;
esac
MARK="network/RUNLOG/.last_done_count"
PREV=$(cat "$MARK" 2>/dev/null || echo "")
STALL="network/RUNLOG/.stall_count"
if [ -n "$PREV" ] && [ "$DONE_N" = "$PREV" ]; then
  S=$(( $(cat "$STALL" 2>/dev/null || echo 0) + 1 )); echo "$S" > "$STALL"
  if [ "$S" -ge 3 ]; then
    note_blocked "NO PROGRESS IN $S CONSECUTIVE SESSIONS" \
      "DONE count has been stuck at $DONE_N across $S wake-ups (~$((S*5)) h)." \
      "The scheduler and sessions are running, but nothing is completing." \
      "Check the newest entries above and \`network/RUNLOG/progress.jsonl\`."
    commit_network "network: no progress in $S sessions - DONE stuck at $DONE_N"
  fi
else
  echo 0 > "$STALL"
fi
echo "$DONE_N" > "$MARK"
