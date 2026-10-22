#!/bin/bash
# Run queue rows back-to-back instead of one per wake.
#
#   nohup bash network/chain.sh [budget_hours] > /tmp/chain.log 2>&1 &
#
# run_next.sh does exactly one row and exits, so between a row finishing and the
# next scheduled wake (every 5 h) the machine sits idle. This closes that gap.
#
# It is deliberately conservative:
#   - calls run_next.sh unmodified, so the driver lock, caps, checkpointing,
#     progress refunds, secret scan and commit/push all behave identically
#   - stops after BUDGET hours (default 12)
#   - stops after 3 consecutive non-productive iterations, so a row that fails
#     instantly cannot spin
#   - stops immediately if the queue is empty (exit 3) or a script is missing (exit 2)
#   - a scheduled wake arriving mid-chain simply finds the lock held and skips,
#     which is the existing designed behaviour
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

BUDGET_H=${1:-12}
END=$(( $(date +%s) + BUDGET_H * 3600 ))
LOG="network/RUNLOG/chain.log"
stamp() { date -u +%Y-%m-%dT%H:%M:%SZ; }
say() { echo "[$(stamp)] $*" | tee -a "$LOG"; }

say "chain starting, budget ${BUDGET_H}h"
idle=0; n=0
while [ "$(date +%s)" -lt "$END" ]; do
  done_before=$(awk -F'\t' 'NR>1 && $3=="DONE"{c++} END{print c+0}' network/RUNLOG/queue.tsv)
  n=$((n+1))
  say "iteration $n (DONE=$done_before)"

  out=$(bash network/run_next.sh 2>&1); rc=$?
  printf '%s\n' "$out" >> "$LOG"

  case "$rc" in
    3) say "queue has no claimable rows (rc=3); stopping"; break ;;
    2) say "next row needs a script that does not exist (rc=2); stopping for a human"; break ;;
  esac

  # A held lock is NOT a failed row - it means a scheduled session or an earlier
  # run_next is mid-experiment. Wait it out without touching the idle counter,
  # or the chain retires itself in under a second (which it did on first run).
  if printf '%s' "$out" | grep -q "another driver is running\|lock race"; then
    n=$((n-1))
    say "another driver holds the lock; waiting 120s"
    sleep 120
    continue
  fi

  done_after=$(awk -F'\t' 'NR>1 && $3=="DONE"{c++} END{print c+0}' network/RUNLOG/queue.tsv)
  if [ "$done_after" -gt "$done_before" ]; then
    say "progress: DONE $done_before -> $done_after"; idle=0
  else
    idle=$((idle+1)); say "no new DONE (rc=$rc); idle streak $idle/3"
    [ "$idle" -ge 3 ] && { say "3 iterations without a completion; stopping"; break; }
  fi
done
say "chain finished after $n iterations"
