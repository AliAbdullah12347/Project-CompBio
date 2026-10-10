#!/bin/bash
# Watchdog for the scheduler itself.
#
# wake.sh's zero-throughput check only runs WHEN a session runs, so it cannot
# detect the scheduler failing to fire at all. That happened on 2026-10-09:
# ProcessType=Background plus StartInterval meant the 23:56Z firing never
# occurred and the run sat idle for 4 hours with no signal anywhere.
#
# Runs hourly. If no wake has been recorded for STALE_H hours and nothing is
# currently running, it kickstarts the main agent and records the fact.

export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:$HOME/.local/bin"
ROOT="$HOME/Project-CompBio"; cd "$ROOT" || exit 1
LOG="network/RUNLOG/watchdog.log"
CRON="network/RUNLOG/cron.log"
AGENT="com.aaylab.networkarm"
STALE_H=6
mkdir -p network/RUNLOG
stamp() { date -u +%Y-%m-%dT%H:%M:%SZ; }

# Something is running: healthy, nothing to do.
if pgrep -f "claude -p" >/dev/null 2>&1 || pgrep -f "run_next.sh" >/dev/null 2>&1 \
   || pgrep -f "exec/R .*network/experiments" >/dev/null 2>&1; then
  exit 0
fi

last=$(grep -o '=== wake [0-9T:-]*Z ===' "$CRON" 2>/dev/null | tail -1 | grep -o '[0-9]\{4\}-[0-9-]*T[0-9:]*Z')
[ -z "$last" ] && { echo "$(stamp) no wake ever recorded; kickstarting" >> "$LOG"; launchctl kickstart "gui/$(id -u)/$AGENT" 2>>"$LOG"; exit 0; }

last_s=$(date -j -u -f %Y-%m-%dT%H:%M:%SZ "$last" +%s 2>/dev/null)
[ -z "$last_s" ] && exit 0
age_h=$(( ( $(date +%s) - last_s ) / 3600 ))

if [ "$age_h" -ge "$STALE_H" ]; then
  echo "$(stamp) last wake $last was ${age_h}h ago (>= ${STALE_H}h) and nothing is running; kickstarting" >> "$LOG"
  launchctl kickstart "gui/$(id -u)/$AGENT" 2>>"$LOG"
  {
    echo "## $(stamp)  SCHEDULER STALL DETECTED AND RECOVERED"
    echo "No wake recorded for ${age_h}h (last: $last) and no session, driver or experiment running."
    echo "The watchdog kickstarted \`$AGENT\`. If this entry repeats, the scheduler is not"
    echo "firing on its own and the launchd configuration needs a human."
    echo
  } >> network/BLOCKED.md
  git add -- network/BLOCKED.md network/RUNLOG/watchdog.log 2>/dev/null
  git -c user.name="network-arm" -c user.email="aabdullah@colgate.edu" \
      commit -q -m "network: watchdog recovered a scheduler stall (${age_h}h idle)" -- network/ 2>/dev/null
  git push -q origin main 2>/dev/null
fi
