#!/bin/bash
# Park the iMac safely: screen off and locked, session KEPT.
#
#   bash ~/Project-CompBio/network/park.sh
#
# Use this instead of logging out. Logging out locks the login keychain, which
# holds the eduroam 802.1X credential, so the machine loses its network at the
# login window and every scheduled session dies - that cost this run 7 days
# across two incidents. See network/notes/2026-10-18-logout-kills-eduroam.md.
#
# Lock Screen is Control-Command-Q. Log Out is Shift-Command-Q. They sit next to
# each other in the Apple menu and differ by one modifier; this script removes
# the chance of hitting the wrong one.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

echo "Pre-flight:"
run_alive=$(pgrep -f "network/run_next.sh" >/dev/null && echo yes || echo no)
echo "  experiment running : $run_alive"
echo "  console user       : $(stat -f%Su /dev/console)"
if curl -s -o /dev/null --max-time 8 https://api.anthropic.com/v1/messages; then
  echo "  network            : up"
else
  echo "  network            : DOWN - fix before walking away"
fi
un=$(git log --oneline origin/main..HEAD 2>/dev/null | wc -l | tr -d ' ')
echo "  unpushed commits   : $un"
[ "$un" != "0" ] && echo "  (they will push on the next successful session)"

echo
echo "Locking screen and sleeping the display. The session stays logged in,"
echo "the keychain stays unlocked, Wi-Fi stays up, jobs keep running."
echo

# Lock first, then blank. pmset alone turns the display off but does NOT lock.
if command -v open >/dev/null 2>&1; then
  open -a ScreenSaverEngine 2>/dev/null \
    || osascript -e 'tell application "System Events" to keystroke "q" using {control down, command down}' 2>/dev/null
fi
sleep 2
pmset displaysleepnow
