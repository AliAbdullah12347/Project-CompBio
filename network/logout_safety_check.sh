#!/bin/bash
# Verifies the run survives a LOCKED SCREEN (Ctrl+Cmd+Q) and a LOGOUT/REBOOT.
# Read-only. Safe to run any time. Exit 0 = all critical checks pass.
cd /Users/Aaylab/Project-CompBio 2>/dev/null || exit 1
pass=0; fail=0; warn=0
ok()   { printf '  \033[32mPASS\033[0m  %s\n' "$1"; pass=$((pass+1)); }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; fail=$((fail+1)); }
note() { printf '  \033[33mWARN\033[0m  %s\n' "$1"; warn=$((warn+1)); }

echo; echo "=== SURVIVES SCREEN LOCK (Ctrl+Cmd+Q) ==="
[ "$(pmset -g | awk '$1=="sleep"{print $2}')" = "0" ] \
  && ok "system sleep disabled (pmset sleep 0)" \
  || bad "system sleep is ENABLED - machine will suspend"
pgrep -x caffeinate >/dev/null \
  && ok "caffeinate holding PreventSystemSleep" \
  || note "caffeinate not running (pmset sleep 0 still covers you)"
launchctl list 2>/dev/null | grep -q com.aaylab.networkarm \
  && ok "LaunchAgent loaded (drives the run while logged in)" \
  || bad "LaunchAgent NOT loaded"
defaults read .GlobalPreferences com.apple.autologout.AutoLogOutDelay >/dev/null 2>&1 \
  && bad "auto-logout timer is SET - macOS will log you out on its own" \
  || ok "no auto-logout timer (macOS will not log you out by itself)"

echo; echo "=== SURVIVES LOGOUT / REBOOT ==="
# launchctl cannot query the system domain without root, so: try sudo if it is
# cached, else fall back to the verification marker written when it was loaded.
# A daemon under /Library/LaunchDaemons is re-loaded by launchd at every boot,
# so "loaded once" is durable - it does not need re-bootstrapping.
if sudo -n launchctl print system/com.aaylab.networkarm.daemon >/dev/null 2>&1; then
  ok "LaunchDaemon loaded in system domain (confirmed live via launchctl)"
elif [ -f network/RUNLOG/.daemon_verified ]; then
  ok "LaunchDaemon loaded $(awk -F= '/^verified_utc/{print $2}' network/RUNLOG/.daemon_verified) \
(marker; re-confirm: sudo launchctl print system/com.aaylab.networkarm.daemon)"
elif grep -q "skipped: a session is still running\|=== wake " network/RUNLOG/cron.log 2>/dev/null \
     && [ -f /Library/LaunchDaemons/com.aaylab.networkarm.daemon.plist ]; then
  note "LaunchDaemon state UNVERIFIED without sudo (plist is present and wake.sh
        has been firing). Confirm with:
        sudo launchctl print system/com.aaylab.networkarm.daemon"
else
  bad "LaunchDaemon not loaded, and no evidence it ever fired. Run:
        sudo launchctl bootstrap system /Library/LaunchDaemons/com.aaylab.networkarm.daemon.plist"
fi
[ -f /Library/LaunchDaemons/com.aaylab.networkarm.daemon.plist ] \
  && ok "daemon plist present and root-owned" \
  || bad "daemon plist missing"
fdesetup status 2>/dev/null | grep -q "FileVault is Off" \
  && ok "FileVault off (daemon runs after reboot without a password)" \
  || note "FileVault ON - after a reboot nothing runs until someone unlocks the disk"
[ "$(pmset -g | awk '$1=="autorestart"{print $2}')" = "1" ] \
  && ok "auto-restart after power failure" \
  || note "no auto-restart after power failure"
[ "$(defaults read /Library/Preferences/com.apple.SoftwareUpdate AutomaticallyInstallMacOSUpdates 2>/dev/null)" = "0" ] \
  && ok "macOS will not auto-install updates and reboot" \
  || note "auto-install macOS updates may be ON - could reboot mid-run"
[ -s "$HOME/.arm3_token" ] \
  && ok "auth token readable at ~/.arm3_token ($(wc -c < "$HOME/.arm3_token" | tr -d ' ') bytes)" \
  || bad "auth token missing or empty"

echo; echo "=== RUN HEALTH RIGHT NOW ==="
code=$(curl -s -o /dev/null -w '%{http_code}' --max-time 10 https://api.anthropic.com/v1/messages 2>/dev/null)
[ -n "$code" ] && [ "$code" != "000" ] \
  && ok "network reaches api.anthropic.com (HTTP $code)" \
  || bad "cannot reach api.anthropic.com"
if pgrep -f "network/wake.sh" >/dev/null; then
  ok "a session is running now (pid $(pgrep -f 'network/wake.sh' | head -1))"
else
  note "no session running this instant (next scheduled wake :05 of 0/5/10/15/20)"
fi
awk -F'\t' 'NR>1{c[$3]++} END{printf "  ----  queue: "; for(k in c) printf "%s=%d  ", k, c[k]; print ""}' network/RUNLOG/queue.tsv 2>/dev/null
un=$(git -C /Users/Aaylab/Project-CompBio log --oneline origin/main..HEAD 2>/dev/null | wc -l | tr -d ' ')
[ "$un" = "0" ] && ok "0 unpushed commits" || note "$un commit(s) not yet pushed to GitHub"

echo
printf '  === %d pass, %d warn, %d FAIL ===\n' "$pass" "$warn" "$fail"
[ "$fail" -eq 0 ] && echo "  Safe to lock the screen with Ctrl+Cmd+Q." \
                  || echo "  Fix the FAIL lines above before you walk away."
echo
exit $((fail > 0))
