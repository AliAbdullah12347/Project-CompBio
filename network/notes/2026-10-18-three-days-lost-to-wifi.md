# Three days lost: every wake started a session into a dead network

**Window:** 2026-10-15 20:23Z to 2026-10-18 20:31Z. **Cost:** 21 of 24 sessions, ~105 h,
zero queue progress. DONE stayed at 1.

## What happened

Every wake fired on schedule and launched a session. Each died ~3 minutes in with
`API Error: Unable to connect to API (FailedToOpenSocket)`, and the `git pull --rebase`
in the same session failed too. Both failing together points at the link, not the API:
the default route is **en1**, the campus Wi-Fi (`IOSkywalkNetworkBSDClient`,
gateway 149.43.192.1).

Verified it was never the token or the environment — with the machine back online,
`claude -p` returns `PONG` both from a normal shell and from launchd's exact minimal
environment (`env -i PATH=... HOME=... USER=... LOGNAME=...`). The token is the same
108-byte file it has always been.

## Why nothing caught it

The throughput watchdog *did* catch it and wrote `NO PROGRESS IN 17/18/19/20/21
CONSECUTIVE SESSIONS` into `BLOCKED.md`. Nobody was there to read it.

The structural problem is that the hourly watchdog kickstarts only when **no wake has
been recorded for 6 h**. A wake that starts a session and fails still writes its
`=== wake ... ===` banner, so to the watchdog the scheduler looked perfectly healthy.
A failing run and a working run were indistinguishable.

## Fix

A connectivity gate in `wake.sh`, inserted **before** the wake banner is written:

- probes `api.anthropic.com` and `github.com`, retries for up to 30 min at 60 s intervals
- if the network never comes back it exits 0 **without launching a session and without
  writing the wake banner**, so no slot is consumed
- because no banner is written, the hourly watchdog now correctly sees "no wake" and
  kickstarts — the two mechanisms finally compose instead of masking each other
- `wake.sh --check` reports reachability too, so the health check surfaces it

## Still open

The real reliability fix is **wired Ethernet**. `en0` is a free Ethernet port; the
machine is running on campus Wi-Fi. A cable would remove this failure mode entirely
rather than retrying around it.
