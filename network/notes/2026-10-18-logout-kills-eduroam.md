# The real cause of the 3-day outage: logging out disconnects eduroam

**Supersedes the diagnosis in** `2026-10-18-three-days-lost-to-wifi.md`, which blamed
flaky campus Wi-Fi. The link was not flaky. It was off because nobody was logged in.

## Evidence

`last` shows the console sessions:

```
aaylab console  Sun Oct 18 16:29   still logged in
aaylab console  Thu Oct 15 16:14 - 17:44  (01:29)   <- logged OUT 17:44
aaylab console  Wed Oct 14 14:48 - 17:26  (02:38)   <- logged OUT 17:26
```

The dead window (2026-10-15 20:23Z through 2026-10-18 20:31Z) sits exactly inside the
gap between the Oct 15 logout and the Oct 18 login.

## Mechanism

The machine is on **eduroam**, which is 802.1X enterprise authentication. The stored
credential is a **user-scope** keychain item:

```
svce = "com.apple.network.eap.user.item.wlan.ssid.eduroam"
```

`eap.USER.item` means it lives in the login keychain. Logging out locks that keychain, so
macOS cannot complete 802.1X at the login window and the machine has **no network at
all**. This is why both symptoms appeared together and why neither looked like an
application fault:

- the LaunchDaemon fired on schedule (system domain — it genuinely does survive logout)
- each session launched and then died with `FailedToOpenSocket`
- `git pull --rebase` in the same session failed too

A LaunchDaemon surviving logout is necessary but not sufficient. **The network has to
survive logout as well**, and on user-scope 802.1X it does not.

## What to do instead of logging out

**Lock the screen** (Control-Command-Q). The session, the keychain, the Wi-Fi and every
running job stay up; a password is still required to get back in.

**To turn the display off while staying logged in:** `pmset displaysleepnow`, or a hot
corner set to "Put Display to Sleep". Display sleep is not logout and not lock — nothing
stops.

Note that **Control-Option-Q is not a macOS shortcut.** Lock Screen is
Control-Command-Q; Log Out is Shift-Command-Q. On a PC keyboard the key marked Alt is
Option and the key marked Win is Command, which makes the two easy to confuse.

## The permanent fix

Plug in Ethernet. `en0` is free, reports `status: inactive`, and Ethernet is already
service-order #1, so a cable takes over automatically with no configuration. Wired campus
ports do not generally need a user-keychain credential, which removes the dependency on
anyone being logged in.

## Detection

`network/logout_safety_check.sh` now fails loudly on this condition: it checks for an
active wired interface, looks for a user-scope 802.1X item, and reports the console user.
