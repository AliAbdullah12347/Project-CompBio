#!/bin/bash
# Fail if any EXPERIMENT uses d$logtpm as expression input instead of expr_for()
# or expr_samples(). Bypassing them skips per-sample offset removal, which
# silently scores a raw test matrix against an offset-removed reference. That
# happened in 36 places across 32 scripts and invalidated base-001's lith draws
# and the whole of base-002-ceil (whose cache was even written under an
# "_offsetrm" key). Benign uses (colnames/rownames/dim) are fine.
#
# Historical analyses/ are reported as warnings: they already ran, on the raw
# matrix, and are documented in findings/2026-10-15-library-offset-artifact.md.
cd "$(dirname "$0")/.." || exit 1
scan() {                       # $1=glob label, $2...=files
  local label="$1"; shift; local n=0
  for f in "$@"; do
    [ -f "$f" ] || continue
    while IFS= read -r raw; do
      local body="${raw#*:}"
      case "$(printf '%s' "$body" | sed 's/^[[:space:]]*//' | cut -c1)" in '#') continue ;; esac
      printf '%s' "$body" | grep -qE 'colnames\(|rownames\(|ncol\(|nrow\(|dim\(' && continue
      echo "  $label $f:${raw%%:*}  $(printf '%s' "$body" | sed 's/^[[:space:]]*//' | cut -c1-72)"
      n=$((n+1))
    done < <(grep -nE '\bd\$logtpm\b' "$f" 2>/dev/null)
  done
  echo "$n"
}
echo "== experiments (must be clean) =="
e=$(scan FAIL network/experiments/*.R); en=$(echo "$e" | tail -1); echo "$e" | sed '$d'
echo "== analyses (historical, warning only) =="
a=$(scan warn network/analyses/*/run.R); an=$(echo "$a" | tail -1); echo "$a" | sed '$d'
echo
echo "experiments: ${en:-0} violation(s)   analyses: ${an:-0} (already documented)"
[ "${en:-0}" -gt 0 ] && { echo "FAILED - use expr_for() or expr_samples()"; exit 1; }
echo "OK: no experiment uses d\$logtpm as expression input"
