#!/bin/bash
# Check every PENDING row for faults that would otherwise burn its full cap
# (up to 3h) on each of 3 attempts before anyone noticed. Read-only.
#
#   bash network/preflight_queue.sh
#
# Catches: missing script, parse error (base-003-remerge sat unparseable at
# queue position 1), `d` not in scope at an assert_expr_mode() call (my own
# mechanical edit of 29 files), and missing external inputs.
cd "$(dirname "$0")/.." || exit 1
fail=0; warn=0; n=0
while IFS=$'\t' read -r id _; do
  n=$((n+1)); f="network/experiments/$id.R"
  if [ ! -f "$f" ]; then echo "  MISSING SCRIPT   $id"; fail=$((fail+1)); continue; fi
  R --vanilla --slave -e "invisible(parse('$f'))" >/dev/null 2>&1 \
    || { echo "  PARSE ERROR      $id"; fail=$((fail+1)); continue; }
  al=$(grep -n "^[[:space:]]*assert_expr_mode(" "$f" | head -1 | cut -d: -f1)
  if [ -n "$al" ]; then
    dl=$(grep -nE "^[[:space:]]*d[[:space:]]*<<?-[[:space:]]*load_project" "$f" | head -1 | cut -d: -f1)
    { [ -z "$dl" ] || [ "$dl" -gt "$al" ]; } && { echo "  d NOT IN SCOPE   $id"; fail=$((fail+1)); }
  fi
  grep -oE '"(de_final|de_v2|de_analysis|mediation|data)/[A-Za-z0-9_/.-]+"' "$f" 2>/dev/null \
    | tr -d '"' | sort -u | while read -r p; do
        [ -e "$p" ] || echo "  missing input    $id -> $p"
      done
  grep -q "DEPENDENCY MISSING" "$f" && warn=$((warn+1))
done < <(awk -F'\t' 'NR>1 && $3=="PENDING"{print $1"\t"$2}' network/RUNLOG/queue.tsv)
echo
echo "  checked $n PENDING rows; $fail hard fault(s)"
echo "  (rows printing 'missing input' that build it themselves, e.g. bmind-046, are fine)"
[ "$fail" -gt 0 ] && exit 1
echo "  OK: no row will fail for a reason detectable before it runs"
