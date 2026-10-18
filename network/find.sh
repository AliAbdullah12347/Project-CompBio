#!/bin/bash
# Search everything in the network arm: catalog entries, finding text, analysis
# code and results. Read-only.
#
#   bash network/find.sh darkturquoise
#   bash network/find.sh --type finding composition
#   bash network/find.sh --files lithium     # just the file list
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

TYPE=""; FILES_ONLY=0
while [[ "${1:-}" == --* ]]; do
  case "$1" in
    --type)  TYPE="$2"; shift 2 ;;
    --files) FILES_ONLY=1; shift ;;
    --help)  sed -n '2,9p' network/find.sh; exit 0 ;;
    *) echo "unknown flag: $1"; exit 2 ;;
  esac
done
Q="${1:-}"
[ -z "$Q" ] && { sed -n '2,9p' network/find.sh; exit 2; }

[ -f network/catalog.csv ] || { echo "no catalog yet — run: Rscript network/build_catalog.R"; exit 1; }

bold() { printf '\033[1m%s\033[0m\n' "$1"; }

if [ "$FILES_ONLY" -eq 0 ]; then
  bold "── catalog entries matching '$Q' ──"
  awk -F'","' -v q="$Q" -v t="$TYPE" '
    NR==1{next}
    {
      # R quotes every field, so strip the outer quotes off the first and last
      for (i = 1; i <= NF; i++) { gsub(/^"/, "", $i); gsub(/"$/, "", $i) }
      if (t != "" && $1 != t) next
      if (tolower($0) ~ tolower(q)) {
        printf "  %-11s %-32s %s\n", $1, substr($2,1,32), substr($3,1,74)
        printf "  %-11s   -> %s\n", "", $8
      }
    }' network/catalog.csv | head -60
  echo
fi

if [ "$FILES_ONLY" -eq 1 ]; then
  grep -rlniI "$Q" network/findings network/notes network/analyses \
       network/RESULTS.md network/METHODS.md 2>/dev/null | sort -u | sed 's|^|  |'
  exit 0
fi

bold "── text matches in findings / notes / analyses ──"
grep -rniI --color=never "$Q" \
     network/findings network/notes network/analyses network/RESULTS.md \
     network/METHODS.md 2>/dev/null \
  | grep -v "catalog.csv" \
  | awk -F: '{ printf "  %s:%s\n", $1, $2; $1="";$2=""; s=$0; sub(/^::/,"",s);
               gsub(/^[ \t]+/,"",s); print "      " substr(s,1,118) }' \
  | head -50

echo
bold "── result files whose contents mention '$Q' ──"
grep -rlI "$Q" network/analyses/*/results network/runs 2>/dev/null | head -20 | sed 's|^|  |'
echo
echo "tip: 'grep -i \"$Q\" network/catalog.csv' for the raw rows; INDEX.md for the overview"
