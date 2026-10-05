#!/usr/bin/env bash
export PATH="/c/Users/hp/R/R-4.6.1/bin/x64:$PATH"
export R_LIBS_USER="C:/Users/hp/R/win-library/4.6"
cd "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
for job in "LI raw" "LI ilr" "BPD raw" "BPD ilr"; do
  set -- $job
  echo "=================== $1 $2  started $(date) ==================="
  Rscript scripts/02_bmind.R "$1" "$2" 2>&1 | grep -v "renamed to normLibSizes"
done
echo "=================== ALL DONE $(date) ==================="
