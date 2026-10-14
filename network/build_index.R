#!/usr/bin/env Rscript
# Regenerates RANKING.md, ranking.csv, RESULTS.md and results_index.csv from the
# per-finding files in network/findings/ and the queue.
#
# These four files are DERIVED. Never hand-edit them: edit a findings/ file and
# re-run this. That is what makes a merge conflict in any of them resolvable by
# `Rscript network/build_index.R` rather than by merging prose. See CONTRIBUTING.md.
setwd(Sys.getenv("ARM3_ROOT", "/Users/Aaylab/Project-CompBio"))
fs <- sort(list.files("network/findings", pattern="[.]md$", full.names=TRUE))
# Match only the numerator, never the "/25" denominator, and normalise the
# Unicode minus sign (U+2212) that markdown prose uses for negative deltas.
grab <- function(x, pat) {
  m <- regmatches(x, regexpr(pat, x))
  if (!length(m)) return(NA_real_)
  v <- gsub("\u2212", "-", m)                 # Unicode minus -> ASCII
  v <- regmatches(v, regexpr("[+-]?[0-9]+([.][0-9]+)?", v))
  if (!length(v)) NA_real_ else as.numeric(v)
}

F <- do.call(rbind, lapply(fs, function(f) {
  tx <- readLines(f, warn=FALSE); hd <- paste(tx, collapse="\n")
  title <- sub("^#\\s*", "", tx[grep("^# ", tx)[1]])
  ev <- grab(hd, "Evidence:\\*\\*\\s*[0-9]+"); iv <- grab(hd, "Interest:\\*\\*\\s*[0-9]+")
  dl <- grab(hd, "\u0394:\\*\\*\\s*[\u2212+-]?[0-9.]+")
  comp <- regmatches(hd, regexpr("E=[0-9] R=[0-9] B=[0-9] D=[0-9] C=[0-9]", hd))
  an <- regmatches(hd, regexpr("\\*\\*Analysis:\\*\\*[^\n]*", hd))
  dt <- regmatches(hd, regexpr("\\*\\*Date:\\*\\*\\s*[0-9-]+", hd))
  data.frame(file=basename(f), title=title,
             date=if(length(dt)) gsub("[^0-9-]","",dt) else "",
             analysis=if(length(an)) trimws(gsub("\\*\\*Analysis:\\*\\*|`","",an)) else "",
             evidence=ev, interest=iv, delta=dl,
             components=if(length(comp)) comp else "", stringsAsFactors=FALSE)
}))
miss <- is.na(F$delta)
if (any(miss)) F$delta[miss] <- round((F$interest[miss]/10 - F$evidence[miss]/25)*10, 1)

A <- F[order(-F$evidence, -F$interest), ]; A$rank <- seq_len(nrow(A))
B <- F[order(-F$interest, -F$evidence), ]; B$rank <- seq_len(nrow(B))

out <- c("# Arm 3 — Live Result Ranking",
"", "> **GENERATED FILE — do not hand-edit.** Rebuild with `Rscript network/build_index.R`.",
"> Edit the per-finding files in `network/findings/` instead. On a merge conflict here,",
"> regenerate rather than merge. See `CONTRIBUTING.md`.", "",
paste0("**Findings:** ", nrow(F), "   **Updated:** ", format(Sys.Date())), "",
"Two rankings. **Table A (Evidence, 0–25)** governs the write-up; it rewards evidential",
"strength so a well-powered null can score top. **Table B (Interest, 0–10)** records how",
"striking a result is. **Δ = (Interest/10 − Evidence/25) × 10**: positive means the story",
"is ahead of the proof.", "",
"## Table A — Evidence", "",
"| # | Evidence | Interest | Δ | finding | components | source |",
"|--:|--:|--:|--:|:--|:--|:--|")
out <- c(out, sprintf("| %d | **%d/25** | %d/10 | %+.1f | [%s](findings/%s) | %s | %s |",
  A$rank, A$evidence, A$interest, A$delta, A$title, A$file, A$components, A$analysis))
out <- c(out, "", "## Table B — Interest", "",
"| # | Interest | Evidence | Δ | finding |", "|--:|--:|--:|--:|:--|")
out <- c(out, sprintf("| %d | **%d/10** | %d/25 | %+.1f | [%s](findings/%s) |",
  B$rank, B$interest, B$evidence, B$delta, B$title, B$file))
writeLines(out, "network/RANKING.md")

write.csv(data.frame(rank_evidence=A$rank, evidence=A$evidence, interest=A$interest,
                     delta=A$delta, finding=A$title, components=A$components,
                     analysis=A$analysis, file=A$file, date=A$date),
          "network/ranking.csv", row.names=FALSE)

q <- tryCatch(read.delim("network/RUNLOG/queue.tsv", stringsAsFactors=FALSE), error=function(e) NULL)
if (!is.null(q)) {
  st <- table(factor(q$status, levels=c("DONE","PARTIAL","RUNNING","PENDING","BLOCKED")))
  ri <- data.frame(run_id=q$id, family=q$family, status=q$status, wall_sec=q$wall_sec, note=q$note)
  write.csv(ri, "network/results_index.csv", row.names=FALSE)
  cat(sprintf("queue: %s\n", paste(names(st), st, sep="=", collapse=" ")))
}
cat(sprintf("regenerated: RANKING.md (%d findings), ranking.csv, results_index.csv\n", nrow(F)))
