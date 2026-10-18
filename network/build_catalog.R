#!/usr/bin/env Rscript
# Regenerates network/INDEX.md and network/catalog.csv — the single entry point
# for finding anything in this arm. GENERATED: never hand-edit either output.
# On a merge conflict in them, regenerate rather than merge (see CONTRIBUTING.md).
#
#   Rscript network/build_catalog.R
#
# catalog.csv is deliberately flat and one-row-per-thing so it greps well:
#   grep ',finding,'    network/catalog.csv
#   grep -i darkturquoise network/catalog.csv
#   awk -F, '$4=="RETRACTED"' network/catalog.csv
options(stringsAsFactors = FALSE)
setwd(normalizePath(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))

esc <- function(x) { x[is.na(x)] <- ""; gsub('"', '""', x) }
rows <- list()
add  <- function(...) rows[[length(rows) + 1]] <<- data.frame(...)

## ---- findings --------------------------------------------------------------
ff <- sort(Sys.glob("network/findings/*.md"))
for (f in ff) {
  tx    <- readLines(f, warn = FALSE)
  title <- sub("^#\\s*", "", tx[grep("^#\\s", tx)[1]])
  ev <- it <- NA_real_
  # NOT gsub("\\D",""): that strips digits out of the "/25" too and yields 2425.
  m <- regmatches(tx, regexpr("\\*\\*Evidence:\\*\\*\\s*([0-9]+)/25", tx))
  if (length(m)) ev <- as.numeric(sub(".*?([0-9]+)/25.*", "\\1", m[1]))
  m <- regmatches(tx, regexpr("\\*\\*Interest:\\*\\*\\s*([0-9]+)/10", tx))
  if (length(m)) it <- as.numeric(sub(".*?([0-9]+)/10.*", "\\1", m[1]))
  dt  <- sub(".*?(\\d{4}-\\d{2}-\\d{2}).*", "\\1", basename(f))
  src <- ""
  s <- grep("^\\*\\*Analysis:\\*\\*", tx, value = TRUE)
  if (length(s)) src <- gsub("`|\\*\\*Analysis:\\*\\*\\s*", "", s[1])
  status <- if (any(grepl("RETRACTED", tx[1:12]))) "RETRACTED"
            else if (any(grepl("CORRECTION", tx[1:14]))) "CORRECTED" else "live"
  add(type = "finding", id = basename(f), title = title, status = status,
      evidence = ev, interest = it, date = dt, path = f, family = "", extra = src)
}

## ---- analyses --------------------------------------------------------------
for (d in sort(Sys.glob("network/analyses/*/"))) {
  id <- basename(sub("/$", "", d))
  n  <- length(list.files(d, pattern = "\\.(csv|rds|txt|json)$", recursive = TRUE))
  rd <- file.path(d, "README.md")
  title <- if (file.exists(rd)) sub("^#\\s*", "", readLines(rd, warn = FALSE)[1]) else id
  crashed <- any(grepl("Execution halted|^Error",
                 unlist(lapply(Sys.glob(file.path(d, "results", "*output*.txt")),
                               function(p) tail(readLines(p, warn = FALSE), 30)))))
  add(type = "analysis", id = id, title = title,
      status = if (crashed) "ERROR" else if (n > 0) "has results" else "no results",
      evidence = NA, interest = NA,
      date = format(file.info(d)$mtime, "%Y-%m-%d"), path = d,
      family = "", extra = paste0(n, " result files"))
}

## ---- queue experiments -----------------------------------------------------
q <- "network/RUNLOG/queue.tsv"
if (file.exists(q)) {
  Q <- read.delim(q, check.names = FALSE)
  for (i in seq_len(nrow(Q)))
    add(type = "experiment", id = Q$id[i], title = substr(Q$spec[i], 1, 120),
        status = Q$status[i], evidence = NA, interest = NA, date = "",
        path = paste0("network/experiments/", Q$id[i], ".R"),
        family = Q$family[i], extra = paste0("cap=", Q$cap_sec[i], "s attempts=", Q$attempts[i]))
}

## ---- notes and run archives ------------------------------------------------
for (f in sort(Sys.glob("network/notes/*.md"))) {
  tx <- readLines(f, warn = FALSE)
  add(type = "note", id = basename(f), title = sub("^#\\s*", "", tx[grep("^#\\s", tx)[1]]),
      status = "", evidence = NA, interest = NA,
      date = sub(".*?(\\d{4}-\\d{2}-\\d{2}).*", "\\1", basename(f)),
      path = f, family = "", extra = "")
}
for (d in sort(c(Sys.glob("network/runs/*/"), Sys.glob("network/runs/_raw_matrix_archive_*/")))) {
  id <- basename(sub("/$", "", d))
  add(type = "run", id = id, title = id,
      status = if (grepl("archive", id)) "ARCHIVED (raw matrix)" else "live",
      evidence = NA, interest = NA,
      date = format(file.info(d)$mtime, "%Y-%m-%d"), path = d, family = "",
      extra = paste0(length(list.files(d, recursive = TRUE)), " files"))
}

C <- do.call(rbind, rows)
C <- C[order(factor(C$type, levels = c("finding","analysis","experiment","run","note")),
             -ifelse(is.na(C$evidence), -1, C$evidence), C$id), ]
# Every column is written as character so EVERY field is quoted. Numeric columns
# would be written bare, which breaks any consumer that splits on '","' - that is
# exactly how find.sh lost the path column.
Cc <- C
for (nm in names(Cc)) Cc[[nm]] <- ifelse(is.na(Cc[[nm]]), "", as.character(Cc[[nm]]))
write.csv(Cc, "network/catalog.csv", row.names = FALSE, quote = TRUE)

## ---- INDEX.md --------------------------------------------------------------
o <- c("# Arm 3 — Master Index",
  "", "> **GENERATED FILE — do not hand-edit.** Rebuild with `Rscript network/build_catalog.R`.",
  "> On a merge conflict here, regenerate rather than merge. See `CONTRIBUTING.md`.",
  "", sprintf("**Updated:** %s · **Catalogued:** %d items", format(Sys.Date()), nrow(C)), "",
  "## How to search", "",
  "```bash",
  "bash network/find.sh darkturquoise      # search everything for a term",
  "grep ',finding,' network/catalog.csv    # list findings",
  "awk -F, '$4==\"RETRACTED\"' network/catalog.csv",
  "```", "",
  "| file | what it is |", "|:--|:--|",
  "| [INDEX.md](INDEX.md) | this file — entry point |",
  "| [catalog.csv](catalog.csv) | every item, one row each, grep-friendly |",
  "| [RANKING.md](RANKING.md) | findings by Evidence and by Interest |",
  "| [RESULTS.md](RESULTS.md) | narrative write-up |",
  "| [METHODS.md](METHODS.md) | decision log — why each choice was made |",
  "| [STATE.md](STATE.md) | what is running right now |",
  "| [CONTRIBUTING.md](CONTRIBUTING.md) | file conventions, merge-conflict rules |", "")

sec <- function(ttl, sub, cols, hdr) {
  d <- C[C$type == sub, , drop = FALSE]
  if (!nrow(d)) return(character(0))
  c(sprintf("## %s (%d)", ttl, nrow(d)), "",
    paste0("| ", paste(hdr, collapse = " | "), " |"),
    paste0("|", paste(rep(":--", length(hdr)), collapse = "|"), "|"),
    apply(d[, cols, drop = FALSE], 1, function(r)
      paste0("| ", paste(ifelse(is.na(r) | r == "", "—", r), collapse = " | "), " |")), "")
}
F <- C[C$type == "finding", ]
F$link <- sprintf("[%s](%s)", F$title, sub("^network/", "", F$path))
F$sc   <- ifelse(is.na(F$evidence), "—", sprintf("%s/25 · %s/10", F$evidence, F$interest))
o <- c(o, sprintf("## Findings (%d — %d live, %d retracted)", nrow(F),
                  sum(F$status != "RETRACTED"), sum(F$status == "RETRACTED")), "",
       "| Evidence · Interest | status | finding | date |", "|:--|:--|:--|:--|",
       apply(F[, c("sc","status","link","date")], 1,
             function(r) paste0("| ", paste(r, collapse = " | "), " |")), "")
o <- c(o, sec("Analyses", "analysis", c("id","status","extra","title"),
              c("id","status","output","title")))
o <- c(o, sec("Notes — infrastructure and defects", "note", c("date","id","title"),
              c("date","file","what happened")))
o <- c(o, sec("Run directories", "run", c("id","status","extra"), c("run","status","size")))
E <- C[C$type == "experiment", ]
o <- c(o, sprintf("## Queue experiments (%d)", nrow(E)), "",
       paste0("Status: ", paste(names(table(E$status)), table(E$status),
                                sep = "=", collapse = "  ")), "",
       "Full detail in `RUNLOG/queue.tsv`; one row per experiment in `catalog.csv`.", "",
       "| id | family | status | spec |", "|:--|:--|:--|:--|",
       apply(E[, c("id","family","status","title")], 1,
             function(r) paste0("| ", paste(r, collapse = " | "), " |")), "")
writeLines(o, "network/INDEX.md")
cat(sprintf("regenerated: INDEX.md + catalog.csv (%d items: %s)\n", nrow(C),
            paste(names(table(C$type)), table(C$type), sep="=", collapse=" ")))
