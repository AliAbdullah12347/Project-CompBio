# ==============================================================================
# doc_render.R -- a tiny document model rendered two ways from ONE source:
#
#   Markdown  parseable, diff-able, lives in the repo
#   HTML      plain semantic HTML that survives copy-paste into Google Docs:
#             open it in a browser, select all, copy, paste. Headings become
#             Google Docs headings (so the outline / navigation pane works),
#             tables become real tables, bold stays bold.
#
# Because both come from the same blocks, the two versions cannot disagree.
# Inline markup in text: **bold**, *italic*, `code`.
# ==============================================================================

H     <- function(level, text) list(type = "h", level = level, text = text)
Par   <- function(...) list(type = "p", text = paste0(...))
UL    <- function(items) list(type = "ul", items = items)
OL    <- function(items) list(type = "ol", items = items)
TB    <- function(df, caption = NULL) list(type = "tb", df = df, caption = caption)
NOTE  <- function(...) list(type = "note", text = paste0(...))
GENES <- function(title, genes, note = NULL) list(type = "genes", title = title, genes = genes, note = note)
CODE  <- function(text) list(type = "code", text = text)
HR    <- function() list(type = "hr")

## ---- inline conversion ------------------------------------------------------
esc_html <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;", s, fixed = TRUE)
  gsub(">", "&gt;", s, fixed = TRUE)
}
inline_html <- function(s) {
  s <- esc_html(s)
  s <- gsub("\\*\\*([^*]+)\\*\\*", "<b>\\1</b>", s, perl = TRUE)
  s <- gsub("`([^`]+)`", "<code>\\1</code>", s, perl = TRUE)
  s <- gsub("(?<![*\\w])\\*([^*\\n]+)\\*(?![*\\w])", "<i>\\1</i>", s, perl = TRUE)
  s
}
md_cell <- function(s) gsub("|", "\\|", s, fixed = TRUE)

fmt_df <- function(df) {
  for (j in seq_along(df)) df[[j]] <- ifelse(is.na(df[[j]]), "—", as.character(df[[j]]))
  df
}

## ---- Markdown ------------------------------------------------------------------
render_md <- function(blocks) {
  out <- character(0)
  for (b in blocks) {
    out <- c(out, switch(b$type,
      h  = c(paste(strrep("#", b$level), b$text), ""),
      p  = c(b$text, ""),
      ul = c(paste0("- ", b$items), ""),
      ol = c(paste0(seq_along(b$items), ". ", b$items), ""),
      note = c(paste0("> ", strsplit(b$text, "\n", fixed = TRUE)[[1]]), ""),
      code = c("```", b$text, "```", ""),
      hr = c("---", ""),
      genes = c(paste0("**", b$title, "** (", length(b$genes), ")"), "",
                if (!is.null(b$note)) c(paste0("*", b$note, "*"), "") else NULL,
                if (length(b$genes)) paste(b$genes, collapse = ", ") else "*none*", ""),
      tb = {
        df <- fmt_df(b$df)
        c(if (!is.null(b$caption)) c(paste0("*", b$caption, "*"), "") else NULL,
          paste0("| ", paste(md_cell(names(df)), collapse = " | "), " |"),
          paste0("|", paste(rep(":--", ncol(df)), collapse = "|"), "|"),
          apply(df, 1, function(r) paste0("| ", paste(md_cell(r), collapse = " | "), " |")),
          "")
      }))
  }
  out
}

## ---- HTML (Google-Docs-friendly) ----------------------------------------------
TD <- "border:1px solid #9a9a9a;padding:3px 7px;vertical-align:top;font-size:10pt"
TH <- paste0(TD, ";background:#ececec;font-weight:bold;text-align:left")
render_html <- function(blocks, title) {
  out <- c("<!DOCTYPE html>", "<html><head><meta charset=\"utf-8\">",
           paste0("<title>", esc_html(title), "</title>"),
           "<style>body{font-family:Arial,Helvetica,sans-serif;font-size:11pt;line-height:1.4;",
           "max-width:960px;margin:24px auto;padding:0 16px;color:#111;background:#fff}",
           "code{font-family:Consolas,'Courier New',monospace;font-size:10pt}",
           "h1{font-size:20pt}h2{font-size:16pt;margin-top:28px}h3{font-size:13pt}h4{font-size:11.5pt}</style>",
           "</head><body>")
  for (b in blocks) {
    out <- c(out, switch(b$type,
      h  = sprintf("<h%d>%s</h%d>", b$level, inline_html(b$text), b$level),
      p  = paste0("<p>", inline_html(b$text), "</p>"),
      ul = c("<ul>", paste0("<li>", inline_html(b$items), "</li>"), "</ul>"),
      ol = c("<ol>", paste0("<li>", inline_html(b$items), "</li>"), "</ol>"),
      note = paste0("<table style=\"border-collapse:collapse;width:100%\"><tr><td style=\"",
                    "border:1px solid #c08a00;background:#fff6dd;padding:6px 10px;font-size:10.5pt\">",
                    gsub("\n", "<br>", inline_html(b$text), fixed = TRUE), "</td></tr></table><p></p>"),
      code = paste0("<pre style=\"font-family:Consolas,'Courier New',monospace;font-size:10pt;",
                    "background:#f4f4f4;padding:6px\">", esc_html(b$text), "</pre>"),
      hr = "<hr>",
      genes = c(paste0("<p><b>", inline_html(b$title), "</b> (", length(b$genes), ")</p>"),
                if (!is.null(b$note)) paste0("<p><i>", inline_html(b$note), "</i></p>") else NULL,
                paste0("<p style=\"font-size:9.5pt\">",
                       if (length(b$genes)) esc_html(paste(b$genes, collapse = ", ")) else "<i>none</i>",
                       "</p>")),
      tb = {
        df <- fmt_df(b$df)
        c(if (!is.null(b$caption)) paste0("<p><i>", inline_html(b$caption), "</i></p>") else NULL,
          "<table style=\"border-collapse:collapse\">",
          paste0("<tr>", paste0("<th style=\"", TH, "\">", inline_html(names(df)), "</th>", collapse = ""), "</tr>"),
          apply(df, 1, function(r) paste0("<tr>", paste0("<td style=\"", TD, "\">", inline_html(r), "</td>",
                                                          collapse = ""), "</tr>")),
          "</table><p></p>")
      }))
  }
  c(out, "</body></html>")
}

## ---- number formatting helpers ------------------------------------------------
# Every formatter returns a real NA for a missing value, so fmt_df can show it as
# an em dash; formatC on its own would print the text "NA" into the table.
na_safe <- function(f) function(x, ...) { out <- f(x, ...); out[is.na(x)] <- NA_character_; out }
f2 <- na_safe(function(x) formatC(x, format = "f", digits = 2))
f3 <- na_safe(function(x) formatC(x, format = "f", digits = 3))
f4 <- na_safe(function(x) formatC(x, format = "f", digits = 4))
fp <- na_safe(function(p) ifelse(p < 1e-4, formatC(p, format = "e", digits = 1),
                                 formatC(p, format = "f", digits = 4)))
fs <- na_safe(function(x, d = 3) formatC(signif(x, d), format = "g", digits = d))
fi <- na_safe(function(x) format(x, big.mark = ",", trim = TRUE))
pct <- na_safe(function(x, d = 1) {
  s <- formatC(100 * x, format = "f", digits = d)
  s <- sub("^-(0(\\.0+)?)$", "\\1", s)          # a value that rounds to zero is "0%", not "-0%"
  paste0(s, "%") })
