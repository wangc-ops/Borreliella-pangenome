########################################################################
## 03_copynumber_erp_mlp_stats.R
## erp/mlp copy numbers in plasmid-sequence carriers: two boxplots plus
## a single multi-sheet workbook of statistics.
## Statistics: Welch ANOVA + Games-Howell on bb/gar/bav; CLD letters
##             remapped so that "a" labels the highest mean.
## afz (n = 6) is display-only: excluded from tests, no letters, axis
## label carries a dagger.
## Requires: rstatix, multcompView, writexl
## Outputs: results/erp_mlp_copynumber_stats.xlsx (4 sheets)
##          results/boxplot_erp_copies.pdf, results/boxplot_mlp_copies.pdf
########################################################################
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
  library(rstatix); library(multcompView); library(writexl)
})

dir.create("results", showWarnings = FALSE)

copy_path   <- "data/antigen_copy_numbers.csv"
master_path <- "data/strain_plasmid_status_master_corrected.csv"

to_pres <- function(x) {
  if (is.logical(x))   return(x)
  if (is.character(x)) return(tolower(trimws(x)) %in% c("1", "true"))
  as.numeric(x) > 0
}

to_sp_code <- function(x) {
  x <- trimws(tolower(as.character(x)))
  if (all(x %in% c("bb", "gar", "afz", "bav"))) return(x)
  out <- rep(NA_character_, length(x))
  out[grepl("burgdorferi", x)] <- "bb"
  out[grepl("garinii", x)]     <- "gar"
  out[grepl("afzelii", x)]     <- "afz"
  out[grepl("bavariensis", x)] <- "bav"
  if (any(is.na(out))) stop("Unmappable species name(s): ",
                            paste(unique(x[is.na(out)]), collapse = ", "))
  out
}

## Remap multcompLetters so that "a" marks the group with the highest mean
remap_by_mean <- function(cld_named, mean_named) {
  ord <- names(sort(-mean_named[names(cld_named)]))
  cur <- cld_named
  for (i in seq_along(ord)) {
    tgt <- letters[i]; g <- ord[i]
    gl  <- strsplit(cur[g], "", fixed = TRUE)[[1]]
    if (tgt %in% gl) next
    from <- setdiff(gl, letters[seq_len(i - 1)])[1]
    cur <- vapply(cur, function(s) {
      ch <- strsplit(s, "", fixed = TRUE)[[1]]
      ch[ch == from] <- "\x01"; ch[ch == tgt] <- from; ch[ch == "\x01"] <- tgt
      paste0(ch, collapse = "")
    }, character(1))
  }
  vapply(strsplit(cur, "", fixed = TRUE),
         function(x) paste0(sort(x), collapse = ""), character(1))
}

cn <- fread(copy_path)
cn[, genome := sub("_genomic$", "", Genome_ID)]
mast <- fread(master_path)
mast[, genome  := sub("_genomic$", "", Genome_ID)]
mast[, sp      := to_sp_code(species)]
mast[, carrier := to_pres(has_plasmid_corrected)]

d <- merge(cn[, .(genome, erp = as.numeric(erp_copies), mlp = as.numeric(mlp_copies))],
           mast[, .(genome, sp, carrier)], by = "genome")
d <- d[carrier == TRUE]
stopifnot(nrow(d) == 135)

long <- rbind(d[, .(genome, sp, trait = "Erp", value = erp)],
              d[, .(genome, sp, trait = "Mlp", value = mlp)])
long[, sp := factor(sp, levels = c("bb", "gar", "afz", "bav"))]

## ---- Descriptive statistics (afz row is display-only, daggered) ----
desc <- long[, .(n       = .N,
                 mean    = round(mean(value), 2),
                 sd      = round(sd(value), 2),
                 median  = as.numeric(median(value)),
                 Q1      = as.numeric(quantile(value, 0.25)),
                 Q3      = as.numeric(quantile(value, 0.75))),
             by = .(trait, sp)]
cat("== Descriptive statistics (afz display-only) ==\n"); print(desc)

## ---- Welch ANOVA + Games-Howell + CLD (bb/gar/bav only) ----
welch_rows <- list(); gh_rows <- list(); cld_rows <- list(); letters_list <- list()
for (tr in c("Erp", "Mlp")) {
  dd <- long[trait == tr & sp %in% c("bb", "gar", "bav")]
  dd[, sp := droplevels(sp)]
  w  <- welch_anova_test(dd, value ~ sp)
  gh <- games_howell_test(dd, value ~ sp)
  cat("\n== ", tr, " Welch ANOVA ==\n", sep = ""); print(as.data.frame(w))
  cat("== ", tr, " Games-Howell (direction judged by group means) ==\n", sep = "")
  print(as.data.frame(gh)[, c("group1", "group2", "estimate", "p.adj")])
  pvec <- gh$p.adj; names(pvec) <- paste(gh$group1, gh$group2, sep = "-")
  cld   <- multcompLetters(pvec)$Letters
  means <- dd[, mean(value), by = sp][, setNames(V1, sp)]
  cld2  <- remap_by_mean(cld, means)
  cat(tr, " CLD before remap:", paste(names(cld), cld, sep = "=", collapse = " "),
      " -> after:", paste(names(cld2), cld2, sep = "=", collapse = " "), "\n")
  welch_rows[[tr]] <- data.frame(trait = tr, as.data.frame(w))
  gh_rows[[tr]]    <- data.frame(trait = tr, as.data.frame(gh))
  cld_rows[[tr]]   <- data.frame(trait = tr, sp = names(cld2), cld = unname(cld2))
  letters_list[[tr]] <- data.frame(sp = names(cld2), cld = unname(cld2),
                                   stringsAsFactors = FALSE)
}

## ---- Single workbook with four sheets ----
write_xlsx(
  list(descriptive  = as.data.frame(desc),
       welch        = as.data.frame(rbindlist(welch_rows)),
       games_howell = as.data.frame(rbindlist(gh_rows)),
       cld          = as.data.frame(rbindlist(cld_rows))),
  "results/erp_mlp_copynumber_stats.xlsx")
cat("Saved: results/erp_mlp_copynumber_stats.xlsx (sheets: descriptive / welch / games_howell / cld)\n")

## ---- Plots (two separate figures) ----
pal    <- c(bb = "#1F618D", gar = "#E74C3C", afz = "#1ABC9C", bav = "#5DADE2")
sp_lab <- c(bb = "B. burgdorferi", gar = "B. garinii",
            afz = "B. afzelii †", bav = "B. bavariensis")

make_panel <- function(tr, ylab) {
  dd <- long[trait == tr]
  L  <- letters_list[[tr]]
  L$sp <- factor(L$sp, levels = levels(dd$sp))
  off <- 0.06 * max(dd$value)
  L$y  <- vapply(as.character(L$sp),
                 function(s) max(dd$value[dd$sp == s]) + off, numeric(1))
  ggplot(dd, aes(sp, value)) +
    geom_boxplot(aes(fill = sp), alpha = 0.85, outlier.shape = NA,
                 fatten = 1.5, linewidth = 0.5, color = "black") +
    geom_jitter(aes(fill = sp), shape = 21, color = "black",
                width = 0.12, size = 1.8, alpha = 0.85, stroke = 0.25) +
    geom_text(data = L, aes(sp, y, label = cld), size = 4.5, fontface = "bold") +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_discrete(labels = sp_lab) +
    labs(x = "Species", y = ylab) +
    theme_classic(base_size = 10) +
    theme(axis.text.x  = element_text(face = "bold.italic", angle = 45,
                                      hjust = 1, size = 9.5, color = "black"),
          axis.text.y  = element_text(color = "black"),
          axis.title   = element_text(face = "bold", size = 12),
          axis.line    = element_line(color = "black", linewidth = 0.6),
          axis.ticks   = element_line(color = "black", linewidth = 0.5),
          panel.grid   = element_blank())
}

p_erp <- make_panel("Erp", "Erp copies per genome")
p_mlp <- make_panel("Mlp", "Mlp copies per genome")
print(p_erp)
print(p_mlp)
ggsave("results/boxplot_erp_copies.pdf", p_erp,
       width = 4.5, height = 4.5, device = cairo_pdf)
ggsave("results/boxplot_mlp_copies.pdf", p_mlp,
       width = 4.5, height = 4.5, device = cairo_pdf)
cat("Saved: results/boxplot_erp_copies.pdf, results/boxplot_mlp_copies.pdf\n")
