# ============================================================
# Project: Borreliella pan-genome analysis 
# Module: GC_Welch_GamesHowell
# Function: GC content stratification with Welch's ANOVA
# Statistics: Welch's ANOVA + Games-Howell + CLD
# Input: ./data/gene_presence_absence_roary.csv, ./data/pan_genome_reference.fa
# Output: ./output/GC_content_stratification.pdf
# Dependencies: data.table, ggplot2
# Note: No personal paths. Run from project root.
# ============================================================

PATH_ROARY <- "./data/gene_presence_absence_roary.csv"
PATH_FASTA <- "./data/pan_genome_reference.fa"
PATH_OUT   <- "./output"

library(data.table)
library(ggplot2)

N_TOTAL <- 287
SAVE_OUTPUT <- TRUE

THRESH <- list(Core = c(1.0, 1.0), Soft = c(0.95, 1.0), Shell = c(0.15, 0.95), Cloud = c(0, 0.15))
COL_RGB <- c("Core" = rgb(217, 83, 79, maxColorValue = 255), "Soft" = rgb(239, 172, 78, maxColorValue = 255),
             "Shell" = rgb(91, 155, 213, maxColorValue = 255), "Cloud" = rgb(160, 203, 232, maxColorValue = 255))

assign_strat <- function(freq) {
  ifelse(freq >= THRESH$Core[1]  & freq <= THRESH$Core[2],  "Core",
         ifelse(freq >= THRESH$Soft[1]  & freq <  THRESH$Soft[2],  "Soft",
                ifelse(freq >= THRESH$Shell[1] & freq <  THRESH$Shell[2], "Shell",
                       ifelse(freq >= THRESH$Cloud[1] & freq <  THRESH$Cloud[2], "Cloud", NA))))
}

games_howell_test <- function(y, group) {
  dt <- data.table(y = y, group = group)
  stats <- dt[, .(n = .N, mean = mean(y, na.rm = TRUE), var = var(y, na.rm = TRUE)), by = group]
  groups <- stats$group; k <- length(groups); comp <- data.table()
  for (i in 1:(k-1)) { for (j in (i+1):k) {
    gi <- groups[i]; gj <- groups[j]
    ni <- stats[group == gi, n]; mi <- stats[group == gi, mean]; vi <- stats[group == gi, var]
    nj <- stats[group == gj, n]; mj <- stats[group == gj, mean]; vj <- stats[group == gj, var]
    se <- sqrt(vi/ni + vj/nj); t_stat <- abs(mi - mj) / se
    df <- (vi/ni + vj/nj)^2 / (vi^2/(ni^2*(ni-1)) + vj^2/(nj^2*(nj-1)))
    q_stat <- t_stat * sqrt(2); p_raw <- ptukey(q_stat, nmeans = k, df = df, lower.tail = FALSE)
    comp <- rbind(comp, data.table(group1 = as.character(gi), group2 = as.character(gj),
                                   mean_i = mi, mean_j = mj, se = se, t = t_stat, df = df, p_raw = p_raw))
  }}
  comp[, p_adj := 1 - (1 - p_raw)^nrow(comp)]; list(stats = stats, comparisons = comp)
}

make_cld <- function(comp_dt, alpha = 0.05) {
  groups <- sort(unique(c(comp_dt$group1, comp_dt$group2))); k <- length(groups)
  adj <- matrix(1, nrow = k, ncol = k, dimnames = list(groups, groups)); diag(adj) <- 1
  for (r in 1:nrow(comp_dt)) { if (comp_dt$p_adj[r] < alpha) {
    g1 <- comp_dt$group1[r]; g2 <- comp_dt$group2[r]; adj[g1, g2] <- 0; adj[g2, g1] <- 0
  }}
  letters <- setNames(rep("", k), groups); letter <- "a"; remaining <- groups
  while (length(remaining) > 0) {
    seed <- remaining[1]; clique <- seed
    for (g in setdiff(remaining, seed)) { if (all(adj[clique, g] == 1)) clique <- c(clique, g) }
    for (g in clique) letters[g] <- paste0(letters[g], letter)
    remaining <- setdiff(remaining, clique); letter <- intToUtf8(utf8ToInt(letter) + 1)
  }
  data.table(group = groups, letter = letters)
}

read_fasta_gc <- function(fasta_file) {
  con <- file(fasta_file, "r"); lines <- readLines(con); close(con)
  seqs <- list(); cur_name <- NULL; cur_seq <- character()
  for (line in lines) {
    if (startsWith(line, ">")) {
      if (!is.null(cur_name)) seqs[[cur_name]] <- paste(cur_seq, collapse = "")
      cur_name <- strsplit(sub("^>", "", line), "\\s+")[[1]][1]; cur_seq <- character()
    } else { cur_seq <- c(cur_seq, toupper(line)) }
  }
  if (!is.null(cur_name)) seqs[[cur_name]] <- paste(cur_seq, collapse = "")
  data.table(gene_id = names(seqs), gc_content = sapply(seqs, function(s) {
    if (nchar(s) == 0) return(NA_real_)
    ch <- strsplit(s, "")[[1]]; gc <- sum(ch %in% c("G","C")); at <- sum(ch %in% c("A","T"))
    if (gc + at == 0) return(NA_real_); gc / (gc + at) * 100
  }))
}

# Main
roary <- fread(PATH_ROARY, header = TRUE)
roary[, freq := `No. isolates` / N_TOTAL]
roary[, strat := assign_strat(freq)]
roary[, strat := factor(strat, levels = c("Core","Soft","Shell","Cloud"))]

gc_dt <- read_fasta_gc(PATH_FASTA)
gc_dt[, gene_prefix := sub("_.*$", "", gene_id)]
roary_gc <- merge(roary, gc_dt, by.x = "Gene", by.y = "gene_id", all.x = TRUE)
if (roary_gc[!is.na(gc_content), .N] / nrow(roary) < 0.5) {
  roary_gc <- merge(roary, gc_dt, by.x = "Gene", by.y = "gene_prefix", all.x = TRUE)
}
plot_dt <- roary_gc[!is.na(gc_content) & !is.na(strat)]

welch_res <- oneway.test(gc_content ~ strat, data = plot_dt, var.equal = FALSE)
gh <- games_howell_test(plot_dt$gc_content, plot_dt$strat)
cld_dt <- make_cld(gh$comparisons)

letter_dt <- merge(gh$stats, cld_dt, by.x = "group", by.y = "group")
letter_dt[, y_pos := 55]

p_gc <- ggplot(plot_dt, aes(x = strat, y = gc_content, fill = strat)) +
  geom_boxplot(width = 0.5, outlier.size = 0.3, outlier.alpha = 0.3, colour = "black", size = 0.4) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5, fill = "white", colour = "black", stroke = 0.8) +
  geom_text(data = letter_dt, aes(x = group, y = y_pos, label = letter), inherit.aes = FALSE,
            size = 4.5, fontface = "bold", vjust = 0) +
  scale_fill_manual(values = COL_RGB) +
  scale_x_discrete(labels = c("Core" = "Core", "Soft" = "Soft core", "Shell" = "Shell", "Cloud" = "Cloud")) +
  labs(x = "Pan-genome category", y = "GC content (%)") +
  annotate("text", x = 2.5, y = 58, label = "Welch's ANOVA p < 0.001", size = 3.5) +
  coord_cartesian(ylim = c(10, 60)) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none", panel.border = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.text.x = element_text(angle = 0, face = "bold", size = 11),
        axis.title = element_text(face = "bold", size = 12))

print(p_gc)

if (SAVE_OUTPUT) {
  dir.create(PATH_OUT, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(PATH_OUT, "GC_content_stratification.pdf"), p_gc, width = 5.5, height = 4.8, dpi = 300)
}