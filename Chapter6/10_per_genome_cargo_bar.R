# ============================================================
# 10. Per-genome plasmid cargo decomposition for 135 strains (Fig 6e)
# Input : data/fig6a_function_composition_135.csv   (from script 02)
#             columns: genome, location, category, N, sp, total, prop
#             NB: total repeats across the 5 category rows of a
#                 genome-location; unique before use
#         data/strain_plasmid_status_master_corrected.csv (diagnostic 5 only;
#             skipped automatically when the columns are absent)
# Output: results/fig6e_integrated_cargo.pdf
# Note  : the per-genome expansion of the Fig 6a mean-composition pattern,
#         supporting the narrative that gar contraction mainly offloads
#         antigens and mobile elements, that the high bav burden is strongly
#         conserved, and that bb forms a continuous gradient
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

DIR_DATA <- "data"
DIR_RES  <- "results"
dir.create(DIR_RES, showWarnings = FALSE)

comp <- fread(file.path(DIR_DATA, "fig6a_function_composition_135.csv"))
comp[, sp := tolower(sp)]
stopifnot(uniqueN(comp$genome) == 135)

cat_levels <- c("Metabolic/Other", "Surface-associated", "Hypothetical",
                "MGE/Plasmid", "Antigen")
# Final palette (identical to Fig 6a)
cat_cols <- c("Metabolic/Other" = "#5DADE2", "Surface-associated" = "#1ABC9C",
              "Hypothetical" = "#E74C3C", "MGE/Plasmid" = "#F39C12",
              "Antigen" = "#7D3C98")
comp[, category := factor(category, levels = cat_levels)]
stopifnot(!any(is.na(comp$category)))

pl <- comp[location == "Plasmid"]

# ---------- sanity anchors (previously verified values) ----------
pl_tot <- unique(pl[, .(genome, sp, total)])
anchor_pg <- c(bb = 421, bav = 399, afz = 340, gar = 211)
pgm <- pl_tot[, .(mean_pg = mean(total)), by = sp]
print(pgm)
for (s in names(anchor_pg)) stopifnot(abs(pgm[sp == s, mean_pg] - anchor_pg[[s]]) < 1)

pool_ag <- pl[category == "Antigen", sum(N)] / pl_tot[, sum(total)]
stopifnot(abs(pool_ag - 0.0643) < 0.002)
gar_ag <- pl[sp == "gar" & category == "Antigen", mean(prop)]
stopifnot(abs(gar_ag - 0.0383) < 0.005)
cat("sanity: plasmid gene counts / pooled antigen 6.43% / gar antigen 3.83% anchors all passed\n")

# ---------- ordering: species bb -> gar -> afz -> bav, descending plasmid
#   gene count within species ----------
sp_levels <- c("bb", "gar", "afz", "bav")
pl_tot[, sp := factor(sp, levels = sp_levels)]
ord <- pl_tot[order(sp, -total), genome]
pl[, genome := factor(genome, levels = ord)]
pl[, spf := factor(sp, levels = sp_levels,
                   labels = c("B. burgdorferi", "B. garinii",
                              "B. afzelii†", "B. bavariensis"))]

# ---------- plot ----------
p <- ggplot(pl, aes(x = genome, y = N, fill = category)) +
  geom_col(width = 0.9, color = NA) +
  facet_wrap(~spf, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = cat_cols, name = "Function",
                    breaks = c("Antigen", "MGE/Plasmid", "Hypothetical",
                               "Surface-associated", "Metabolic/Other")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Plasmid-borne genes (N)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "gray95", color = "black"),
        strip.text  = element_text(face = "italic", size = 11),
        panel.grid  = element_blank(),
        panel.spacing = unit(0.6, "lines"),
        legend.position = "bottom")
print(p)
ggsave(file.path(DIR_RES, "fig6e_integrated_cargo.pdf"),
       p, width = 13.5, height = 4.8)

# ---------- narrative diagnostics + final value locks (main-text run values) ----------
cat("\n=== diagnostic 1: plasmid burden distribution ===\n")
burden <- pl_tot[, .(n      = .N,
                     median = as.numeric(median(total)),
                     q1     = as.numeric(quantile(total, 0.25)),
                     q3     = as.numeric(quantile(total, 0.75)),
                     min    = as.numeric(min(total)),
                     max    = as.numeric(max(total))), by = sp]
print(burden)
# Main-text locks: bb median 418 (26-696); gar median 193 (108-320);
#   bav median 423, 21 of 22 strains concentrated in 361-446, outlier-poor strain 102
stopifnot(burden[sp == "bb", median] == 418,
          burden[sp == "bb", min] == 26, burden[sp == "bb", max] == 696,
          abs(burden[sp == "gar", median] - 193) <= 0.5,   # n=16 even, median may be x.5
          burden[sp == "gar", min] == 108, burden[sp == "gar", max] == 320,
          abs(burden[sp == "bav", median] - 423) <= 0.5,   # n=22 even, median may be x.5
          burden[sp == "bav", min] == 102,
          pl_tot[sp == "bav", sum(total >= 361 & total <= 446)] == 21)

cat("\n=== diagnostic 2: compartment means (gar/bb ratios: antigen 30%, MGE 49%, surface/hypothetical 66%) ===\n")
catN <- pl[, .(meanN = mean(N), sdN = sd(N)), by = .(sp, category)]
print(dcast(catN, category ~ sp, value.var = "meanN"))
mg <- setNames(catN[sp == "gar", meanN], catN[sp == "gar", category])
mb <- setNames(catN[sp == "bb",  meanN], catN[sp == "bb",  category])
stopifnot(abs(mg[["Antigen"]] - 8.3) < 0.1,
          abs(mg[["Antigen"]]          / mb[["Antigen"]]          - 0.30) < 0.015,
          abs(mg[["MGE/Plasmid"]]      / mb[["MGE/Plasmid"]]      - 0.49) < 0.015,
          abs(mg[["Surface-associated"]]/ mb[["Surface-associated"]] - 0.66) < 0.015,
          abs(mg[["Hypothetical"]]     / mb[["Hypothetical"]]     - 0.66) < 0.015)

cat("\n=== diagnostic 3: bavariensis plasmid gene counts, descending ===\n")
print(pl_tot[sp == "bav", total][order(-pl_tot[sp == "bav", total])])

cat("\n=== diagnostic 4: bb burden, top/bottom 10 strains ===\n")
bb_sorted <- pl_tot[sp == "bb", total][order(-pl_tot[sp == "bb", total])]
print(list(head = head(bb_sorted, 10), tail = tail(bb_sorted, 10)))

cat("\n=== diagnostic 5: identity check of plasmid-poor strains (<150) ===\n")
master <- fread(file.path(DIR_DATA, "strain_plasmid_status_master_corrected.csv"))
poor <- pl_tot[total < 150]
if (all(c("n_contigs", "asm_level") %in% names(master))) {
  print(merge(poor, master[, .(Genome_ID, n_contigs, asm_level)],
              by.x = "genome", by.y = "Genome_ID", all.x = TRUE))
} else {
  cat("master lacks n_contigs/asm_level columns; identity check skipped (can be filled from Table S1)\n")
  print(poor)
}
cat("batch genomes (GCF_02657/02664 bare-chromosome batch) present:",
    any(grepl("^GCF_02657|^GCF_02664", poor$genome)), "(should be FALSE)\n")
cat("=== 10 done ===\n")
