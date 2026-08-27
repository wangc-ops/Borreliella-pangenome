########################################################################
## 05_erp_cp32_correlation.R
## erp copy number vs cp32-like plasmid count (carrier subset).
## Spearman correlation + bootstrap 95% CI (2000 replicates, set.seed(42));
## tested in bb/gar/bav only. The afz panel is display-only (no test
## annotation) and its strip label carries a dagger.
## Outputs: results/erp_cp32_spearman.csv, results/scatter_erp_vs_cp32.pdf
########################################################################
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

dir.create("results", showWarnings = FALSE)

copy_path   <- "data/antigen_copy_numbers.csv"
master_path <- "data/strain_plasmid_status_master_corrected.csv"
cp_path     <- "data/cp32_like_burden.csv"

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

cn <- fread(copy_path);  cn[, genome := sub("_genomic$", "", Genome_ID)]
mast <- fread(master_path)
mast[, genome  := sub("_genomic$", "", Genome_ID)]
mast[, sp      := to_sp_code(species)]
mast[, carrier := to_pres(has_plasmid_corrected)]
cp <- fread(cp_path)
stopifnot(all(c("genome_id", "cp32_count") %in% names(cp)))
cp[, genome := sub("_genomic$", "", genome_id)]

dd <- merge(cn[, .(genome, erp = as.numeric(erp_copies))],
            mast[, .(genome, sp, carrier)], by = "genome")
dd <- merge(dd[carrier == TRUE], cp[, .(genome, cp32 = as.numeric(cp32_count))],
            by = "genome")
stopifnot(nrow(dd) == 135)
dd[, sp := factor(sp, levels = c("bb", "gar", "afz", "bav"))]

set.seed(42)
res <- rbindlist(lapply(c("bb", "gar", "bav"), function(s) {
  x <- dd$erp[dd$sp == s]; y <- dd$cp32[dd$sp == s]
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  boot <- replicate(2000, {
    i <- sample.int(length(x), replace = TRUE)
    suppressWarnings(cor(x[i], y[i], method = "spearman"))
  })
  ci <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)
  data.table(sp = s, n = length(x), rho = round(unname(ct$estimate), 3),
             CI_lo = round(ci[1], 3), CI_hi = round(ci[2], 3), P = ct$p.value)
}))
cat("== erp x cp32 Spearman (carrier subset, bb/gar/bav) ==\n"); print(res)
fwrite(res, "results/erp_cp32_spearman.csv")

p_lab <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
res[, lab := sprintf("ρ = %.3f\n%s", rho, vapply(P, p_lab, character(1)))]
res[, sp := factor(sp, levels = levels(dd$sp))]

pal    <- c(bb = "#1F618D", gar = "#E74C3C", afz = "#1ABC9C", bav = "#5DADE2")
sp_lab <- c(bb = "B. burgdorferi", gar = "B. garinii",
            afz = "B. afzelii †", bav = "B. bavariensis")

p <- ggplot(dd, aes(cp32, erp, color = sp)) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
  geom_text(data = res, aes(label = lab), x = -Inf, y = Inf,
            hjust = -0.1, vjust = 1.3, size = 3, inherit.aes = FALSE) +
  facet_wrap(~ sp, nrow = 1, labeller = as_labeller(sp_lab)) +
  scale_color_manual(values = pal, guide = "none") +
  labs(x = "cp32-like plasmid count", y = "erp copies per genome") +
  theme_classic(base_size = 10) +
  theme(strip.text = element_text(face = "italic"),
        axis.line  = element_line(color = "black"),
        panel.grid = element_blank())
print(p)
ggsave("results/scatter_erp_vs_cp32.pdf", p,
       width = 10, height = 3.2, device = cairo_pdf)
cat("Saved: results/scatter_erp_vs_cp32.pdf, results/erp_cp32_spearman.csv\n")
