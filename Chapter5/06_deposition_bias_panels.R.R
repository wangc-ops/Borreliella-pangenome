########################################################################
## 06_deposition_bias_panels.R
## Plasmid-sequence deposition bias in the full genome set (two panels).
## Depends on: 01_build_pav_matrices.R
##   (results/antigen_pav_matrix_all.csv, results/antigen_pav_matrix_carriers.csv)
## Panel A: deposited vs not-deposited genome counts per species
## Panel B: family detection rates in the full set vs the carrier subset
## Sanity check printed: the full-set ospC detection rate must equal the
## carrier fraction in every species.
## Output: results/deposition_bias_panels.pdf
########################################################################
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

dir.create("results", showWarnings = FALSE)

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

w_all      <- fread("results/antigen_pav_matrix_all.csv")
w_carriers <- fread("results/antigen_pav_matrix_carriers.csv")
for (dt in list(w_all, w_carriers)) {
  spcol <- intersect(names(dt), c("sp", "species"))[1]
  setnames(dt, spcol, "sp")
  dt[, sp := to_sp_code(sp)]
  dt[, sp := factor(sp, levels = c("bb", "gar", "afz", "bav"))]
}
fams <- c("OspA","OspB","OspC","OspD","DbpA","DbpB","Erp","Mlp",
          "CspA_Z","BBK32","VlsE","BmpA_B","P66")

## Sanity check: carrier fraction vs full-set ospC detection rate
chk <- w_all[, .(n = .N, carriers = sum(carrier),
                 carrier_frac  = round(mean(carrier), 3),
                 ospC_rate_all = round(mean(OspC), 3)), by = sp]
cat("== Check: ospC_rate_all must equal carrier_frac ==\n"); print(chk)

pal    <- c(bb = "#1F618D", gar = "#E74C3C", afz = "#1ABC9C", bav = "#5DADE2")
sp_lab <- c(bb = "B. burgdorferi", gar = "B. garinii",
            afz = "B. afzelii", bav = "B. bavariensis")

## Panel A
stat <- w_all[, .(n = .N), by = .(sp, carrier)]
stat[, Status := factor(ifelse(carrier == 1, "Deposited", "Not deposited"),
                        levels = c("Deposited", "Not deposited"))]
pA <- ggplot(stat, aes(sp, n, fill = Status)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c("Deposited" = "#7FB3D5", "Not deposited" = "#E59866")) +
  scale_x_discrete(labels = sp_lab) +
  labs(x = NULL, y = "Genomes") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(face = "italic"),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank())

## Panel B
rr <- rbindlist(lapply(fams, function(f)
  rbind(w_all[,      .(family = f, set = "all",      rate = 100 * mean(get(f))), by = sp],
        w_carriers[, .(family = f, set = "carriers", rate = 100 * mean(get(f))), by = sp])))
rrw <- dcast(rr, family + sp ~ set, value.var = "rate")
setnames(rrw, c("all", "carriers"), c("rate_all", "rate_carriers"))
rrw[, loc := ifelse(family %in% c("BmpA_B", "P66"), "Chromosome", "Plasmid")]
pB <- ggplot(rrw, aes(rate_all, rate_carriers, color = sp, shape = loc)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = pal, labels = sp_lab, name = NULL) +
  scale_shape_manual(values = c(Chromosome = 17, Plasmid = 16), name = NULL) +
  labs(x = "Detection rate in the full genome set (%)",
       y = "Detection rate in plasmid-sequence carriers (%)") +
  theme_classic(base_size = 10) +
  theme(legend.text = element_text(face = "italic"),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank())

print(pA | pB)
ggsave("results/deposition_bias_panels.pdf", pA | pB,
       width = 14, height = 5, device = cairo_pdf)
cat("Saved: results/deposition_bias_panels.pdf\n")
