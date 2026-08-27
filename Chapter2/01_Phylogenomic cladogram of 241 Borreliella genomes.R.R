# Upload version: circular cladogram with dual outer rings
# Inner ring: Hypothetical protein proportion (purple gradient)
# Outer ring: Genome size (blue gradient)
# Requirements: ape, phangorn, ggtree, ggplot2, dplyr, readxl, ggnewscale

library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggnewscale)

# -------------------- INPUT PATHS (edit before running) --------------------
tree_path <- "PATH/TO/tree_241_four_species.treefile"
s1_path   <- "PATH/TO/Supplementary_Table_S1.xls"
gff_dir   <- "PATH/TO/gff3_four_species"

# -------------------- species order & colors --------------------
sp_colors <- c(
  "Borreliella burgdorferi" = "#1F618D",
  "Borreliella garinii"     = "#E74C3C",
  "Borreliella afzelii"     = "#1ABC9C",
  "Borreliella bavariensis" = "#5DADE2"
)
sp_levels <- names(sp_colors)

# -------------------- 1. Load metadata --------------------
s1 <- read_excel(s1_path, sheet = 1)
s1 <- as.data.frame(s1, stringsAsFactors = FALSE)
s1_four <- s1 %>% 
  filter(species %in% sp_levels) %>%
  mutate(Genome_ID = sub("_genomic$", "", Genome_ID))

# -------------------- 2. Tree: midpoint root + cladogram --------------------
tree_raw <- read.tree(tree_path)
tree_raw <- midpoint(tree_raw)
tree <- tree_raw
tree$edge.length <- rep(1, nrow(tree$edge))

# Match species annotation
tip_clean <- sub("_genomic$", "", tree$tip.label)
match_idx <- match(tip_clean, s1_four$Genome_ID)
if(any(is.na(match_idx))) {
  tip_pfx <- sub("^(GC[FA]_\\d+\\.\\d+).*", "\\1", tip_clean)
  ref_pfx <- sub("^(GC[FA]_\\d+\\.\\d+).*", "\\1", s1_four$Genome_ID)
  m2 <- match(tip_pfx, ref_pfx)
  match_idx[is.na(match_idx)] <- m2[is.na(match_idx)]
}
if(any(is.na(match_idx))) stop("Unmatched tips found")

tip_anno <- data.frame(
  label   = tree$tip.label,
  species = factor(s1_four$species[match_idx], levels = sp_levels),
  stringsAsFactors = FALSE
)

# Genome stats
genome_stats <- s1_four %>%
  select(Genome_ID, Genome_Size = QUAST_Length_Mb)

# -------------------- 3. Extract Hypothetical protein % --------------------
extract_hypo <- function(gid) {
  gff_path <- file.path(gff_dir, paste0(gid, "_genomic.gff3"))
  if (!file.exists(gff_path)) return(NA)
  lines <- readLines(gff_path)
  feat_lines <- lines[!grepl("^#", lines)]
  if (length(feat_lines) == 0) return(NA)
  types <- sapply(strsplit(feat_lines, "\t"), function(x) if(length(x) >= 3) x[3] else NA)
  attrs <- sapply(strsplit(feat_lines, "\t"), function(x) if(length(x) >= 9) x[9] else "")
  cds_n <- sum(types == "CDS", na.rm = TRUE)
  hypo_n <- sum(grepl("product=hypothetical protein", attrs, ignore.case = TRUE))
  if (cds_n == 0) return(NA)
  hypo_n / cds_n * 100
}

hypo_vec <- sapply(s1_four$Genome_ID, extract_hypo)
names(hypo_vec) <- s1_four$Genome_ID

# -------------------- 4. Bootstrap labels --------------------
tree_data <- fortify(tree)
node_burg <- MRCA(tree, tree$tip.label[tip_anno$species == "Borreliella burgdorferi"])
node_gar  <- MRCA(tree, tree$tip.label[tip_anno$species == "Borreliella garinii"])
node_afz  <- MRCA(tree, tree$tip.label[tip_anno$species == "Borreliella afzelii"])
node_bav  <- MRCA(tree, tree$tip.label[tip_anno$species == "Borreliella bavariensis"])

parent_nodes <- sapply(c(node_burg, node_gar, node_afz, node_bav), function(n) {
  tree_data$parent[tree_data$node == n]
})
all_show_nodes <- unique(c(node_burg, node_gar, node_afz, node_bav, parent_nodes))
all_show_nodes <- all_show_nodes[all_show_nodes > length(tree$tip.label)]

bs_data_raw <- fortify(tree_raw)
bs_data <- bs_data_raw %>%
  filter(!isTip, !is.na(label), node %in% all_show_nodes) %>%
  mutate(bs = suppressWarnings(as.numeric(label))) %>%
  filter(!is.na(bs), bs >= 70)

# -------------------- 5. Draw cladogram --------------------
p <- ggtree(tree, layout = "circular", size = 0.25) %<+% tip_anno +
  aes(color = species) +
  scale_color_manual(values = sp_colors, breaks = sp_levels, name = "Species") +
  geom_tippoint(aes(color = species), size = 0.9, alpha = 0.9) +
  geom_text2(data = bs_data, aes(subset = node %in% all_show_nodes, label = bs),
             size = 2.0, color = "#333333", vjust = -1.0) +
  theme(
    legend.position = "none", plot.margin = margin(20, 20, 20, 20),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(), axis.line = element_blank(),
    axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()
  ) +
  xlim(0, max(tree_data$x) * 1.28)

# -------------------- 6. Outer ring data --------------------
tip_df <- p$data[p$data$isTip, c("x", "y", "label")]
tip_df$genome_size <- genome_stats$Genome_Size[match(sub("_genomic$", "", tip_df$label), genome_stats$Genome_ID)]
tip_df$hypo_ratio  <- hypo_vec[sub("_genomic$", "", tip_df$label)]
tip_df$species     <- tip_anno$species[match(tip_df$label, tip_anno$label)]
tip_df <- tip_df[!is.na(tip_df$genome_size) & !is.na(tip_df$hypo_ratio) & !is.na(tip_df$species), ]
tip_df <- tip_df[order(tip_df$y), ]

n <- nrow(tip_df)
y_sorted <- tip_df$y
if (n > 1) {
  gaps <- diff(y_sorted)
  tip_df$ymin <- y_sorted - c(gaps[1], gaps) / 2
  tip_df$ymax <- y_sorted + c(gaps, gaps[length(gaps)]) / 2
} else {
  tip_df$ymin <- y_sorted - 0.5
  tip_df$ymax <- y_sorted + 0.5
}

r_max <- max(tip_df$x)

# Inner ring: Hypothetical protein % (fixed width)
offset1 <- r_max * 0.03
width1  <- r_max * 0.08
tip_df$xmin1 <- r_max + offset1
tip_df$xmax1 <- r_max + offset1 + width1

# Outer ring: Genome size (fixed width)
gap <- r_max * 0.012
width2 <- r_max * 0.08
tip_df$xmin2 <- tip_df$xmax1 + gap
tip_df$xmax2 <- tip_df$xmin2 + width2

# -------------------- 7. Draw dual rings --------------------
p2 <- p + 
  # Inner ring: Hypothetical protein % (purple)
  geom_rect(data = tip_df,
            aes(xmin = xmin1, xmax = xmax1, ymin = ymin, ymax = ymax, fill = hypo_ratio),
            inherit.aes = FALSE, color = "white", size = 0.15) +
  scale_fill_gradient(
    low = "#F3E5F5",
    high = "#7B1FA2",
    name = "Hypothetical\nprotein (%)"
  ) +
  # Outer ring: Genome size (blue)
  new_scale_fill() +
  geom_rect(data = tip_df,
            aes(xmin = xmin2, xmax = xmax2, ymin = ymin, ymax = ymax, fill = genome_size),
            inherit.aes = FALSE, color = "white", size = 0.15) +
  scale_fill_gradient(
    low = "#E3F2FD",
    high = "#1565C0",
    name = "Genome size\n(Mb)"
  ) +
  # Compact legend on right
  theme(
    legend.position = "right",
    legend.justification = "top",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.margin = margin(2, 2, 2, 5),
    legend.spacing.y = unit(2, "mm"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7, face = "italic"),
    legend.key.size = unit(3, "mm"),
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.background = element_rect(fill = "white", colour = "grey80", size = 0.2),
    plot.margin = margin(15, 50, 15, 15)
  ) +
  guides(
    color = guide_legend(title = "Species", order = 1, override.aes = list(size = 2))
  )

print(p2)
cat("\nPreview done. Inner ring = purple (Hypothetical %), Outer ring = blue (Genome size).\n")