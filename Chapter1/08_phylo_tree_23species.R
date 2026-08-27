# ============================================================================
# Phylogenetic tree of 23 Borreliella species
# Description: Pruned from 287-genome core-gene tree to 23 representative
#              genomes (one per species). No legend. Minimal comments.
# Requirements: R >= 4.0, ggtree >= 3.0, tidyverse >= 1.3, ape, phangorn

# ============================================================================

library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)
library(dplyr)
library(readxl)

# --- Configuration: modify these paths to match your local environment ---
PATH_S1   <- "./data/Supplementary_Table_S1.xls"
PATH_TREE <- "./data/tree_287_species.treefile"
PATH_OUT  <- "./output/phylogenetic_tree_23species.pdf"

# --- Color palette ---
# Four pathogenic species follow the manuscript color scheme;
# remaining species are shown in light gray.
COLORS <- c(
  "Borreliella burgdorferi" = "#1F618D",
  "Borreliella garinii"     = "#D35400",
  "Borreliella afzelii"     = "#C0392B",
  "Borreliella bavariensis" = "#117A65",
  "Other species"           = "#BDC3C7"
)

# --- Core workflow ---
s1 <- read_excel(PATH_S1)
tree <- midpoint(read.tree(PATH_TREE))

species_counts <- s1 %>% count(species, name = "n_total")

best_ref <- s1 %>%
  group_by(species) %>%
  mutate(
    ref_score = (assembly_level == "Complete Genome") * 100 +
      (is_refseq == TRUE) * 50 +
      ifelse(is.na(QUAST_N50_kb), 0, QUAST_N50_kb) * 0.1 +
      ifelse(is.na(`BUSCO_Complete_%`), 0, `BUSCO_Complete_%`) * 0.5
  ) %>%
  slice_max(order_by = ref_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(species_counts, by = "species") %>%
  select(species, Genome_ID, n_total)

ref_tips <- paste0(best_ref$Genome_ID, "_genomic")
ref_tips <- ref_tips[ref_tips %in% tree$tip.label]
sub_tree <- keep.tip(tree, ref_tips)

tip_clean <- sub("_genomic$", "", sub_tree$tip.label)
match_idx <- match(tip_clean, best_ref$Genome_ID)
if(any(is.na(match_idx))) stop("Tip-label matching failed.")

sub_species <- best_ref$species[match_idx]
sub_n       <- best_ref$n_total[match_idx]

tip_anno <- data.frame(
  label   = sub_tree$tip.label,
  species = sub_species,
  n_total = sub_n,
  stringsAsFactors = FALSE
) %>%
  mutate(
    species_group = factor(
      case_when(
        species %in% c("Borreliella burgdorferi", "Borreliella garinii",
                       "Borreliella afzelii", "Borreliella bavariensis") ~ species,
        TRUE ~ "Other species"
      ),
      levels = c("Borreliella burgdorferi", "Borreliella garinii",
                 "Borreliella afzelii", "Borreliella bavariensis",
                 "Other species")
    ),
    label_text = case_when(
      species == "Borreliella burgdorferi" ~ paste0("B. burgdorferi (n = ", n_total, ")"),
      species == "Borreliella garinii"     ~ paste0("B. garinii (n = ", n_total, ")"),
      species == "Borreliella afzelii"     ~ paste0("B. afzelii (n = ", n_total, ")"),
      species == "Borreliella bavariensis" ~ paste0("B. bavariensis (n = ", n_total, ")"),
      TRUE ~ gsub("Borreliella ", "B. ", species)
    ),
    label_face = ifelse(species_group != "Other species", "bold", "plain")
  )

tree_data <- fortify(sub_tree)
bs_data <- tree_data %>%
  filter(!isTip, !is.na(label)) %>%
  mutate(bs = suppressWarnings(as.numeric(label))) %>%
  filter(!is.na(bs), bs >= 70)

p <- ggtree(sub_tree, layout = "rectangular", size = 0.5) %<+% tip_anno +
  aes(color = species_group) +
  scale_color_manual(values = COLORS, name = NULL) +
  geom_text(data = bs_data, aes(x = x, y = y, label = bs),
            hjust = 1.25, vjust = -0.5, size = 2.5, color = "#555555",
            inherit.aes = FALSE) +
  geom_tiplab(aes(label = label_text, fontface = label_face),
              offset = 0.003, size = 3.3, align = TRUE, linesize = 0.2) +
  geom_treescale(x = 0.005, y = 1, linesize = 0.5, fontsize = 3) +
  theme_tree2() +
  theme(legend.position = "none",
        plot.margin = margin(10, 120, 10, 10)) +
  xlim(0, max(tree_data$x) * 1.35)

ggsave(PATH_OUT, p, width = 10, height = 10, dpi = 300)