# ============================================================
# 04. Plasmid cluster category distribution by species
#     (Core / Shared / Species-specific composition of 483 clusters)
# Input : data/Supplementary_Table_S1.xls
#         data/strain_plasmid_status_master_corrected.csv
#         data/plasmid_clusters_cluster_nochrom.tsv (483 plasmid clusters, chromosomes removed)
# Output: results/cluster_category_stack.pdf
#         results/plasmid_cluster_distribution.csv   (Table S11 source)
# ============================================================
pkgs <- c("tidyverse", "readxl")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(tidyverse); library(readxl)

dir.create("results", showWarnings = FALSE)

sp_levels <- c("Borreliella burgdorferi", "Borreliella garinii",
               "Borreliella afzelii", "Borreliella bavariensis")
cat_colors <- c("Core (4 species)"     = "#D5DBDB",
                "Shared (2-3 species)" = "#F1948A",
                "Species-specific"     = "#85C1E9")

# half-up rounding to 1 decimal (R's round() is banker's rounding)
round1 <- function(x) floor(x * 10 + 0.5) / 10

master <- read_csv("data/strain_plasmid_status_master_corrected.csv", show_col_types = FALSE) %>%
  mutate(Genome_ID = sub("_genomic$", "", Genome_ID))
plasmid_ids <- master %>% filter(has_plasmid_corrected == TRUE) %>% pull(Genome_ID)

s1 <- read_excel("data/Supplementary_Table_S1.xls") %>%
  mutate(Genome_ID = sub("_genomic$", "", Genome_ID),
         species = case_when(
           grepl("burgdorferi", species, ignore.case = TRUE) ~ "Borreliella burgdorferi",
           grepl("garinii",     species, ignore.case = TRUE) ~ "Borreliella garinii",
           grepl("afzelii",     species, ignore.case = TRUE) ~ "Borreliella afzelii",
           grepl("bavariensis", species, ignore.case = TRUE) ~ "Borreliella bavariensis",
           TRUE ~ NA_character_)) %>%
  filter(!is.na(species)) %>% select(Genome_ID, species)

clusters <- read_tsv("data/plasmid_clusters_cluster_nochrom.tsv",
                     col_names = c("rep_id", "member_id"), show_col_types = FALSE) %>%
  mutate(Genome_ID = sub("_genomic$", "", str_extract(member_id, "^[^|]+")))

# representative-sequence size parsed from the "length=" tag in rep_id
rep_sizes <- clusters %>%
  distinct(rep_id) %>%
  mutate(size_kb = round1(as.numeric(str_extract(rep_id, "(?<=length=)[0-9]+")) / 1000))

pav_sp <- clusters %>%
  distinct(Genome_ID, rep_id) %>%
  mutate(present = 1) %>%
  pivot_wider(id_cols = Genome_ID, names_from = rep_id,
              values_from = present, values_fill = 0) %>%
  inner_join(s1, by = "Genome_ID") %>%
  filter(Genome_ID %in% plasmid_ids)
stopifnot(nrow(pav_sp) == 135, ncol(pav_sp) - 2 == 483)

# --- cluster categories by number of species detected ---
cluster_meta <- pav_sp %>%
  pivot_longer(-c(Genome_ID, species), names_to = "rep_id", values_to = "present") %>%
  group_by(rep_id, species) %>%
  summarise(has_it = any(present == 1), .groups = "drop") %>%
  group_by(rep_id) %>%
  summarise(n_species = sum(has_it), .groups = "drop") %>%
  mutate(category = case_when(
    n_species == 4 ~ "Core (4 species)",
    n_species == 1 ~ "Species-specific",
    TRUE ~ "Shared (2-3 species)"),
    category = factor(category, levels = c("Core (4 species)",
                                           "Shared (2-3 species)",
                                           "Species-specific")))
cat("Cluster category totals:\n"); print(cluster_meta %>% count(category))

# --- per-cluster distribution table (Table S11) ---
# rates: detection rate among plasmid carriers of each species
# (denominators: Bbu 91 / Bga 16 / Baf 6 / Bba 22), half-up rounded to 1 decimal
cluster_distribution <- pav_sp %>%
  pivot_longer(-c(Genome_ID, species), names_to = "rep_id", values_to = "present") %>%
  group_by(rep_id, species) %>%
  summarise(rate = round1(mean(present) * 100), .groups = "drop") %>%
  pivot_wider(names_from = species, values_from = rate, values_fill = 0) %>%
  left_join(cluster_meta, by = "rep_id") %>%
  left_join(rep_sizes, by = "rep_id") %>%
  mutate(
    species_present = paste(
      ifelse(`Borreliella burgdorferi` > 0, "Bbu", ""),
      ifelse(`Borreliella garinii` > 0, "Bga", ""),
      ifelse(`Borreliella afzelii` > 0, "Baf", ""),
      ifelse(`Borreliella bavariensis` > 0, "Bba", ""), sep = "-"),
    species_present = gsub("-+", "-", species_present),
    species_present = gsub("^-|-$", "", species_present)) %>%
  select(rep_id, size_kb,
         `Borreliella burgdorferi`, `Borreliella garinii`,
         `Borreliella afzelii`, `Borreliella bavariensis`,
         n_species, category, species_present) %>%
  arrange(desc(n_species), desc(size_kb))

# anchors: 483 clusters, none shared by all 4 species, exactly 5 shared by 2-3
stopifnot(nrow(cluster_distribution) == 483,
          sum(cluster_distribution$n_species == 4) == 0,
          sum(cluster_distribution$n_species >= 2) == 5)

write_csv(cluster_distribution, "results/plasmid_cluster_distribution.csv")

cat("Shared cluster detail:\n")
cluster_distribution %>% filter(category == "Shared (2-3 species)") %>%
  select(rep_id, size_kb, n_species, species_present,
         `Borreliella burgdorferi`, `Borreliella garinii`,
         `Borreliella afzelii`, `Borreliella bavariensis`) %>%
  print(width = Inf)

# --- per-genome category means (for stacked bar) ---
cat_summary <- pav_sp %>%
  pivot_longer(-c(Genome_ID, species), names_to = "rep_id", values_to = "present") %>%
  left_join(cluster_meta, by = "rep_id") %>%
  group_by(Genome_ID, species, category) %>%
  summarise(n = sum(present), .groups = "drop") %>%
  mutate(species = factor(species, levels = sp_levels)) %>%
  group_by(species, category) %>%
  summarise(mean_n = mean(n), se = sd(n) / sqrt(n()), .groups = "drop") %>%
  filter(mean_n > 0)
print(cat_summary)

fig <- ggplot(cat_summary, aes(x = species, y = mean_n, fill = category)) +
  geom_bar(stat = "identity", position = "stack",
           color = "black", linewidth = 0.4, width = 0.65) +
  scale_fill_manual(values = cat_colors, drop = FALSE) +
  scale_x_discrete(labels = c("Borreliella burgdorferi" = "B. burgdorferi",
                              "Borreliella garinii"     = "B. garinii",
                              "Borreliella afzelii"     = "B. afzelii †",
                              "Borreliella bavariensis" = "B. bavariensis")) +
  labs(x = "Species", y = "Mean cluster count per genome", fill = "Cluster category") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.ticks = element_line(color = "black", linewidth = 0.8),
        axis.ticks.length = unit(2.5, "pt"),
        axis.text.x = element_text(angle = 30, hjust = 1, face = "italic", size = 11, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12, color = "black", face = "bold"),
        axis.title.y = element_text(size = 12, color = "black"),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        plot.margin = margin(10, 10, 10, 10))
print(fig)
ggsave("results/cluster_category_stack.pdf", fig, width = 8, height = 5, dpi = 300)
cat("Done: cluster category distribution\n")