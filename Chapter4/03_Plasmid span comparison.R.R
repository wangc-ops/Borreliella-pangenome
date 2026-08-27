# ============================================================
# 03. Plasmid span comparison (total plasmid length per genome, kb)
# Input : data/Supplementary_Table_S1.xls
#         data/strain_plasmid_status_master_corrected.csv
#         data/plasmid_contig_lengths_filtered.tsv (>5 kb plasmid contigs, no chromosomes)
# Output: results/plasmid_span.pdf
#         results/plasmid_span_summary.csv
# Note  : analyses conditional on 135 plasmid-carrying genomes;
#         B. afzelii (n=6) shown descriptively only (dagger), excluded from tests
# ============================================================
pkgs <- c("tidyverse", "readxl", "rstatix", "multcompView")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(tidyverse); library(readxl); library(rstatix); library(multcompView)

dir.create("results", showWarnings = FALSE)

sp_levels <- c("Borreliella burgdorferi", "Borreliella garinii",
               "Borreliella afzelii", "Borreliella bavariensis")
sp_test3  <- c("Borreliella burgdorferi", "Borreliella garinii",
               "Borreliella bavariensis")
sp_colors <- c("Borreliella burgdorferi" = "#1F618D",
               "Borreliella garinii"     = "#E74C3C",
               "Borreliella afzelii"     = "#1ABC9C",
               "Borreliella bavariensis" = "#5DADE2")

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

# --- contig lengths: auto-detect id column and length column ---
lens <- read_tsv("data/plasmid_contig_lengths_filtered.tsv", show_col_types = FALSE)
id_col  <- names(lens)[sapply(lens, function(x) is.character(x) &&
                                any(grepl("GCF_", x)))][1]
len_col <- names(lens)[sapply(lens, is.numeric)][1]
stopifnot(!is.na(id_col), !is.na(len_col), max(lens[[len_col]]) < 120000)

span_df <- lens %>%
  transmute(Genome_ID = sub("_genomic$", "", str_extract(.data[[id_col]], "^[^|]+")),
            len = .data[[len_col]]) %>%
  inner_join(s1, by = "Genome_ID") %>%
  filter(Genome_ID %in% plasmid_ids) %>%
  group_by(Genome_ID, species) %>%
  summarise(span_kb = sum(len) / 1000, .groups = "drop") %>%
  mutate(species = factor(species, levels = sp_levels))
stopifnot(nrow(span_df) == 135)

span_sum <- span_df %>%
  group_by(species) %>%
  summarise(n = n(), median = median(span_kb),
            q1 = quantile(span_kb, 0.25), q3 = quantile(span_kb, 0.75), .groups = "drop")
print(span_sum)

# --- Welch ANOVA + Games-Howell + CLD (bb/gar/bav only) ---
d <- span_df %>% filter(species %in% sp_test3) %>%
  mutate(species = factor(as.character(species), levels = sp_test3))
welch_res <- welch_anova_test(d, span_kb ~ species)
gh <- games_howell_test(d, span_kb ~ species)
print(welch_res); print(gh)

p_mat <- matrix(1, 3, 3, dimnames = list(sp_test3, sp_test3))
for (i in seq_len(nrow(gh))) {
  a <- as.character(gh$group1[i]); b <- as.character(gh$group2[i])
  p_mat[a, b] <- gh$p.adj[i]; p_mat[b, a] <- gh$p.adj[i]
}
raw_cld <- multcompLetters(p_mat, threshold = 0.05)$Letters
med_order <- d %>% group_by(species) %>% summarise(m = median(span_kb), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(species) %>% as.character()
uniq <- unique(unlist(strsplit(raw_cld[med_order], "")))
remap <- setNames(letters[seq_along(uniq)], uniq)
cld_letters <- sapply(raw_cld, function(x)
  paste(sort(remap[strsplit(x, "")[[1]]]), collapse = ""))
print(cld_letters)

write_csv(span_sum %>% mutate(cld = cld_letters[as.character(species)]),
          "results/plasmid_span_summary.csv")

# --- plot ---
cld_df <- span_df %>%
  filter(as.character(species) %in% names(cld_letters)) %>%
  group_by(species) %>%
  summarise(ypos = max(span_kb) + 0.06 * diff(range(span_df$span_kb)), .groups = "drop") %>%
  mutate(letter = cld_letters[as.character(species)])

fig <- ggplot(span_df, aes(x = species, y = span_kb, fill = species)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA, width = 0.6,
               color = "black", linewidth = 0.5) +
  geom_jitter(width = 0.12, size = 1.8, shape = 21,
              color = "black", stroke = 0.3, alpha = 0.7) +
  geom_text(data = cld_df, aes(x = species, y = ypos, label = letter),
            inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_fill_manual(values = sp_colors) +
  scale_x_discrete(labels = c("Borreliella burgdorferi" = "B. burgdorferi",
                              "Borreliella garinii"     = "B. garinii",
                              "Borreliella afzelii"     = "B. afzelii †",
                              "Borreliella bavariensis" = "B. bavariensis")) +
  labs(x = "Species", y = "Total plasmid span per genome (kb)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.ticks = element_line(color = "black", linewidth = 0.8),
        axis.ticks.length = unit(2.5, "pt"),
        axis.text.x = element_text(angle = 30, hjust = 1, face = "italic", size = 11, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12, color = "black", face = "bold"),
        axis.title.y = element_text(size = 12, color = "black"),
        legend.position = "none",
        plot.margin = margin(10, 10, 10, 10))
print(fig)
ggsave("results/plasmid_span.pdf", fig, width = 7, height = 5.5, dpi = 300)
cat("Done: plasmid span comparison\n")