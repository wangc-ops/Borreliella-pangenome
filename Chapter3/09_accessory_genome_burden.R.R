# ============================================================
# Accessory genome burden per genome in plasmid-carrying genomes
# Accessory = non-core gene families (Soft-core + Shell + Cloud; Core >= 99%)
# Statistics: Welch's ANOVA + Games-Howell on three species;
# B. afzelii (n = 6) shown descriptively only
# ============================================================
library(data.table); library(ggplot2); library(dplyr)
library(rstatix); library(multcompView)

pav_file     <- "path/to/gene_presence_absence_roary.csv"
status_file  <- "path/to/plasmid_status_table.csv"   # Genome_ID, species, has_plasmid
out_dir      <- "path/to/output"

# --- Family frequency & accessory definition ---
pav <- fread(pav_file)
genome_cols <- names(pav)[15:ncol(pav)]
mat <- as.matrix(pav[, ..genome_cols])
mat_bin <- (mat != "" & mat != "0" & !is.na(mat)) * 1L
freq <- rowSums(mat_bin) / ncol(mat_bin)
accessory_fam <- freq < 0.99

burden <- data.table(Genome_ID = sub("_genomic$", "", genome_cols),
                     burden = colSums(mat_bin[accessory_fam, , drop = FALSE]))
status <- fread(status_file)
burden <- merge(burden, status[, .(Genome_ID, species, has_plasmid)], by = "Genome_ID")

# --- Plasmid-carrier subset ---
sp_levels <- c("Borreliella burgdorferi", "Borreliella garinii",
               "Borreliella afzelii", "Borreliella bavariensis")
sp_cols <- c("Borreliella burgdorferi" = "#1F618D",
             "Borreliella garinii"     = "#E74C3C",
             "Borreliella afzelii"     = "#1ABC9C",
             "Borreliella bavariensis" = "#5DADE2")
stat_sp <- sp_levels[c(1, 2, 4)]   # exclude B. afzelii (n = 6) from tests

sub <- burden[has_plasmid == TRUE]
sub[, species := factor(species, levels = sp_levels)]
ds <- sub[species %in% stat_sp]; ds[, species := droplevels(species)]

wt <- welch_anova_test(ds, burden ~ species)
gh <- games_howell_test(ds, burden ~ species)
cld <- multcompLetters(setNames(gh$p.adj, paste(gh$group1, gh$group2, sep = "-")),
                       threshold = 0.05)$Letters
print(wt); print(as.data.frame(gh)); print(cld)
write.csv(as.data.frame(gh), file.path(out_dir, "accessory_burden_games_howell.csv"),
          row.names = FALSE)

# --- Plot ---
cld_df <- data.frame(species = factor(sp_levels, levels = sp_levels),
                     cld = c(cld[stat_sp[1]], cld[stat_sp[2]], "†", cld[stat_sp[3]]))
y_max <- sub %>% group_by(species) %>%
  summarise(y_max = max(burden), .groups = "drop")
cld_df <- cld_df %>% left_join(y_max, by = "species")

p <- ggplot(sub, aes(x = species, y = burden, fill = species)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA, width = 0.55,
               linewidth = 0.4, color = "black") +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.5,
              shape = 21, color = "black", stroke = 0.3) +
  geom_text(data = cld_df, aes(x = species, y = y_max, label = cld),
            inherit.aes = FALSE, vjust = -0.8, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = sp_cols) +
  scale_x_discrete(labels = function(x) gsub("Borreliella ", "B. ", x)) +
  labs(x = NULL, y = "Accessory gene families per genome",
       caption = "Plasmid-carrying genomes only. † B. afzelii (n = 6) shown descriptively,\nexcluded from statistical tests (Welch's ANOVA with Games-Howell post hoc).") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1, face = "italic",
                                   size = 10, color = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        plot.caption = element_text(size = 7.5, hjust = 0, color = "grey30"))

ggsave(file.path(out_dir, "accessory_burden_plasmid_carriers.pdf"), p,
       width = 4.5, height = 4.8, dpi = 300)