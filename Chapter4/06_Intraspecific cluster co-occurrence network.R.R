# ============================================================
# 06. Intraspecific plasmid cluster co-occurrence network
#     (B. burgdorferi, Jaccard >= 0.3 edges, walktrap modules; n = 91)
# Input : data/Supplementary_Table_S1.xls
#         data/strain_plasmid_status_master_corrected.csv
#         data/plasmid_clusters_cluster_nochrom.tsv (483 plasmid clusters, chromosomes removed)
# Output: results/bb_jaccard_network.pdf
#         results/bb_network_edges_nodes.xlsx (two sheets: edges / nodes)
# ============================================================
pkgs <- c("tidyverse", "readxl", "igraph", "writexl")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(tidyverse); library(readxl); library(igraph)

dir.create("results", showWarnings = FALSE)

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

# --- B. burgdorferi x cluster PAV matrix (91 plasmid-carrying genomes) ---
pav_bb <- clusters %>%
  inner_join(s1, by = "Genome_ID") %>%
  filter(species == "Borreliella burgdorferi",
         Genome_ID %in% plasmid_ids) %>%
  distinct(Genome_ID, rep_id) %>%
  mutate(present = 1) %>%
  pivot_wider(id_cols = Genome_ID, names_from = rep_id,
              values_from = present, values_fill = 0)
stopifnot(nrow(pav_bb) == 91)

# --- cross-check: clusters with >10% detection among the 91 genomes ---
bb_rates <- pav_bb %>%
  pivot_longer(-Genome_ID, names_to = "rep_id", values_to = "present") %>%
  group_by(rep_id) %>% summarise(rate = mean(present), .groups = "drop")
cat("Clusters with >10% detection in B. burgdorferi:", sum(bb_rates$rate > 0.10), "\n")

# --- Jaccard similarity between genomes ---
X <- as.matrix(pav_bb %>% select(-Genome_ID))
rownames(X) <- pav_bb$Genome_ID
inter_m <- tcrossprod(X)
sums <- rowSums(X)
union_m <- outer(sums, sums, "+") - inter_m
jaccard <- inter_m / union_m
jaccard[union_m == 0] <- 0

# --- edges: Jaccard >= 0.3, upper triangle (diagonal excluded) ---
edge_idx <- which(jaccard >= 0.3 & row(jaccard) < col(jaccard), arr.ind = TRUE)
edges <- tibble(genome_1 = rownames(jaccard)[edge_idx[, 1]],
                genome_2 = rownames(jaccard)[edge_idx[, 2]],
                jaccard  = jaccard[edge_idx])
cat("Edges (Jaccard >= 0.3):", nrow(edges), "\n")

# --- network + walktrap modules (all 91 nodes kept; isolates are meaningful) ---
g <- graph_from_data_frame(edges, directed = FALSE,
                           vertices = tibble(name = rownames(jaccard)))
wt <- cluster_walktrap(g)
V(g)$module <- membership(wt)
V(g)$degree <- degree(g)

cat("Walktrap modules:", length(wt), " modularity:", round(modularity(wt), 3), "\n")
cat("Isolated nodes (no edge >= 0.3):", sum(degree(g) == 0), "\n")
print(table(membership(wt)))

# --- output tables (single workbook, two sheets) ---
s_edges <- edges %>%
  mutate(module_1 = V(g)$module[match(genome_1, V(g)$name)],
         module_2 = V(g)$module[match(genome_2, V(g)$name)],
         same_module = module_1 == module_2)
s_nodes <- tibble(Genome_ID = V(g)$name,
                  n_clusters = rowSums(X)[V(g)$name],
                  degree = V(g)$degree, module = V(g)$module) %>%
  arrange(module, desc(degree))
writexl::write_xlsx(list(edges = s_edges, nodes = s_nodes),
                    "results/bb_network_edges_nodes.xlsx")

# --- plot: singletons merged in grey, legend lists multi-genome modules only ---
multi_mods <- as.character(sort(as.numeric(names(which(table(V(g)$module) >= 2)))))
V(g)$mod_plot <- ifelse(as.character(V(g)$module) %in% multi_mods,
                        as.character(V(g)$module), "Singleton")
pal <- scales::hue_pal()(length(multi_mods))
mod_colors <- c(setNames(pal, multi_mods), "Singleton" = "grey75")
n_iso <- sum(degree(g) == 0)

pdf("results/bb_jaccard_network.pdf", width = 8.5, height = 7)
set.seed(42)
plot(g, layout = layout_with_fr(g),
     vertex.color = mod_colors[V(g)$mod_plot],
     vertex.size = ifelse(V(g)$mod_plot == "Singleton", 3.5, 4 + V(g)$degree * 0.4),
     vertex.label = NA, vertex.frame.color = "black",
     edge.width = 0.5, edge.color = adjustcolor("grey40", alpha.f = 0.5),
     main = "B. burgdorferi plasmid-cluster Jaccard network (Jaccard >= 0.3, n = 91)")
legend("topright",
       legend = c(paste("Module", multi_mods, "(n =", table(V(g)$module)[multi_mods], ")"),
                  paste0("Singletons (n = ", n_iso, ")")),
       col = mod_colors[c(multi_mods, "Singleton")],
       pch = 16, pt.cex = 1.3, bty = "n", cex = 0.75)
dev.off()
cat("Done: intraspecific plasmid cluster co-occurrence network\n")