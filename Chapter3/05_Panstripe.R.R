# ============================================================================
# Panstripe pangenome evolutionary rate analysis

# Purpose: Infer gene gain/loss rates for four Borreliella pathogenic species
#          based on phylogenetic branch lengths
# 
# Dependencies: R >= 4.0, tidyverse, ape, panstripe, patchwork, readxl
# 
# Input files:
#   1. gene_presence_absence.csv  -- Panaroo PAV matrix output
#   2. tree.nwk                   -- Rooted phylogenetic tree (Newick, with branch lengths)
#   3. metadata.xlsx              -- Sample metadata (Genome_ID, species)
# 
# Important notes:
#   - The panstripe pa argument expects rows = genomes, columns = genes in the
#     source code, contrary to the package documentation (rows = genes). This
#     code performs the necessary transpose.
#   - Zero-length branches are jittered with 1e-7 to avoid GLM divergence.
#   - Some species require fallback to quasipoisson or gaussian families.
# ============================================================================

library(tidyverse)
library(ape)
library(panstripe)
library(patchwork)

# --- 1. Path configuration (modify as needed) ---
pav_path      <- "data/panaroo/gene_presence_absence.csv"
tree_path     <- "data/phylogeny/tree_four_species.nwk"
metadata_path <- "data/metadata.xlsx"
out_dir       <- "output/panstripe"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- 2. Read sample metadata ---
# Required columns: Genome_ID (genome identifier), species (species name)
meta <- readxl::read_excel(metadata_path) %>%
  select(Genome_ID, species) %>%
  mutate(
    Genome_ID = sub("_genomic$", "", Genome_ID),
    species = gsub("Borrelia", "Borreliella", species),
    species = factor(species, levels = c(
      "Borreliella burgdorferi", "Borreliella garinii",
      "Borreliella afzelii", "Borreliella bavariensis"
    ))
  )

# --- 3. Read PAV matrix and convert to 0/1 ---
# Panaroo format: rows = gene families, columns = genomes, values = gene sequence or empty
pav_raw <- read_csv(pav_path, show_col_types = FALSE)

meta_cols <- c("Gene", "Non-unique Gene name", "Annotation",
               "No. isolates", "No. sequences", "Avg group size nuc",
               "Genome Fragment", "Order within Fragment",
               "Accessory Fragment", "Accessory Order with Fragment",
               "Scaffold", "Position within scaffold",
               "Minimum Start", "Maximum End",
               "Average Start", "Average End",
               "Strand", "Annotation ID", "Pangenome fraction",
               "Maximun Group Size", "Minimum Group Size",
               "Average Group Size", "Presence/Absence Rtab")

genome_cols <- setdiff(names(pav_raw), meta_cols)
genome_cols_clean <- sub("_genomic$", "", genome_cols)

pav_mat <- pav_raw %>%
  select(all_of(genome_cols)) %>%
  mutate(across(everything(), ~ as.integer(!is.na(.) & . != ""))) %>%
  as.matrix()
rownames(pav_mat) <- pav_raw$Gene
colnames(pav_mat) <- genome_cols_clean

# --- 4. Read and process phylogenetic tree ---
tree <- read.tree(tree_path)

# Root using the longest branch as outgroup
tip_edges <- which(tree$edge[, 2] <= length(tree$tip.label))
longest_tip <- tree$tip.label[tree$edge[tip_edges, 2][which.max(tree$edge.length[tip_edges])]]
tree <- root(tree, outgroup = longest_tip, resolve.root = TRUE)
tree <- ladderize(tree)

# Clean tip labels to match PAV column names
tree$tip.label <- sub("_genomic$", "", tree$tip.label)
tree <- keep.tip(tree, intersect(tree$tip.label, colnames(pav_mat)))

# --- 5. Globally align PAV matrix to tree ---
pav_aligned <- pav_mat[, tree$tip.label, drop = FALSE]
pav_aligned[is.na(pav_aligned)] <- 0

# --- 6. Per-species Panstripe analysis function ---
run_panstripe_species <- function(sp, pav_aligned, tree, meta) {
  
  gids <- tree$tip.label[tree$tip.label %in% meta$Genome_ID[meta$species == sp]]
  cat("\n", sp, ": n =", length(gids), "\n")
  
  if (length(gids) < 5) return(NULL)
  
  tree_sp <- keep.tip(tree, gids)
  if (!is.binary(tree_sp)) tree_sp <- multi2di(tree_sp, random = FALSE)
  
  # Jitter zero-length branches
  zero_len <- tree_sp$edge.length == 0
  if (any(zero_len)) tree_sp$edge.length[zero_len] <- 1e-7
  
  # Extract PAV and transpose (critical: source code expects rows = genomes, cols = genes)
  col_idx <- match(tree_sp$tip.label, colnames(pav_aligned))
  pa_raw <- pav_aligned[, col_idx, drop = FALSE]
  pa_sp <- t(pa_raw)
  
  stopifnot(nrow(pa_sp) == length(tree_sp$tip.label))
  stopifnot(all(rownames(pa_sp) == tree_sp$tip.label))
  
  # Four-level fallback fitting strategy
  fit <- NULL; method_used <- ""
  
  # L1: Default Poisson
  fit <- tryCatch(panstripe(pa = pa_sp, tree = tree_sp, nboot = 1000), error = function(e) NULL)
  if (!is.null(fit)) method_used <- "default"
  
  # L2: Quasipoisson (handles overdispersion)
  if (is.null(fit)) {
    fit <- tryCatch(panstripe(pa = pa_sp, tree = tree_sp, nboot = 1000,
                              fit_method = "glm", family = "quasipoisson"), error = function(e) NULL)
    if (!is.null(fit)) method_used <- "glm_quasipoisson"
  }
  
  # L3: Gaussian (most stable)
  if (is.null(fit)) {
    fit <- tryCatch(panstripe(pa = pa_sp, tree = tree_sp, nboot = 1000,
                              fit_method = "glm", family = "gaussian"), error = function(e) NULL)
    if (!is.null(fit)) method_used <- "glm_gaussian"
  }
  
  # L4: No bootstrap (point estimates only)
  if (is.null(fit)) {
    fit <- tryCatch(panstripe(pa = pa_sp, tree = tree_sp, nboot = 0), error = function(e) NULL)
    if (!is.null(fit)) method_used <- "default_noboot"
  }
  
  if (is.null(fit)) return(NULL)
  
  cat("  Success:", method_used, "\n")
  
  smry <- fit$summary; data_df <- fit$data
  
  get_term <- function(tn) ifelse(tn %in% smry$term, smry$estimate[smry$term == tn], NA)
  get_p    <- function(tn) ifelse(tn %in% smry$term, smry$p.value[smry$term == tn], NA)
  
  list(
    species = sp, fit = fit, method = method_used, n_genomes = length(gids),
    core_coef = get_term("core"), core_p = get_p("core"),
    tip_coef = get_term("istip"), tip_p = get_p("istip"),
    depth_coef = get_term("depth"),
    standardized_rate = sum(data_df$acc) / sum(data_df$core)
  )
}

# --- 7. Run analysis for all four species ---
species_list <- levels(meta$species)
sp_results <- list()
for (sp in species_list) sp_results[[sp]] <- run_panstripe_species(sp, pav_aligned, tree, meta)
sp_results <- compact(sp_results)

# --- 8. Save results ---
res_df <- map_dfr(sp_results, ~ data.frame(
  species = .x$species, method = .x$method, n_genomes = .x$n_genomes,
  core_coef = .x$core_coef, core_p = .x$core_p,
  tip_coef = .x$tip_coef, tip_p = .x$tip_p,
  depth_coef = .x$depth_coef, standardized_rate = .x$standardized_rate
))

write_csv(res_df, file.path(out_dir, "panstripe_four_species_summary.csv"))

# --- 9. Generate publication-quality visualizations ---
sp_colors <- c(
  "Borreliella burgdorferi" = "#1F618D",
  "Borreliella garinii" = "#E74C3C",
  "Borreliella afzelii" = "#1ABC9C",
  "Borreliella bavariensis" = "#5DADE2"
)

# Parameter estimate plots (per species)
plot_list <- map(sp_results, ~ {
  plot_pangenome_params(.x$fit) + 
    ggtitle(.x$species) +
    theme(plot.title = element_text(face = "italic", size = 12))
})

walk2(plot_list, sp_results, ~ {
  ggsave(file.path(out_dir, paste0("panstripe_params_", .y$species, ".png")),
         .x, width = 8, height = 5, dpi = 300)
})

# Combined parameter plot (2x2 layout)
if (length(plot_list) >= 2) {
  combined <- wrap_plots(plot_list, ncol = 2) + 
    plot_annotation(title = "Panstripe parameter estimates across four *Borreliella* species")
  ggsave(file.path(out_dir, "panstripe_params_combined.png"), combined,
         width = 14, height = 10, dpi = 300)
}

# Cumulative curve plots
walk2(sp_results, sp_results, ~ {
  p_cum <- plot_pangenome_cumulative(.x$fit) + 
    ggtitle(.x$species) +
    theme(plot.title = element_text(face = "italic"))
  ggsave(file.path(out_dir, paste0("panstripe_cumulative_", .x$species, ".png")),
         p_cum, width = 7, height = 5, dpi = 300)
})

# Gain/loss tree plots
walk2(sp_results, sp_results, ~ {
  p_gl <- plot_gain_loss(.x$fit) + 
    ggtitle(.x$species) +
    theme(plot.title = element_text(face = "italic"))
  ggsave(file.path(out_dir, paste0("panstripe_gainloss_", .x$species, ".png")),
         p_gl, width = 12, height = 10, dpi = 300)
})

cat("\n=== Complete. Output directory:", out_dir, "===\n")