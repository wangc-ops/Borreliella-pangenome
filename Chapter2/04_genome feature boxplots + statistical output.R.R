library(ggplot2)
library(dplyr)
library(readxl)
library(rstatix)

# -------------------- PATHS (edit as needed) --------------------
s1_path <- "C:/Users/YourName/Desktop/Supplementary_Table_S1.xls"
gff_dir <- "C:/YourPath/gff3_four_species"
fna_dir <- "C:/YourPath/fna"

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

s1_meta <- s1 %>%
  filter(species %in% sp_levels) %>%
  mutate(
    Genome_ID = sub("_genomic.*$", "", Genome_ID),
    species = factor(species, levels = sp_levels)
  ) %>%
  select(Genome_ID, species)

# -------------------- 2. Extract GFF3 features --------------------
extract_gff <- function(gid) {
  gff_path <- file.path(gff_dir, paste0(gid, "_genomic.gff3"))
  if (!file.exists(gff_path)) {
    return(c(CDS_N = NA, CDS_len = NA, MGE_N = NA, Annotated_N = NA, Hypothetical_N = NA))
  }
  df <- read.delim(gff_path, comment.char = "#", header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE,
                   col.names = c("seqid","source","type","start","end","score","strand","phase","attrs"))
  if (nrow(df) == 0) {
    return(c(CDS_N = NA, CDS_len = NA, MGE_N = NA, Annotated_N = NA, Hypothetical_N = NA))
  }
  is_cds <- df$type == "CDS"
  cds_n <- sum(is_cds, na.rm = TRUE)
  cds_len <- sum(df$end[is_cds] - df$start[is_cds] + 1, na.rm = TRUE)
  attrs <- df$attrs[is_cds]
  
  mge_pat <- "transposase|integrase|recombinase|resolvase|invertase|mobilization|conjugative|partition|toxin-antitoxin|insertion sequence|transposon"
  mge_n <- sum(grepl(mge_pat, attrs, ignore.case = TRUE))
  
  anno_pat <- "Dbxref=.*(COG:|KEGG:|GO:|EC:|Pfam:)"
  anno_n <- sum(grepl(anno_pat, attrs, ignore.case = TRUE))
  
  hypo_n <- sum(grepl("product=hypothetical protein", attrs, ignore.case = TRUE))
  
  c(CDS_N = cds_n, CDS_len = cds_len, MGE_N = mge_n, Annotated_N = anno_n, Hypothetical_N = hypo_n)
}

cat("Extracting GFF3...\n")
gff_res <- lapply(s1_meta$Genome_ID, extract_gff)
gff_df <- do.call(rbind, gff_res) %>% as.data.frame()
gff_df$Genome_ID <- s1_meta$Genome_ID

# -------------------- 3. Calculate FASTA stats --------------------
calc_fasta <- function(fp) {
  if (!file.exists(fp)) return(c(genome_size = NA, GC = NA, AT_rich = NA))
  lines <- readLines(fp)
  seq <- paste(lines[!grepl("^>", lines)], collapse = "")
  seq_clean <- gsub("[^A-Za-z]", "", seq)
  seq_upper <- toupper(seq_clean)
  n <- nchar(seq_upper)
  if (n == 0) return(c(genome_size = NA, GC = NA, AT_rich = NA))
  
  bytes <- charToRaw(seq_upper)
  gc_pct <- sum(bytes %in% charToRaw("GC")) / n * 100
  
  window <- 50
  at_rich <- NA_real_
  if (n >= window) {
    is_at <- as.integer((bytes == charToRaw("A")) | (bytes == charToRaw("T")))
    cs <- cumsum(is_at)
    ws <- cs[window:n] - c(0, cs[1:(n - window)])
    at_rich <- mean(ws / window > 0.90) * 100
  }
  c(genome_size = n, GC = gc_pct, AT_rich = at_rich)
}

cat("Calculating FASTA stats...\n")
fasta_res <- lapply(s1_meta$Genome_ID, function(gid) {
  cand <- file.path(fna_dir, paste0(gid, c("_genomic.fna", ".fna")))
  fp <- cand[file.exists(cand)][1]
  if (is.na(fp)) return(c(genome_size = NA, GC = NA, AT_rich = NA))
  calc_fasta(fp)
})
fasta_df <- do.call(rbind, fasta_res) %>% as.data.frame()
fasta_df$Genome_ID <- s1_meta$Genome_ID

# -------------------- 4. Merge & compute ratios --------------------
plot_data <- s1_meta %>%
  left_join(gff_df, by = "Genome_ID") %>%
  left_join(fasta_df, by = "Genome_ID") %>%
  mutate(
    Genome_Size_Mb = genome_size / 1e6,
    IGR_ratio = (genome_size - CDS_len) / genome_size * 100,
    MGE_ratio = MGE_N / CDS_N * 100,
    Annotated_ratio = Annotated_N / CDS_N * 100,
    Hypothetical_Ratio = Hypothetical_N / CDS_N * 100
  )

cat("Data ready:", nrow(plot_data), "genomes\n")

# -------------------- 5. Statistics: Welch + Games-Howell + CLD --------------------
get_cld <- function(data, y_var) {
  means <- data %>% 
    group_by(species) %>% 
    summarise(m = mean(!!sym(y_var), na.rm = TRUE), .groups = "drop") %>% 
    arrange(desc(m))
  
  f <- as.formula(paste(y_var, "~ species"))
  gh <- data %>% games_howell_test(f)
  
  letters <- character(length(sp_levels))
  names(letters) <- sp_levels
  cur <- "a"
  letters[as.character(means$species[1])] <- cur
  
  for (i in 2:nrow(means)) {
    sp <- as.character(means$species[i])
    assigned <- as.character(means$species[1:(i-1)])
    same <- FALSE
    for (prev in assigned) {
      pair <- gh %>% filter((group1 == sp & group2 == prev) | (group1 == prev & group2 == sp))
      if (nrow(pair) > 0 && pair$p.adj > 0.05) {
        letters[sp] <- letters[prev]
        same <- TRUE
        break
      }
    }
    if (!same) {
      cur <- intToUtf8(utf8ToInt(cur) + 1)
      letters[sp] <- cur
    }
  }
  letters[sp_levels]
}

main_metrics <- list(
  "Genome size (Mbp)" = "Genome_Size_Mb",
  "GC content (%)" = "GC",
  "Hypothetical protein (%)" = "Hypothetical_Ratio"
)

supp_metrics <- list(
  "IGR proportion (%)" = "IGR_ratio",
  "AT-rich low-complexity (%)" = "AT_rich",
  "MGE gene proportion (%)" = "MGE_ratio",
  "Functionally annotated gene (%)" = "Annotated_ratio"
)

all_metrics <- c(main_metrics, supp_metrics)

cld_list <- list()
for (m in names(all_metrics)) {
  col <- all_metrics[[m]]
  f <- as.formula(paste(col, "~ species"))
  aov_res <- plot_data %>% welch_anova_test(f)
  cld <- get_cld(plot_data, col)
  cld_list[[col]] <- cld
  cat(sprintf("\n%s: F(%.1f,%.1f) = %.2f, P = %.3g | CLD: %s\n", 
              m, aov_res$DFn, aov_res$DFd, aov_res$statistic, aov_res$p,
              paste(cld, collapse = ", ")))
}

# -------------------- 6. Plotting function (preview only) --------------------
draw_boxplot <- function(data, y_var, y_label, cld_letters) {
  cld_df <- data.frame(
    species = factor(sp_levels, levels = sp_levels),
    cld = cld_letters,
    stringsAsFactors = FALSE
  )
  y_max_vals <- data %>% group_by(species) %>% summarise(y_max = max(!!sym(y_var), na.rm = TRUE), .groups = "drop")
  cld_df <- cld_df %>% left_join(y_max_vals, by = "species")
  
  ggplot(data, aes(x = species, y = !!sym(y_var), fill = species)) +
    geom_boxplot(alpha = 0.85, outlier.shape = NA, width = 0.55, linewidth = 0.4, color = "black") +
    geom_jitter(width = 0.12, size = 1.8, alpha = 0.5, shape = 21, color = "black", stroke = 0.3) +
    scale_fill_manual(values = sp_colors) +
    scale_x_discrete(limits = sp_levels, labels = function(x) gsub("Borreliella ", "B. ", x)) +
    labs(x = "Species", y = y_label) +
    geom_text(data = cld_df, aes(x = species, y = y_max, label = cld), inherit.aes = FALSE,
              vjust = -0.8, size = 4.5, fontface = "bold", color = "black") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(2, "mm"),
      axis.text.x = element_text(angle = 30, hjust = 1, face = "italic", size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 11, face = "bold", color = "black"),
      axis.title.y = element_text(size = 11, face = "bold", color = "black"),
      legend.position = "none",
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
}

# -------------------- 7. Main text plots (Fig 2d-f) --------------------
cat("\n--- Main text: Fig 2d-f preview ---\n")
for (i in seq_along(main_metrics)) {
  m <- names(main_metrics)[i]
  col <- main_metrics[[i]]
  cld <- cld_list[[col]]
  cat("\n>>> Preview:", m, "| CLD:", paste(cld, collapse = " "), "\n")
  p <- draw_boxplot(plot_data, col, m, cld)
  print(p)
  if (i < length(main_metrics)) invisible(readline("Press Enter for next..."))
}

# -------------------- 8. Supplementary plots (Fig S3a-d) --------------------
cat("\n--- Supplementary: Fig S3a-d preview ---\n")
for (i in seq_along(supp_metrics)) {
  m <- names(supp_metrics)[i]
  col <- supp_metrics[[i]]
  cld <- cld_list[[col]]
  cat("\n>>> Preview:", m, "| CLD:", paste(cld, collapse = " "), "\n")
  p <- draw_boxplot(plot_data, col, m, cld)
  print(p)
  if (i < length(supp_metrics)) invisible(readline("Press Enter for next..."))
}

cat("\nAll previews done.\n")