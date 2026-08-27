# 去染色体 + 去表头行补丁（nochrom 文件的生成依据）
library(readr); library(dplyr)
f <- "plasmid_clusters_cluster_nochrom.tsv"   # 先由 stage1 文件剔除 median_len>500kb 的 cluster 生成
x <- read_tsv(f, col_names = c("rep_id", "member_id"), show_col_types = FALSE) %>%
  filter(rep_id != "rep_id")
write_tsv(x, f, col_names = FALSE)
cat("clusters:", n_distinct(x$rep_id), "(expected 483)\n")
