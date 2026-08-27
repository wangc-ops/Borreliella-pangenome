Pipeline commands
Command-line workflow corresponding to Materials and Methods Sections 2.1–2.3 of the manuscript. All paths are placeholders; adjust to local directories. R analysis scripts are organized in the chapter folders; R package versions are recorded in sessionInfo.txt.

2.1 Dataset construction, quality control and annotation
CheckM v1.2.2 — completeness and contamination
checkm lineage_wf -t 16 -x fna genomes/ checkm_out/
QUAST v5.2.0 — N50 and contig counts
quast -o quast_out -t 16 genomes/*.fna
BUSCO v6.0.0 — gene-space completeness (bacterial lineage dataset)
busco -i genome.fna -m genome -l bacteria_odb12 -o busco_out -c 16
# Adjust the lineage dataset name to the locally installed BUSCO dataset
Bakta v1.11.1 — standardized re-annotation (Bakta database 6.0, default parameters)
bakta --db /path/to/bakta_db_6.0 --output bakta_out --prefix genome genome.fna
FastANI v1.34 — all-vs-all average nucleotide identity
fastANI --ql genome_list.txt --rl genome_list.txt \
        --fragLen 1000 --minFraction 0.2 -t 16 \
        -o fastani_all_vs_all.out

2.2 Pan-genome construction and phylogenetics
Panaroo v1.3.0 — pan-genome construction
# Run separately for the 287 genus-level genomes and the 241 four-species genomes, same parameters
panaroo -i bakta_gff3/*.gff3 -o panaroo_out \
        --clean-mode strict --threshold 0.7 --len_dif_percent 0.5 \
        -a core --aligner mafft -t 16
IQ-TREE v3.0.1 — maximum-likelihood phylogenetic trees [28]
# Genus level: concatenated alignment of 749 single-copy core genes (23 species)
iqtree3 -s core_genes_concat.fasta -m GTR+I+G -B 1000 -T 16

# Four species: core-genome SNP alignment (241 strains)
iqtree3 -s core_snp_alignment.fasta -m GTR+I+G -B 1000 -T 16
Panstripe v0.4.0 — accessory-genome accumulation dynamics
R package; see the Panstripe script in Chapter3/.

2.3 Plasmidome, antigen reservoir and plasmid-borne accessory genome
MMseqs2 v18.8cc5c — de novo clustering of plasmid contigs
# Preprocessing: extract contigs >5 kb per genome, exclude the longest contig
# (chromosome), and concatenate into a single fasta
mmseqs createdb plasmid_contigs.fna contigs_db
mmseqs easy-cluster contigs_db cluster_db mmseqs_tmp --min-seq-id 0.7 -c 0.8
mmseqs createtsv contigs_db contigs_db cluster_db plasmid_clusters_cluster.tsv
# Downstream (R/Python): clusters with representative sequences >500 kb were
# removed, yielding 483 plasmid clusters
HMMER v3.4 — PFam32 cp32/lp28 marker screening
hmmsearch --tblout pfam32_hits.tbl PFam32_cp32_lp28_markers.hmm panaroo_representative_proteins.faa
# Marker HMM library built from PFam32 partitioning proteins curated per Casjens et al. 2018 [4]
DIAMOND v2.1.18 — antigen family assignment
# Reference: 13 literature-curated antigen families [4, 32], downloaded from NCBI
diamond makedb --in antigen_reference.faa -d antigen_ref
diamond blastp -q panaroo_representative_proteins.faa -d antigen_ref \
      -o antigen_hits.tsv -e 1e-10 --query-cover 70 -k 1
# Best reciprocal hit assignment is performed downstream in R/Python (see the scripts in `Chapter5/`)
VFDB virulence factor screening
diamond makedb --in VFDB_2026.faa -d vfdb
diamond blastp -q panaroo_representative_proteins.faa -d vfdb \
      -o vfdb_hits.tsv -e 1e-10 --query-cover 70 -k 1

2.4 Statistical analysis
All statistical analyses were performed in R v4.6.1 (vegan 2.7-5, rstatix 1.0.0, multcompView 0.1-11, and others). See the numbered R scripts in each chapter folder and sessionInfo.txt for the full runtime environment.