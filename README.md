# Comparative genomics of Borreliella: analysis code and intermediate data

Analysis code and intermediate data supporting the manuscript:

**Pan-genomic analysis reveals plasmid-driven diversity and plasticity of Borreliella in Lyme disease pathogens**

## Data source

All 287 genome assemblies analysed in this study were obtained from the NCBI Genome database (January 2026) using the taxonomy query "Borreliella", excluding metagenome-assembled genomes (MAGs). Individual accession numbers, BioSample/BioProject identifiers, host and geographic metadata, and quality-control results (CheckM/QUAST/BUSCO) are provided in Table S1 of the supplementary material. No new sequencing data were generated in this study.

## Repository structure

- `Chapter1/` – genus-level pan-genome structure of the 287 genomes
- `Chapter2/` – phylogeny and genomic features of the four major pathogenic species
- `Chapter3/` – pan-genome tiers, rarefaction and Panstripe accumulation analyses
- `Chapter4/` – interspecific differentiation of the plasmidome
- `Chapter5/` – antigen repertoire, copy-number dosage and deposition-bias analyses
- `Chapter6/` – replicon-scale organisation and functional compartmentation
- `Supplementary_tables/` – Tables S1–S20 as submitted

Each folder contains the numbered R scripts and the input files needed to run them; scripts use relative paths.

## Software and versions

| Software | Version |
|---|---|
| R | 4.6.1 |
| vegan (R package) | 2.7-5 |
| rstatix (R package) | 1.0.0 |
| multcompView (R package) | 0.1-11 |
| CheckM | 1.2.2 |
| QUAST | 5.2.0 |
| BUSCO | 6.0.0 |
| Bakta | 1.11.1 (database 6.0) |
| FastANI | 1.34 |
| Panaroo | 1.3.0 |
| Panstripe | 0.4.0 |
| IQ-TREE | 3.0.1 |
| MMseqs2 | 18.8cc5c |
| HMMER | 3.4 |
| DIAMOND | 2.1.18 |
| VFDB | 2026 release |

Exact R package versions and platform details are recorded in `sessionInfo.txt`.

## Reproducing the analyses

1. Download the 287 genome assemblies from NCBI using the accession numbers in Table S1.
2. Run the command-line pipeline (genome annotation, pan-genome construction, phylogenetics, plasmid clustering, antigen assignment). The exact commands and parameters are given in the Materials and Methods of the manuscript (Sections 2.1–2.3); key commands are also collected in `pipeline_commands_en.md`.
3. Run the R scripts within each chapter folder in numerical order. Scripts use relative paths; place the required input files in the same folder before running.

## License

- Code: MIT License
- Data tables: CC0

## Contact

[Chong Wang], [wangchong@hainanu.edu.cn]
