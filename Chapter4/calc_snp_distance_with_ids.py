#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""带基因组ID的种内/种间核心基因SNP距离（numpy向量化）"""
import pandas as pd
import numpy as np
from Bio import AlignIO
import itertools

aln_path = "/mnt/c/W-数据分析/服务器加载数据/bakta_analysis/panaroo_four_species_complete/core_gene_alignment_filtered.aln"
s1_path  = "/mnt/c/Users/Lenovo/Desktop/过程文件及原始代码/附表/Supplementary_Table_S1.xls"
out_path = "/mnt/c/Users/Lenovo/Desktop/3.4_workflow/snp_distance_pairwise_with_ids.csv"

s1 = pd.read_excel(s1_path)
s1['Genome_ID'] = s1['Genome_ID'].astype(str).str.strip()
valid_species = ['Borreliella burgdorferi', 'Borreliella garinii',
                 'Borreliella afzelii', 'Borreliella bavariensis']
genome_to_species = {str(k).strip(): v for k, v in
                     zip(s1['Genome_ID'], s1['species'].astype(str).str.strip())
                     if v in valid_species}

print("Reading alignment...")
aln = AlignIO.read(aln_path, "fasta")
print(f"{len(aln)} seqs x {aln.get_alignment_length()} bp")

base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
def encode(seq):
    return np.array([base_map.get(b.upper(), -1) for b in seq], dtype=np.int8)

# 按物种组织：(gid, encoded_seq)
sp_data = {sp: {'ids': [], 'seqs': []} for sp in valid_species}
for record in aln:
    gid = record.id.replace("_genomic.fna", "").replace("_genomic", "")
    sp = genome_to_species.get(gid)
    if sp:
        sp_data[sp]['ids'].append(gid)
        sp_data[sp]['seqs'].append(encode(str(record.seq)))

results = []

def pair_dist(mat, ids, tag, sp_label):
    n = len(ids)
    for i, j in itertools.combinations(range(n), 2):
        a, b = mat[i], mat[j]
        valid = (a >= 0) & (b >= 0)
        nv = int(valid.sum())
        if nv == 0:
            continue
        nd = int((valid & (a != b)).sum())
        results.append({'genome_1': ids[i], 'genome_2': ids[j],
                        'species': sp_label, 'comparison': tag,
                        'snp_distance': nd / nv,
                        'diff_sites': nd, 'valid_sites': nv})

# 种内
for sp in valid_species:
    ids = sp_data[sp]['ids']
    mat = sp_data[sp]['seqs']
    print(f"Intra {sp}: {len(ids)} genomes, {len(ids)*(len(ids)-1)//2} pairs")
    pair_dist(mat, ids, 'Intra-species', sp)

# 种间（全配对，不抽样）
for sp1, sp2 in itertools.combinations(valid_species, 2):
    ids1, ids2 = sp_data[sp1]['ids'], sp_data[sp2]['ids']
    print(f"Inter {sp1} vs {sp2}: {len(ids1)*len(ids2)} pairs")
    for i in range(len(ids1)):
        a = sp_data[sp1]['seqs'][i]
        for j in range(len(ids2)):
            b = sp_data[sp2]['seqs'][j]
            valid = (a >= 0) & (b >= 0)
            nv = int(valid.sum())
            if nv == 0:
                continue
            nd = int((valid & (a != b)).sum())
            results.append({'genome_1': ids1[i], 'genome_2': ids2[j],
                            'species': f"{sp1} vs {sp2}",
                            'comparison': 'Inter-species',
                            'snp_distance': nd / nv,
                            'diff_sites': nd, 'valid_sites': nv})

res = pd.DataFrame(results)
res.to_csv(out_path, index=False)
print(f"\nSaved: {out_path} ({len(res)} pairs)")
print(res.groupby(['comparison', 'species'])['snp_distance'].median())
