#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/1/7 8:44 PM
    @Usage: python3 genePAV.py sample_gene.cov genecov_cutoff transcriptcov_cutoff cds_cutoff PAV.tsv
"""
# python3 gene_coverage_v0.py sample.bed chr.gff GenebodyCov_Cutoff atleastCdsCov_Cutoff SampleTag > sample_gene.cov


import sys
import pandas as pd

try:
    covfile = sys.argv[1]
    genecov_cutoff = float(sys.argv[2])
    transcriptcov_cutoff = float(sys.argv[3])
    cdscov_cutoff = float(sys.argv[4])
    pavfile = sys.argv[5]
except Exception as e:
    print('''
    python3 genePAV.py sample_gene.cov genecov_cutoff transcriptcov_cutoff cds_cutoff PAV.tsv
    e.g. python3 genePAV.py merged.cov 0 0 0.8 PAV.tsv
    
    # before this
    # python3 gene_coverage_v0.py sample.bed chr.gff GenebodyCov_Cutoff atleastCdsCov_Cutoff SampleTag > sample_gene.cov
    # cat *_gene.cov > merged.cov
    ''')

pav_dict = dict()

with open(covfile) as f:
    for line in f:
        temp = line.rstrip().split('\t')
        gene = temp[0]
        sample = temp[3]
        if sample not in pav_dict:
            pav_dict[sample] = dict()
        gene_cov = float(temp[4])
        transcript_cov = float(temp[5])
        cds_cov = float(temp[6])
        gene_exist = 0
        if gene_cov >= genecov_cutoff and transcript_cov >= transcriptcov_cutoff and cds_cov >= cdscov_cutoff:
            gene_exist = 1
        pav_dict[sample][gene] = gene_exist

pav_df = pd.DataFrame(pav_dict)
pav_df.fillna(value=0, inplace=True)
pav_df.to_csv(pavfile, sep='\t', float_format='%.0f')
