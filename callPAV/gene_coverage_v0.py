#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/3/12 2:50 PM
    @Usage: python3 gene_coverage_v0.py sample.bed chr.gff GenebodyCov_Cutoff atleastCdsCov_Cutoff SampleTag > sample_gene.cov
"""
# bedtools genomecov -bga -split -ibam sample.bam > sample.bed

import sys


def string2dict(long_string, sep=';', eq='=', rm_quote=False):
    if rm_quote:
        long_string = long_string.replace('"', '').replace("'", '')
    long_string = long_string.replace('; ', ';')
    out_dict = dict()
    tmp = long_string.rstrip(sep).split(sep)
    for i in tmp:
        key, value = i.split(eq)
        out_dict[key] = value
    return out_dict


def union_length(interval_list):
    if interval_list == []:
        return 0
    interval_list = sorted(interval_list)
    list_length = len(interval_list)
    count_length = 0
    for i in range(1, list_length):
        if interval_list[i-1][1] > interval_list[i][1]:
            interval_list[i] = (interval_list[i][0], interval_list[i-1][1])
        count_length += interval_list[i-1][1] - interval_list[i-1][0] + 1 - max(0, interval_list[i-1][1] - interval_list[i][0]+1)
    count_length += interval_list[list_length-1][1] - interval_list[list_length-1][0] + 1
    return count_length


def intersect_interval(interval1, interval2):
    interval1 = sorted(interval1)
    interval2 = sorted(interval2)
    i, j = 0, 0
    res = []
    while i < len(interval1) and j < len(interval2):
        a1, a2 = interval1[i][0], interval1[i][1]
        b1, b2 = interval2[j][0], interval2[j][1]
        if b2 >= a1 and a2 >= b1:
            res.append([max(a1, b1), min(a2, b2)])
        if b2 < a2:
            j += 1
        else:
            i += 1
    return res


def readgff(gff, usedchr='WGC'):
    gene_dict = dict()
    with open(gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            # if not line.startswith(usedchr):
            #    continue
            temp = line.rstrip().split()
            chrn = temp[0]
            annotype = temp[2]
            details = temp[8]
            atrribute = string2dict(details)
            # geneinform = # temp[8].split(';')
            idname = atrribute['ID']  # geneinform[0].split('=')[1]
            start_pos = int(temp[3])
            end_pos = int(temp[4])
            if annotype == "gene":
                # start, end, cds
                geneidname = idname
                if 'gene_name' in atrribute:
                    genename = atrribute['gene_name']
                else:
                    genename = geneidname
                gene_dict[geneidname] = [chrn, start_pos, end_pos, [], [], genename]
            if annotype == "mRNA" or annotype == "transcript":
                geneidname = atrribute['Parent']  # idname.split(':')[0].split('-')[0]
                gene_dict[geneidname][4] = [(start_pos, end_pos)]
                    # if not geneidname in gene_dict:
                    #    gene_dict[geneidname] =  [chrn,start_pos,end_pos,[],[trans_s,trans_e]]

            if annotype == "CDS":
                # if temp[1] == "maker":
                # geneidname = idname.split(':')[0].split('-')[0]
                    # if geneidname in gene_dict:
                     #    gene_dict[geneidname][4] = geneidname
                # else:
                    # geneidname = temp[8].split('gene_id=')[1].split(";")[0]
                #gene_dict[geneidname][5] = geneidname  # temp[8].split('gene_name=')[1].split(";")[0]
                gene_dict[geneidname][3] += [(start_pos, end_pos)]
    return gene_dict


if __name__ == "__main__":
    try:
        bedtoolsout = sys.argv[1]
        gfffile = sys.argv[2]
        genecov_cutoff = float(sys.argv[3])
        cdscov_cutoff = float(sys.argv[4])
        atleast_cov = 1
        sample = sys.argv[5]
    except Exception as e:
        print('''python3 gene_coverage_v0.py sample.bed chr.gff GenebodyCov_Cutoff atleastCdsCov_Cutoff SampleTag > sample_gene.cov
        # bedtools genomecov -bga -split -ibam sample.bam > sample.bed
        ''')

    gffOBJ = readgff(gfffile)
    #print(gffOBJ)

    chrn_covregion = dict()
    with open(bedtoolsout) as f:
        for line in f:
            temp = line.rstrip().split('\t')
            cov = float(temp[3]) # int error if 1.1e6

            chrn = temp[0]
            start_pos = int(temp[1])+1
            end_pos = int(temp[2])
            if chrn not in chrn_covregion:
                chrn_covregion[chrn] = []
            if cov < atleast_cov:
                continue

            # update chr_covregion
            try:
                if start_pos == chrn_covregion[chrn][-1][1]+1:
                    chrn_covregion[chrn][-1][1] = end_pos
                else:
                    chrn_covregion[chrn] += [[start_pos, end_pos]]
            except Exception as e:
                chrn_covregion[chrn] += [[start_pos, end_pos]]
    
    print("#Gene_id\tChr\tGene_name\tSample\tGene_body_coverage\tTranscript_body_coverage\tCDS_body_coverage")
    for key in gffOBJ:
        chrn = gffOBJ[key][0]
        if not chrn in chrn_covregion:
            continue
        if gffOBJ[key][3] == []:
            continue
        gene_interval = [(gffOBJ[key][1], gffOBJ[key][2])]
        gene_cov_interval = intersect_interval(gene_interval, chrn_covregion[chrn])
        gene_cov = union_length(gene_cov_interval) / union_length(gene_interval)
        transcript_interval = gffOBJ[key][4]
        if transcript_interval == gene_interval:
            transcript_cov = gene_cov
        else:
            transcript_cov_interval = intersect_interval(transcript_interval, chrn_covregion[chrn])
            transcript_length = union_length(transcript_interval)
            if transcript_length == 0:
                transcript_cov = 0
            else:
                transcript_cov = union_length(transcript_cov_interval) / transcript_length
        #print("Gene: "+str(gene_interval))
        #print("CDS: "+str(gffOBJ[key][3]))
        #print("Map: "+str(chr_covregion[chr]))
        cds_interval = gffOBJ[key][3]
        cds_cov_interval = intersect_interval(cds_interval, chrn_covregion[chrn])
        cds_length = union_length(gffOBJ[key][3])
        if cds_length == 0:
            cds_cov = 0
        else:
            cds_cov = union_length(cds_cov_interval) / cds_length
        print('\t'.join([key,gffOBJ[key][0], gffOBJ[key][5], sample, str(gene_cov), str(transcript_cov), str(cds_cov)]))
        #print(key + "\t" + gffOBJ[key][0] + "\t" + gffOBJ[key][4] + "\t" + sample + "\t" + str(gene_cov)+ "\t" + str(transcript_cov) + "\t" + str(cds_cov))
