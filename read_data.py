#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# datetime: 2019/6/27 16:41

import os
import configparser
from pyfaidx import Fasta
from .pyhgvs.utils import read_transcripts
from .utils import read_morbidmap, read_pathogenic_site, read_pvs1_levels, create_bed_dict, read_gene_alias

BinPath = os.path.split(os.path.realpath(__file__))[0]

config = configparser.ConfigParser()
config.read(BinPath+'/config.ini')
pathogenic_dict, pathogenic_dict2 = read_pathogenic_site(BinPath+'/'+config['DEFAULT']['pathogenic_ref'])

domain_bed = create_bed_dict(BinPath+'/'+config['DEFAULT']['domain'])
hotspot_bed = create_bed_dict(BinPath+'/'+config['DEFAULT']['hotspot'])
curated_region = create_bed_dict(BinPath+'/'+config['DEFAULT']['curated_region'])
exon_lof_popmax = create_bed_dict(BinPath+'/'+config['DEFAULT']['exon_lof_popmax'])
pvs1_levels = read_pvs1_levels(BinPath+'/'+config['DEFAULT']['pvs1levels'])
gene_alias = read_gene_alias(BinPath+'/'+config['DEFAULT']['gene_alias'])

genome = Fasta(BinPath+'/'+config['DEFAULT']['ref'])
with open(BinPath+'/'+config['DEFAULT']['trans']) as gpefile:
    transcripts = read_transcripts(gpefile)

gene_trans = {}
trans_gene = {}
with open(BinPath+'/'+config['DEFAULT']['gene_trans']) as f:
    for line in f:
        record = line.strip().split("\t")
        gene,trans = record[0], record[1]
        gene_trans[gene] = trans
        trans_gene[trans] = gene
