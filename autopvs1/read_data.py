#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# datetime: 2019/6/27 16:41

# import tabix
import os
import configparser
from pyfaidx import Fasta
from .pyhgvs.utils import read_transcripts
from .utils import read_morbidmap, read_pathogenic_site, \
    read_ba1_exception, read_pvs1_levels, create_bed_dict

BinPath = os.path.split(os.path.realpath(__file__))[0]

config = configparser.ConfigParser()
config.read(BinPath+'/config.ini')
omim_dict = read_morbidmap(BinPath+'/'+config['DEFAULT']['morbidmap'])
pathogenic_dict, pathogenic_dict2 = read_pathogenic_site(BinPath+'/'+config['DEFAULT']['pathogenic_ref'])
ba1_exception = read_ba1_exception(BinPath+'/'+config['DEFAULT']['ba1_exception'])
pvs1_levels = read_pvs1_levels(BinPath+'/'+config['DEFAULT']['pvs1levels'])

domain_bed = create_bed_dict(BinPath+'/'+config['DEFAULT']['domain'])
hotspot_bed = create_bed_dict(BinPath+'/'+config['DEFAULT']['hotspot'])
curated_region = create_bed_dict(BinPath+'/'+config['DEFAULT']['curated_region'])
exon_lof_frequent = create_bed_dict(BinPath+'/'+config['DEFAULT']['exon_lof_frequent'])


genome = Fasta(BinPath+'/'+config['DEFAULT']['ref'])
with open(BinPath+'/'+config['DEFAULT']['trans']) as gpefile:
    transcripts = read_transcripts(gpefile)

gene_trans = {}
trans_gene = {}
with open(BinPath+'/'+config['DEFAULT']['trans']) as f:
    for line in f:
        record = line.strip().split("\t")
        gene = record[12]
        trans = record[1]
        gene_trans[gene] = trans
        trans_gene[trans] = gene
