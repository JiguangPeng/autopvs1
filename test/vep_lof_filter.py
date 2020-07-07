#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# datetime: 2019/7/12 8:56

import sys

from collections import namedtuple

VAR = namedtuple('VAR', ('varid', 'gene', 'trans', 'canonical', 'pick', 'record'))

refgene = '../data/refGenePlus_20191020.gpe'
gene_trans = {}
trans_gene = {}

with open(refgene) as f:
    for line in f:
        record = line.strip().split("\t")
        gene = record[12]
        trans = record[1]
        gene_trans[gene] = trans
        trans_gene[trans] = gene


vep_lof_list = ['frameshift_variant',
                'stop_gained',
                'splice_donor_variant',
                'splice_acceptor_variant',
                'start_lost']


all_fh = open(sys.argv[1]+".all", 'w')
lof_fh = open(sys.argv[1]+".lof", 'w')

var_hash = {}

with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith("##"):
            continue
        elif line.startswith("#"):
            all_fh.write(line)
            lof_fh.write(line)
            continue
        r = line.strip().split("\t")
        r[1] = r[1] if r[1] != '-' else trans_gene.get(r[2], "-")
        r[1] = r[1] if r[1] != '-' else trans_gene.get(r[2].split(".")[0], "NA")
        var = VAR(r[0], r[1], r[2], r[3], r[4], r)
        if var.varid in var_hash:
            var_hash[var.varid].append(var)
        else:
            var_hash[var.varid] = [var]

for varid in var_hash:
    final_choose = []
    for var_anno in var_hash[varid]:
        if gene_trans.get(var_anno.gene) == var_anno.trans:
            final_choose.append(var_anno)

    if len(final_choose) == 0:
        for var_anno in var_hash[varid]:
            if (gene_trans.get(var_anno.gene, "na").split(".")[0] ==
                    var_anno.trans.split(".")[0]):
                final_choose.append(var_anno)

    final = ''
    if len(final_choose) > 1:
        for var_anno in final_choose:
            if var_anno.pick == '1':
                final = var_anno
        if not final:
            final = final_choose[0]
    elif len(final_choose) == 1:
        final = final_choose[0]
    else:
        final = var_hash[varid][0]

    all_fh.write("\t".join(final.record) + "\n")
    if any([i in final.record[5] for i in vep_lof_list]):
        lof_fh.write("\t".join(final.record) + "\n")





