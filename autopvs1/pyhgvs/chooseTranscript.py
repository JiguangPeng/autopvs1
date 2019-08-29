#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# datetime: 2019/7/16 15:55

clinvar = open('Gene.Transcript.Clinvar')

clinvar_trans_gene = {}
for line in clinvar:
    record = line.strip().split("\t")
    clinvar_trans_gene[record[1]] = record[0]


