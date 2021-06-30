#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# created: 2019/6/27 19:13

import os
import re
import sys
import random
import string
from collections import namedtuple

from .pvs1 import PVS1
from .cnv import PVS1CNV, CNVRecord
from .read_data import trans_gene, gene_trans, gene_alias, vep_cache
from .read_data import transcripts_hg19, transcripts_hg38, genome_hg19, genome_hg38
from .utils import vep2vcf, get_transcript, vep_consequence_trans, VCFRecord


lof_type = ['frameshift', 'nonsense', 'splice-5', 'splice-3', 'init-loss']
vep_lof_list = ['frameshift', 'stop_gained', 'splice_donor', 'splice_acceptor', 'start_lost']
VAR = namedtuple('VAR', ('varid', 'gene', 'trans', 'canonical', 'pick', 'record'))


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


class AutoPVS1:
    """Run AutoPVS1"""

    def __init__(self, vcfrecord, genome_version, user_trans=None):
        self.vcfrecord = self.get_vcfrecord(vcfrecord)
        self.chrom = self.vcfrecord.chrom
        self.pos = self.vcfrecord.pos
        self.ref = self.vcfrecord.ref
        self.alt = self.vcfrecord.alt
        self.vcfrecord = VCFRecord(self.chrom, self.pos, self.ref, self.alt)
        self.user_trans = user_trans
        
        if genome_version in ['hg19', 'GRCh37']:
            self.genome_version = 'hg19'
            self.vep_assembly = 'GRCh37'
        elif genome_version in ['hg38', 'GRCh38']:
            self.genome_version = 'hg38'
            self.vep_assembly = 'GRCh38'
        else:
            raise ValueError("Genome version must be hg19/GRCh37 or hg38/GRCh38.")
            self.genome_version = 'hg38'
            self.vep_assembly = 'GRCh38'

        self.id = id_generator()
        self.vep_input = '/tmp/vep.{0}.vcf'.format(self.id)
        self.vep_output = '/tmp/vep.{0}.tab'.format(self.id)

        self.vep_variation = 'na'
        self.vep_symbol = 'na'
        self.vep_trans = 'na'
        self.vep_canonical = 'na'
        self.vep_pick = 'na'
        self.vep_consequence = 'na'
        self.hgvs_c = 'na'
        self.hgvs_p = 'na'
        self.hgvs_g = 'na'
        self.vep_exon = 'na'
        self.vep_intron = 'na'
        self.vep_run()
        self.vep_filter()

        if self.genome_version == 'hg19':
            self.transcript = get_transcript(self.vep_trans, transcripts_hg19)
        else:
            self.transcript = get_transcript(self.vep_trans, transcripts_hg38)

        self.consequence = vep_consequence_trans(self.vep_consequence)
        self.islof = self.consequence in lof_type
        if self.islof:
            self.pvs1 = self.run_pvs1()

        os.system('rm ' + self.vep_input + ' ' + self.vep_output)

    @staticmethod
    def get_vcfrecord(vcf):
        if all(hasattr(vcf, attr) for attr in ['chrom', 'pos', 'ref', 'alt']):
            return vcf
        elif re.match(r'^((?:chr)?[\dXYMT]{1,2})-(\d+)-([ATCG]+)-([ATCG]+)$', vcf, re.I):
            v = vcf.split("-")
            v[0] = re.sub('chr', '', v[0], flags=re.I)
            return VCFRecord(v[0], int(v[1]), v[2], v[3])
        else:
            raise TypeError("Wrong VCF Record")

    def vep_run(self):
        print(self.chrom, self.pos, '.', self.ref, self.alt, '.', 'PASS', ',',
              sep="\t", file=open(self.vep_input, 'w'))
        vepcommand = '''
            vep --offline --refseq --use_given_ref \
            --dir_cache ''' + vep_cache + ''' \
            --species "homo_sapiens" \
            --assembly ''' + self.vep_assembly + ''' \
            --fork 4 \
            --canonical \
            --flag_pick \
            --hgvs --hgvsg --symbol \
            --distance 500 \
            --exclude_predicted \
            --lookup_ref \
            --numbers \
            --force \
            --input_file ''' + self.vep_input + '''\
            --output_file ''' + self.vep_output + ''' --no_stats \
            --tab --fields "Uploaded_variation,SYMBOL,Feature,CANONICAL,PICK,Consequence,HGVSc,HGVSp,HGVSg,EXON,INTRON"
        '''
        os.system(vepcommand)

    def vep_filter(self):
        var_dict = {}
        with open(self.vep_output) as fh:
            for line in fh:
                if line.startswith("##"):
                    continue
                elif line.startswith("#"):
                    header = line.strip().split("\t")
                    header[0] = header[0].replace("#", "")
                    continue
                records = line.strip().split("\t")
                info = dict(zip(header, records))
                if info['SYMBOL'] == '-':
                    info['SYMBOL'] = trans_gene.get(info['Feature'], "NA")
                if info['SYMBOL'] == 'NA':
                    info['SYMBOL'] = trans_gene.get(info['Feature'].split(".")[0], "NA")
                if info['SYMBOL'] in gene_alias:
                    info['SYMBOL'] = gene_alias.get(info['SYMBOL'])
                var = VAR(info['Uploaded_variation'], info['SYMBOL'],
                          info['Feature'], info['CANONICAL'], info['PICK'], info)
                if var.varid in var_dict:
                    var_dict[var.varid].append(var)
                else:
                    var_dict[var.varid] = [var]

        for varid in var_dict:
            final_choose = []
            for var_anno in var_dict[varid]:
                if (gene_trans.get(var_anno.gene) == var_anno.trans or
                    gene_trans.get(var_anno.gene, "na").split(".")[0] == var_anno.trans.split(".")[0]):
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
                final = var_dict[varid][0]

            if self.user_trans:
                for var_anno in var_dict[varid]:
                    if var_anno.trans.split(".")[0] == self.user_trans.split(".")[0]:
                        final = var_anno
                        break

            self.vep_variation = final.record['Uploaded_variation']
            self.vep_symbol = final.record['SYMBOL']
            self.vep_trans = final.record['Feature']
            self.vep_canonical = final.record['CANONICAL']
            self.vep_pick = final.record['PICK']
            self.vep_consequence = final.record['Consequence'].replace('_variant', '')
            self.hgvs_c = final.record['HGVSc']
            self.hgvs_p = final.record['HGVSp'].replace('%3D', '=')
            self.hgvs_g = final.record['HGVSg']
            self.vep_exon = final.record['EXON']
            self.vep_intron = final.record['INTRON']

    def run_pvs1(self):
        pvs1 = PVS1(self.vcfrecord, self.consequence, self.hgvs_c, 
                    self.hgvs_p, self.transcript, self.genome_version)
        return pvs1


class AutoPVS1CNV:
    """Auto PVS1 for CNV"""
    cnvpattern1 = re.compile(r'(del|dup|tdup|ntdup)\((chr)?(\d+|X|Y):(\d+)-(\d+)\)', re.I)
    cnvpattern2 = re.compile(r'(chr)?(\d+|X|Y)-(\d+)-(\d+)-(del|dup|tdup|ntdup)', re.I)

    def __init__(self, cnv, genome_version, user_trans=None):
        cnvmatch1 = self.cnvpattern1.match(cnv)
        cnvmatch2 = self.cnvpattern2.match(cnv)
        if cnvmatch1:
            self.cnvtype = cnvmatch1.group(1).upper()
            self.chrom = cnvmatch1.group(3)
            self.start = int(cnvmatch1.group(4))
            self.end = int(cnvmatch1.group(5))
        elif cnvmatch2:
            self.chrom = cnvmatch2.group(2)
            self.start = int(cnvmatch2.group(3))
            self.end = int(cnvmatch2.group(4))
            self.cnvtype = cnvmatch2.group(5).upper()
        else:
            self.cnvtype = None

        if genome_version in ['hg19', 'GRCh37']:
            self.genome_version = 'hg19'
            self.vep_assembly = 'GRCh37'
        elif genome_version in ['hg38', 'GRCh38']:
            self.genome_version = 'hg38'
            self.vep_assembly = 'GRCh38'
        else:
            raise ValueError("wrong genome version, use hg38/GRCh38 by default")
            self.genome_version = 'hg38'
            self.vep_assembly = 'GRCh38'

        if self.cnvtype:
            self.cnvvariant = '-'.join([self.chrom, str(self.start), str(self.end), self.cnvtype])
            self.cnvrecord = CNVRecord(self.chrom, self.start, self.end, self.cnvtype)
            self.user_trans = user_trans
            self.vep_variation = 'na'
            self.vep_symbol = 'na'
            self.vep_trans = 'na'
            self.vep_canonical = 'na'
            self.vep_pick = 'na'
            self.vep_consequence = 'na'
            self.vep_exon = 'na'
            self.vep_intron = 'na'
            self.vep_cds_pos = 'na'

            self.id = id_generator()
            self.vep_input = '/tmp/vep.{0}.vcf'.format(self.id)
            self.vep_output = '/tmp/vep.{0}.tab'.format(self.id)
            self.vep_run()
            self.vep_filter()
            if self.genome_version == 'hg19':
                self.transcript = get_transcript(self.vep_trans, transcripts_hg19)
            else:
                self.transcript = get_transcript(self.vep_trans, transcripts_hg38)
            self.pvs1 = self.runpvs1()

        os.system('rm ' + self.vep_input + ' ' + self.vep_output)

    def vep_run(self):
        print(self.chrom, self.start, self.end, self.cnvtype, file=open(self.vep_input, 'w'))
        vepcommand = '''
            vep --offline --refseq --use_given_ref \
                --dir_cache ''' + vep_cache + ''' \
                --species "homo_sapiens" \
                --assembly ''' + self.vep_assembly + ''' \
                --fork 1 \
                --hgvs --hgvsg --canonical --symbol \
                --distance 0 \
                --exclude_predicted \
                --flag_pick \
                --lookup_ref \
                --force \
                --input_file ''' + self.vep_input + ''' \
                --output_file ''' + self.vep_output + ''' --no_stats \
                --numbers \
                --tab --fields "Uploaded_variation,SYMBOL,Feature,CANONICAL,PICK,EXON,INTRON,Consequence,CDS_position"
        '''
        os.system(vepcommand)

    def vep_filter(self):
        var_dict = {}
        with open(self.vep_output) as fh:
            for line in fh:
                if line.startswith("##"):
                    continue
                elif line.startswith("#"):
                    header = line.strip().split("\t")
                    header[0] = header[0].replace("#", "")
                    continue
                records = line.strip().split("\t")
                info = dict(zip(header, records))
                if info['SYMBOL'] == '-':
                    info['SYMBOL'] = trans_gene.get(info['Feature'], "-")
                if info['SYMBOL'] == '-':
                    info['SYMBOL'] = trans_gene.get(info['Feature'].split(".")[0], "NA")
                var = VAR(info['Uploaded_variation'], info['SYMBOL'],
                          info['Feature'], info['CANONICAL'], info['PICK'], info)
                if var.varid in var_dict:
                    var_dict[var.varid].append(var)
                else:
                    var_dict[var.varid] = [var]

        for varid in var_dict:
            final_choose = []
            for var_anno in var_dict[varid]:
                if gene_trans.get(var_anno.gene) == var_anno.trans:
                    final_choose.append(var_anno)

            if len(final_choose) == 0:
                for var_anno in var_dict[varid]:
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
                final = var_dict[varid][0]

            if self.user_trans:
                for var_anno in var_dict[varid]:
                    if var_anno.trans.split(".")[0] == self.user_trans.split(".")[0]:
                        final = var_anno
                        break

            self.vep_variation = final.record['Uploaded_variation']
            self.vep_symbol = final.record['SYMBOL']
            self.vep_trans = final.record['Feature']
            self.vep_canonical = final.record['CANONICAL']
            self.vep_pick = final.record['PICK']
            self.vep_consequence = final.record['Consequence'].replace('_variant', '')
            self.vep_exon = final.record['EXON']
            self.vep_intron = final.record['INTRON']
            self.vep_cds_pos = final.record['CDS_position']

    def runpvs1(self):
        pvs1 = PVS1CNV(self.cnvrecord,
                       self.vep_consequence,
                       self.transcript,
                       self.genome_version)
        return pvs1


def main():
    genome_version = sys.argv[1]
    anno_fh = open(sys.argv[2])
    header = list()
    for line in anno_fh:
        if line.strip().startswith("##"):
            continue
        elif line.strip().startswith("#"):
            header = line.strip().split("\t")
            header[0] = header[0].replace("#", "")
            continue
        records = line.strip().split("\t")
        if len(records) == len(header):
            info = dict(zip(header, records))
        else:
            raise Exception("Inconsistent length for line and header!")

        if genome_version == 'hg19':
            vcfrecord = vep2vcf(info['Uploaded_variation'], genome_hg19)
            transcript = get_transcript(info['Feature'], transcripts_hg19)
        else:
            vcfrecord = vep2vcf(info['Uploaded_variation'], genome_hg38)
            transcript = get_transcript(info['Feature'], transcripts_hg38)
        
        consequence = vep_consequence_trans(info['Consequence'])
        vcf_id = "-".join([vcfrecord.chrom, str(vcfrecord.pos), vcfrecord.ref, vcfrecord.alt])
        
        if consequence in lof_type and transcript:
            lof_pvs1 = PVS1(vcfrecord, consequence, info['HGVSc'], info['HGVSp'], transcript, genome_version)
            trans_name = lof_pvs1.transcript.full_name

            print(vcf_id,
                  info['SYMBOL'],
                  info['Feature'],
                  trans_name,
                  lof_pvs1.consequence,
                  lof_pvs1.strength_raw.name,
                  lof_pvs1.strength.name,
                  lof_pvs1.criterion,
                  sep="\t")
        elif consequence in lof_type:
            print(vcf_id,
                  info['SYMBOL'],
                  info['Feature'],
                  'not_canonical',
                  consequence,
                  'Unmet',
                  'Unmet',
                  'na',
                  sep="\t")
        else:
            print(vcf_id,
                  info['SYMBOL'],
                  info['Feature'],
                  'not_lof',
                  consequence,
                  'Unmet',
                  'Unmet',
                  'na',
                  sep="\t")

        '''
        if info['Function'] in ['splice-5', 'splice-3']:
            pvs1 = PVS1(vcfrecord, info['Function'], info['pHGVS1'], info['MIMInheritance'], transcript)
            splice = Splicing(vcfrecord, transcript)
            print(vcf_id,
                  pvs1.function,
                  pvs1.strength_raw.name,
                  pvs1.strength.name,
                  splice.is_undergo_NMD,
                  splice.has_cryptic_splice_site,
                  splice.is_exon_skipping,
                  splice.preserves_reading_frame,
                  splice.is_critical_to_protein_func,
                  splice.variant_removes_10_percent_of_protein,
                  splice.is_critical_to_protein_func_detail,
                  sep="\t")
        '''


if __name__ == '__main__':
    main()
