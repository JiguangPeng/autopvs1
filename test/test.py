#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# datetime: 2019/6/27 19:13

import sys
from collections import namedtuple
sys.path.append('../../')
from autopvs1 import PVS1
from autopvs1.read_data import transcripts, genome, trans_gene, gene_trans
from autopvs1.utils import vep2vcf, get_transcript, vep_consequence_trans, VCFRecord

lof_type = ['frameshift', 'nonsense', 'splice-5', 'splice-3', 'init-loss']
vep_lof_list = ['frameshift_variant',
                'stop_gained',
                'splice_donor_variant',
                'splice_acceptor_variant',
                'start_lost']

VAR = namedtuple('VAR', ('varid', 'gene', 'trans', 'canonical', 'pick', 'record'))


def main():
    anno_fh = open(sys.argv[1])
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

        if 'VCFRecord' in info:
            chrom, pos, ref, alt = info['VCFRecord'].split('-')
            vcfrecord = VCFRecord(chrom, pos, ref, alt)
        elif 'Uploaded_variation' in info:
            vcfrecord = vep2vcf(info['Uploaded_variation'], genome)
        else:
            raise IOError
        vcf_id = "-".join([vcfrecord.chrom, str(vcfrecord.pos), vcfrecord.ref, vcfrecord.alt])

        transcript = get_transcript(info['Feature'], transcripts)
        consequence = vep_consequence_trans(info['Consequence'])

        if consequence in lof_type and transcript:
            lof_pvs1 = PVS1(vcfrecord, consequence, info['HGVSc'], info['HGVSp'], transcript)

            print(vcf_id,
                  info['SYMBOL'],
                  info['Feature'],
                  'canonical',
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
         if consequence == 'init-loss' and transcript:
            lof_pvs1 = PVS1(vcfrecord, consequence, info['HGVSc'], info['HGVSp'], transcript)

            print(vcf_id,
                  info['SYMBOL'],
                  info['Feature'],
                  lof_pvs1.transcript.full_name,
                  lof_pvs1.consequence,
                  lof_pvs1.strength_raw.name,
                  lof_pvs1.strength.name,
                  lof_pvs1.criterion,
                  lof_pvs1.transcript.cds_length,
                  lof_pvs1.altcodon,
                  lof_pvs1.init_path,
                  sep="\t")
        
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
