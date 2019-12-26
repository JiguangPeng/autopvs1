#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# datetime: 2019/6/27 21:54


from collections import namedtuple

from .read_data import domain_bed, hotspot_bed, curated_region, pathogenic_dict, exon_lof_frequent, pvs1_levels
from .strength import Strength
from .utils import contained_in_bed

CNVRecord = namedtuple('CNVRecord', ('chrom', 'start', 'end', 'type'))


class PVS1CNV:
    """
    CNV Dectection
    """
    def __init__(self, cnvrecord, consequence, transcript):
        self.chrom = cnvrecord.chrom
        self.start = cnvrecord.start
        self.end = cnvrecord.end
        self.type = cnvrecord.type
        self.consequence = consequence
        self.transcript = transcript
        self.overlap_exon_index = []
        self.overlap_exon = []
        self.overlap_exon_len = []
        self.overlap_cds_len = []
        if self.transcript:
            self.cnv_exon_overlap()
            if self.type in ['DEL', 'del', 'Del']:
                self.strength_raw, self.criterion = self.verify_DEL()
            else:
                self.strength_raw, self.criterion = self.verify_DUP()
        else:
            self.strength_raw = Strength.Unmet
            self.criterion = 'NA'
        self.strength = self.strength_adjust()

    def cnv_exon_overlap(self):
        for i, (exon_start, exon_end) in enumerate(self.transcript.exonlist):
            overlap_exon = min(self.end, exon_end) - max(self.start, exon_start)
            if len(self.transcript.cdslist) > i:
                cds_start, cds_end = self.transcript.cdslist[i]
                overlap_cds = min(self.end, cds_end) - max(self.start, cds_start)
                overlap_cds = overlap_cds if overlap_cds > 0 else 0
            else:
                overlap_cds = 0
            if overlap_exon > 0:
                self.overlap_exon_index.append(i)
                self.overlap_exon.append((exon_start, exon_end))
                self.overlap_exon_len.append(overlap_exon)
                self.overlap_cds_len.append(overlap_cds)
        if self.transcript.strand == '-':
            self.overlap_exon_index = [len(self.transcript.exon_sizes) - i - 1
                                       for i in self.overlap_exon_index]
            self.overlap_exon_index.reverse()
            self.overlap_exon.reverse()
            self.overlap_exon_len.reverse()
            self.overlap_cds_len.reverse()

    def strength_adjust(self):
        """
        Criteria for LoF disease mechanism
        1) Follow PVS1 Flowchart if 'L0'
        2) Decrease final strength by one level (i.e. VeryStrong to Strong) if 'L1'
        3) Decrease final strength by two levels (i.e. VeryStrong to Moderate) if 'L2'
        4) If there is no evidence that LOF variants cause disease, PVS1 should not
           be applied at any strength level.
        :return: Strength Level
        """
        if self.transcript is None:
            return Strength.Unmet

        gene_name = self.transcript.gene.name
        if gene_name in pvs1_levels:
            if pvs1_levels[gene_name] == 'L0':
                return self.strength_raw
            elif pvs1_levels[gene_name] == 'L1':
                return self.strength_raw.downgrade(1)
            elif pvs1_levels[gene_name] == 'L2':
                return self.strength_raw.downgrade(2)
            else:
                return Strength.Unmet
        else:
            return Strength.Unmet

    @property
    def del_preserve_reading_frame(self):
        return sum(self.overlap_cds_len) % 3 == 0

    @property
    def new_stop_codon(self):
        if not self.del_preserve_reading_frame:
            index = self.overlap_exon_index[-1]
            return sum(self.transcript.cds_sizes[0:index])
        else:
            return self.transcript.cds_length

    @property
    def del_undergo_NMD(self):
        if len(self.transcript.cds_sizes) == 1:
            return True
        elif len(self.transcript.cds_sizes) >= 2:
            if len([i for i in self.transcript.cds_sizes if i > 0]) == 1:
                return True
            else:
                nmd_cutoff = sum(self.transcript.cds_sizes[:-1]) - min(50, self.transcript.cds_sizes[-2])
        else:
            return

        if self.new_stop_codon * 3 > nmd_cutoff:
            return False
        else:
            return True

    @property
    def is_critical_to_protein_func(self):
        return self.functional_region[0]

    @property
    def func_desc(self):
        """
        Description for functional region
        :return:
        """
        return self.functional_region[1]

    @property
    def functional_region(self):
        """
        Truncated/altered region is critical to protein function.
        :return:
        """
        if self.transcript.gene.name == 'CDH1':
            is_func = self.get_pHGVS_termination <= 836
            desc = 'Truncations in NMD-resistant zone located upstream the most 3′ well-characterized ' \
                   'pathogenic variant c.2506G>T (p.Glu836Ter).'
            return is_func, desc

        chrom = self.chrom if 'chr' not in self.chrom else self.chrom.replace('chr', '')
        if self.transcript.strand == '+':
            start = self.start - 1
            end = self.transcript.cds_position.chrom_stop
        else:
            start = self.transcript.cds_position.chrom_start
            end = self.start
        in_domain = contained_in_bed(domain_bed, chrom, start, end)
        in_hotspot = contained_in_bed(hotspot_bed, chrom, start, end)
        in_curated_region = contained_in_bed(curated_region, chrom, start, end)

        is_func, desc = False, ''

        if in_curated_region:
            is_func = True
            desc = 'Expert curated region: {0}. '.format(in_curated_region[1])
        elif in_hotspot:
            is_func = True
            genomic_position, tag, missense_total, missense_PLP, missense_BLB = in_hotspot[1].split('|')
            desc = 'mutational hotspot: {0} pathogenic missense variant and ' \
                   '{1} benign missense variant in {2}. '.format(missense_PLP, missense_BLB, genomic_position)
        if in_domain:
            (domain_name, amino_acids, genomic_position, tag,
             missense_total, missense_PLP, missense_BLB) = in_domain[1].split('|')
            is_func = True if tag == 'WELL' else False
            if missense_total == 0:
                desc += 'No missense variant found in domain: {0} ({1}).'.format(domain_name, amino_acids)
            else:
                desc += '{2} pathogenic missense variant and {3} benign missense variant ' \
                        'found in domain: {0} ({1}).'.format(domain_name, amino_acids, missense_PLP, missense_BLB)
        if not in_hotspot and not in_domain:
            desc += 'No mutational hotspot or functional domain found.'

        return is_func, desc


    @property
    def exon_LoFs_are_frequent_in_pop(self):
        """
        LoF variants in this exon are not frequent in the general population.
        """
        return contained_in_bed(exon_lof_frequent, self.chrom, self.start, self.end)

    @property
    def LoF_removes_more_than_10_percent_of_protein(self):
        """
        CNV removes >10% of protein.
        """
        cnv_cds_length = sum(self.overlap_cds_len)
        if cnv_cds_length / self.transcript.cds_length > 0.1:
            return True
        else:
            return False

    @property
    def is_biologically_relevant(self):
        """
        Trancript is biologically-relevant or not.
        """
        if self.transcript.name:
            return True
        else:
            return False

    def verify_DEL(self):
        if len(self.overlap_exon_len) == len(self.transcript.exon_sizes):
            return Strength.VeryStrong, 'DEL1'
        elif not self.del_preserve_reading_frame:
            if self.del_undergo_NMD:
                if self.is_biologically_relevant:
                    return Strength.VeryStrong, 'DEL2'
                else:
                    return Strength.Unmet, 'DEL3'
            elif self.is_critical_to_protein_func:
                return Strength.Strong, 'DEL4'
            elif self.exon_LoFs_are_frequent_in_pop or not self.is_biologically_relevant:
                return Strength.Unmet, 'DEL5'
            elif self.LoF_removes_more_than_10_percent_of_protein:
                return Strength.Strong, 'DEL6'
            else:
                return Strength.Moderate, 'DEL7'
        elif self.is_critical_to_protein_func:
            return Strength.Strong, 'DEL8'
        elif self.exon_LoFs_are_frequent_in_pop or not self.is_biologically_relevant:
            return Strength.Unmet, 'DEL9'
        elif self.LoF_removes_more_than_10_percent_of_protein:
            return Strength.Strong, 'DEL10'
        else:
            return Strength.Moderate, 'DEL11'

    def verify_DUP(self):
        """
        Duplication: >=1 exon in size and must be completely contained with gene
        TDUP: proven in tadem
        DUP: presumed in tandem
        NTDUP: preven not in tandem
        :return: Strength
        """
        if self.type == 'TDUP':
            return Strength.VeryStrong, 'DUP1'
        elif self.type == 'DUP':
            return Strength.Strong, 'DUP3'
        else:
            return Strength.Unmet, 'DUP5'
