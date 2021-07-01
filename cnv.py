#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# created: 2019/6/27 21:54


from collections import namedtuple
from .strength import Strength
from .utils import contained_in_bed
from .read_data import pvs1_levels
from .read_data import domain_hg19, hotspot_hg19, curated_region_hg19, exon_lof_popmax_hg19
from .read_data import domain_hg38, hotspot_hg38, curated_region_hg38, exon_lof_popmax_hg38


CNVRecord = namedtuple('CNVRecord', ('chrom', 'start', 'end', 'type'))


class PVS1CNV:
    """
    CNV Dectection
    """
    def __init__(self, cnvrecord, consequence, transcript, genome_version):
        self.chrom = cnvrecord.chrom
        self.start = cnvrecord.start
        self.end = cnvrecord.end
        self.type = cnvrecord.type
        self.consequence = consequence
        self.transcript = transcript
        
        if genome_version in ['hg19', 'GRCh37']:
            self.domain = domain_hg19
            self.hotspot = hotspot_hg19
            self.curated_region = curated_region_hg19
            self.exon_lof_popmax = exon_lof_popmax_hg19
        else:
            self.domain = domain_hg38
            self.hotspot = hotspot_hg38
            self.curated_region = curated_region_hg38
            self.exon_lof_popmax = exon_lof_popmax_hg38

        self.overlap_exon_index = []
        self.overlap_exon = []
        self.overlap_exon_len = []
        self.overlap_cds_index = []
        self.overlap_cds = []
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
            if overlap_exon > 0:
                self.overlap_exon_index.append(i)
                self.overlap_exon.append((exon_start, exon_end))
                self.overlap_exon_len.append(overlap_exon)

        for j,  (cds_start, cds_end) in enumerate(self.transcript.cdslist):
            overlap_cds = min(self.end, cds_end) - max(self.start, cds_start)
            if overlap_cds > 0:
                self.overlap_cds_index.append(j)
                self.overlap_cds.append((cds_start, cds_end))
                self.overlap_cds_len.append(overlap_cds)

        if self.transcript.strand == '-':
            self.overlap_exon_index = [len(self.transcript.exon_sizes) - i - 1
                                       for i in self.overlap_exon_index]
            self.overlap_exon_index.reverse()
            self.overlap_exon.reverse()
            self.overlap_exon_len.reverse()
            self.overlap_cds_index = [len(self.transcript.cds_sizes) - i - 1
                                      for i in self.overlap_cds_index]
            self.overlap_cds_index.reverse()
            self.overlap_cds.reverse()
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

        chrom = self.chrom if 'chr' not in self.chrom else self.chrom.replace('chr', '')
        if self.transcript.strand == '+':
            start = self.start - 1
            end = self.transcript.cds_position.chrom_stop
        else:
            start = self.transcript.cds_position.chrom_start
            end = self.start
        in_domain = contained_in_bed(self.domain, chrom, start, end)
        in_hotspot = contained_in_bed(self.hotspot, chrom, start, end)
        in_curated_region = contained_in_bed(self.curated_region, chrom, start, end)

        is_func, desc = False, ''

        if in_curated_region:
            is_func = True
            desc = 'Expert curated region: {0}. '.format(in_curated_region[1])
        elif in_hotspot:
            is_func = True
            genomic_position, tag, missense_total, missense_PLP, missense_BLB = in_hotspot[1].split('|')
            desc = 'Mutational hotspot: <b>{0}</b> pathogenic missense variants and <b>{1}</b> benign ' \
                   'missense variant in <b>{2}</b>. '.format(missense_PLP, missense_BLB, genomic_position)
        if in_domain:
            (domain_name, amino_acids, genomic_position, tag,
             missense_total, missense_PLP, missense_BLB, block_position) = in_domain[1].split('|')
            # is_func = True if tag == 'WELL' else False
            if (int(missense_BLB) == 0 and int(missense_PLP) >= 5) or \
                    (int(missense_BLB) > 0 and int(missense_PLP) / int(missense_BLB) >= 10):
                is_func = True
            if missense_total == 0:
                desc += 'No ClinVar missense variant is found in <b>{0}</b> domain ({1}).'.format(domain_name, amino_acids)
            else:
                uniport_id = amino_acids.split(' ')[-1]
                amino_acid = ' '.join(amino_acids.split(' ')[:-1])
                desc += '<b>{0}</b> ClinVar pathogenic missense variant(s) and <b>{1}</b> benign missense variant(s) ' \
                        'are found in <b>{2}</b> domain ({3} <a href="https://www.uniprot.org/uniprot/{4}">{4}</a>).' \
                    .format(missense_PLP, missense_BLB, domain_name, amino_acid, uniport_id)
        if not in_hotspot and not in_domain:
            desc += 'Neither mutational hotspot nor functional domain is found.'

        return is_func, desc

    @property
    def exon_LoFs_are_frequent_in_pop(self):
        """
        LoF variants in this exon are frequent in the general population.
        discard inheritance consideration

        single LOF variant frequency >= 0.001
        multiple LOF variants frequency >= 0.005
        """
        return self.exon_lof_popmax_info[0]

    @property
    def exon_lof_popmax_desc(self):
        return self.exon_lof_popmax_info[1]

    @property
    def exon_lof_popmax_info(self):
        """
        NM_153818.1.4|1-2338205-G-A:1.19e-04|1-2338230-C-CT:5.62e-05
        :return:
        """
        in_exon_lof = contained_in_bed(self.exon_lof_popmax, self.chrom, self.start, self.end)
        if in_exon_lof:
            lof_list = in_exon_lof[1].split('|')
            lof_num = len(lof_list) - 1
            sum_freq = sum([float(i.split(':')[1]) for i in lof_list[1:]])
            max_lof, max_freq = lof_list[1].split(':')
            transcript, version, exon = lof_list[0].split('.')
            # if float(max_freq) > 0.001 or sum_freq > 0.005:
            if float(max_freq) > 0.001:
                desc = 'Maximum LOF population frequency in exon {0} of {1} is ' \
                       '<a href="https://gnomad.broadinstitute.org/variant/{2}">{3}</a>, ' \
                       'higher than the threshold (0.1%) we pre-defined.'.format(
                        exon, transcript + '.' + version, max_lof, max_freq)
                return True, desc
            else:
                desc = 'Maximum LOF population frequency in exon {0} of {1} is ' \
                       '<a href="https://gnomad.broadinstitute.org/variant/{2}">{3}</a>, ' \
                       'lower than the threshold (0.1%) we pre-defined.'.format(
                        exon, transcript + '.' + version, max_lof, max_freq)
                return False, desc
        else:
            return False, 'No LOF variant is found or the LOF variant dosen\'t exist in gnomAD.'

    @property
    def LoF_removes_more_than_10_percent_of_protein(self):
        """
        CNV removes >10% of protein.
        """
        cnv_cds_length = sum(self.overlap_cds_len)
        if self.transcript.cds_length > 0 and cnv_cds_length / self.transcript.cds_length > 0.1:
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
