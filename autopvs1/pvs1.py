#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# datetime: 2019/6/27 21:54

import re
from .strength import Strength
from .splicing import Splicing
from .utils import contained_in_bed
from .read_data import domain_bed, hotspot_bed, curated_region, exon_lof_frequent, \
    pathogenic_dict, pathogenic_dict2, pvs1_levels, genome


class PVS1:
    """
    very strong strength level for pathogenicity in the ACMG/AMP guideline:
    PVS1 - Null variant in a gene with established LOF as a disease mechanism;
    1) NMD classify: not occurring in the 3' most exon or the 3'-most 50 bp of the penultimate exon
    2) Cryptic Splice Site: nearby (+/- 20 nts) strong consensus splice sequence that may reconstitute in-frame splicing
    3) Truncated/altered region is critical to protein function:
        a) functional domain OR
        b) non-truncating pathogenic variants
    4) initial condon: have other functional trancripts?
       Pathogenic variant(s) upstream of closest potential in-frame start condon
    """

    def __init__(self, vcfrecord, consequence, cHGVS, pHGVS, transcript):
        self.consequence = consequence
        self.vcfrecord = vcfrecord
        self.chrom = vcfrecord.chrom
        self.pos = vcfrecord.pos
        self.ref = vcfrecord.ref
        self.alt = vcfrecord.alt
        self.cHGVS = cHGVS
        self.pHGVS = pHGVS
        self.transcript = transcript
        self.altcodon = 'na'
        self.init_path = 0
        self.criterion = 'na'
        self.strength_raw = self.verify_PVS1
        self.strength = self.adjust_PVS1

    @property
    def verify_PVS1(self):
        """
        verify_PVS1
        :return: Strength
        """
        if self.transcript is None or self.transcript.cds_length == 0:
            self.criterion = 'NF0'
            return Strength.Unmet

        if self.consequence in ['nonsense', 'frameshift']:
            if self.transcript.gene.name == 'PTEN':
                if self.get_pHGVS_terminatiton < 374:
                    self.criterion = 'PTEN'
                    return Strength.VeryStrong
            if self.is_nmd_target or self.transcript.exon_count == 1:
                if self.is_biologically_relevant:
                    self.criterion = 'NF1'
                    strength_raw = Strength.VeryStrong
                else:
                    strength_raw = Strength.Unmet
                    self.criterion = 'NF2'
            else:
                if self.is_critical_to_protein_func:
                    self.criterion = 'NF3'
                    strength_raw = Strength.Strong
                else:
                    if self.exon_LoFs_are_frequent_in_pop or \
                            not self.is_biologically_relevant:
                        self.criterion = 'NF4'
                        strength_raw = Strength.Unmet
                    elif self.LoF_removes_more_than_10_percent_of_protein:
                        self.criterion = 'NF5'
                        strength_raw = Strength.Strong
                    else:
                        self.criterion = 'NF6'
                        strength_raw = Strength.Moderate

        elif self.consequence in ['splice-5', 'splice-3']:
            if self.transcript.gene.name == 'PTEN':
                match = re.search(r'c\.(\d+)', self.cHGVS)
                if match and int(match.group(1))/3 < 374:
                    self.criterion = 'PTEN'
                    return Strength.VeryStrong
            splice = Splicing(self.vcfrecord, self.transcript)
            if splice.preserves_reading_frame:
                if splice.is_critical_to_protein_func:
                    self.criterion = 'SS10'
                    strength_raw = Strength.Strong
                elif self.exon_LoFs_are_frequent_in_pop or \
                        not self.is_biologically_relevant:
                    self.criterion = 'SS7'
                    strength_raw = Strength.Unmet
                elif splice.variant_removes_10_percent_of_protein:
                    self.criterion = 'SS8'
                    strength_raw = Strength.Strong
                else:
                    self.criterion = 'SS9'
                    strength_raw = Strength.Moderate
            elif splice.is_undergo_NMD:
                if self.is_biologically_relevant:
                    self.criterion = 'SS1'
                    strength_raw = Strength.VeryStrong
                else:
                    self.criterion = 'SS2'
                    strength_raw = Strength.Unmet
            else:
                if splice.is_critical_to_protein_func:
                    self.criterion = 'SS3'
                    strength_raw = Strength.Strong
                elif self.exon_LoFs_are_frequent_in_pop or \
                        not self.is_biologically_relevant:
                    self.criterion = 'SS4'
                    strength_raw = Strength.Unmet
                elif splice.variant_removes_10_percent_of_protein:
                    self.criterion = 'SS5'
                    strength_raw = Strength.Strong
                else:
                    self.criterion = 'SS6'
                    strength_raw = Strength.Moderate

        elif self.consequence == 'init-loss':
            strength_raw, self.altcodon, self.init_path, self.criterion = self.PVS1_start_codon()

        else:
            self.criterion = 'IC5'
            strength_raw = Strength.Unmet

        return strength_raw

    @property
    def is_critical_to_protein_func(self):
        """
        Truncated/altered region is critical to protein function.
        """
        if self.transcript.gene.name == 'CDH1':
            return self.get_pHGVS_terminatiton <= 836

        chrom = self.chrom if 'chr' not in self.chrom else self.chrom.replace('chr', '')
        if self.transcript.strand == '+':
            start = int(self.pos) - 1
            end = int(self.transcript.cds_position.chrom_stop)
        else:
            start = int(self.transcript.cds_position.chrom_start)
            end = int(self.pos)

        in_domain = contained_in_bed(domain_bed, chrom, start, end)
        in_hotspot = contained_in_bed(hotspot_bed, chrom, start, end)
        in_curated_region = contained_in_bed(curated_region, chrom, start, end)

        variant_score = 0
        for pos in range(start, end + 1):
            if chrom + ":" + str(pos) in pathogenic_dict:
                variant_score += pathogenic_dict[chrom + ":" + str(pos)]

        # return in_domain or in_hotspot or in_curated_region or variant_score >= 3
        return in_domain or in_hotspot or in_curated_region

    @property
    def adjust_PVS1(self):
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

    def PVS1_start_codon(self):
        """
        variant score:
        1 score   'practice_guideline': 4,
        1 score   'reviewed_by_expert_panel': 3,
        1 score   'criteria_provided,_multiple_submitters,_no_conflicts': 2,
        1/2 score 'criteria_provided,_conflicting_interpretations': 1,
        1/2 score 'criteria_provided,_single_submitter': 1,
        1/3 score 'no_interpretation_for_the_single_variant': 0,
        1/3 score 'no_assertion_criteria_provided': 0,
        1/3 score 'no_assertion_provided': 0

        strength level:
        1 variant with star 2/3/4
        2 variants with star 1
        3 variants with star 0
        :return: Strength, varaiant score
        """
        c_pos, g_pos, interval = self.transcript.closest_potential_start_codon(genome)
        chrom = self.chrom if 'chr' not in self.chrom else self.chrom.replace('chr', '')
        if not c_pos:
            return Strength.VeryStrong, 'na', 0, 'IC0'
        variant_score = 0
        for start, end in interval:
            for pos in range(start, end):
                if chrom + ":" + str(pos) in pathogenic_dict2:
                    variant_score += pathogenic_dict2[chrom + ":" + str(pos)]
        if variant_score < 1:
            return Strength.Supporting, c_pos, variant_score, 'IC4'
        elif 1 <= variant_score <= 3:
            return Strength.Moderate, c_pos, variant_score, 'IC3'
        elif 3 < variant_score <= 6:
            return Strength.Strong, c_pos, variant_score, 'IC2'
        else:
            return Strength.VeryStrong, c_pos, variant_score, 'IC1'

    @property
    def functional_transcript_use_alter_start_codon(self):
        """
        Different functional transcript uses alternative start codon
        :return: Boolean Value
        """
        # TODO: different functional transcript uses alternative start codon
        if self.transcript.name:
            return False
        else:
            return True

    @property
    def is_nmd_target(self):
        """
        Fuction: Nonsense-mediated decay (NMD) classification
        See more information about NMD: https://en.wikipedia.org/wiki/Nonsense-mediated_decay
        NMD classify: not occurring in the 3′ most exon or the 3′-most 50 bp of the penultimate exon
        :return: is NMD target or not
        """
        new_stop_codon = self.get_pHGVS_terminatiton
        if len(self.transcript.cds_sizes) == 1:
            return True
        elif len(self.transcript.cds_sizes) >= 2:
            if len([i for i in self.transcript.cds_sizes if i > 0]) == 1:
                return True
            else:
                nmd_cutoff = sum(self.transcript.cds_sizes[:-1]) - min(50, self.transcript.cds_sizes[-2])
        else:
            return

        if new_stop_codon * 3 > nmd_cutoff:
            return False
        else:
            return True

    @property
    def get_pHGVS_terminatiton(self):
        """
        Get termination position from pHGVS:
        # NM_031475.2:p.Gln98*
        # NM_031475.2:p.Ala586Glyfs*73
        # NP_000305.3:p.Arg378SerfsTer5

        # p.Arg97Glyfs*26 (alternatively p.Arg97GlyfsTer26, or short p.Arg97fs)
        """
        if 'fs' in self.pHGVS:
            pattern1 = re.compile(r'p\.\D+(\d+)\D+fs(\*|X|Ter)(\d+)')
            match1 = pattern1.search(self.pHGVS)
            pattern2 = re.compile(r'p\.\D+(\d+)fs')
            match2 = pattern2.search(self.pHGVS)

            if match1:
                if int(match1.group(1)) / (self.transcript.cds_length/3) > 0.5:
                    terminatiton = int(match1.group(1)) + int(match1.group(3))
                else:
                    terminatiton = int((self.transcript.cds_length/3)/2)
            elif match2:
                terminatiton = int(match2.group(1))
            else:
                terminatiton = -1

        elif '*' in self.pHGVS or 'X' in self.pHGVS or 'Ter' in self.pHGVS:
            pattern = re.compile(r'p\.\D+(\d+)(\*|X|Ter)')
            match = pattern.search(self.pHGVS)
            terminatiton = int(match.group(1)) if match else -1
        else:
            terminatiton = -1
        return terminatiton

    @property
    def is_biologically_relevant(self):
        """
        Trancript is biologically-relevant or not.
        """
        # TODO: construct a biologically-relevant transcript list
        # currently use the most biologically-relevant transcript
        # limitation: some gene may have multiple biologically-relevant transcript
        if self.transcript.name:
            return True
        else:
            return False

    @property
    def exon_LoFs_are_frequent_in_pop(self):
        """
        LoF variants in this exon are frequent in the general population.
        discard inheritance consideration

        PM2, BA1, BS1, BS2 = verify_PM2_BA1_BS1_BS2(self.vcfrecord, self.inheritance)
        return (BA1 == Strength.VeryStrong or
                BS1 == Strength.Strong or
                BS2 == Strength.Strong)
        """
        start = int(self.pos) - 5
        end = int(self.pos) + 5
        return contained_in_bed(exon_lof_frequent, self.chrom, start, end)

    @property
    def LoF_removes_more_than_10_percent_of_protein(self):
        """
        LoF variant removes >10% of protein.
        """
        pattern = re.compile(r'p\.\D+(\d+)(\D+fs)?(\*|X|Ter)(\d+)?')
        match = pattern.search(self.pHGVS)
        codon_offset = int(match.group(1)) if match else -1
        codon_length = self.transcript.cds_length / 3
        if codon_offset > 0 and (codon_length - codon_offset) / codon_length > 0.1:
            return True
        else:
            return False
