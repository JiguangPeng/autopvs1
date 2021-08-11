#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# datetime: 2019/6/27 21:54

import re
from .strength import Strength
from .splicing import Splicing
from .utils import contained_in_bed
from .read_data import pvs1_levels
from .read_data import genome_hg19, domain_hg19, hotspot_hg19, curated_region_hg19, exon_lof_popmax_hg19, pathogenic_hg19
from .read_data import genome_hg38, domain_hg38, hotspot_hg38, curated_region_hg38, exon_lof_popmax_hg38, pathogenic_hg38


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

    def __init__(self, vcfrecord, consequence, cHGVS, pHGVS, transcript, genome_version):
        self.consequence = consequence
        self.vcfrecord = vcfrecord
        self.chrom = vcfrecord.chrom
        self.pos = vcfrecord.pos
        self.ref = vcfrecord.ref
        self.alt = vcfrecord.alt
        self.cHGVS = cHGVS
        self.pHGVS = pHGVS
        self.transcript = transcript

        if genome_version in ['hg19', 'GRCh37']:
            self.genome_version = 'hg19'
            self.vep_assembly = 'GRCh37'
            self.genome = genome_hg19
            self.domain = domain_hg19
            self.hotspot = hotspot_hg19
            self.curated_region = curated_region_hg19
            self.exon_lof_popmax = exon_lof_popmax_hg19
            self.pathogenic_dict = pathogenic_hg19
        else:
            self.genome_version = 'hg38'
            self.vep_assembly = 'GRCh38'
            self.genome = genome_hg38
            self.domain = domain_hg38
            self.hotspot = hotspot_hg38
            self.curated_region = curated_region_hg38
            self.exon_lof_popmax = exon_lof_popmax_hg38
            self.pathogenic_dict = pathogenic_hg38

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
                if self.get_pHGVS_termination < 374:
                    self.criterion = 'PTEN'
                    return Strength.VeryStrong
            if self.is_nmd_target:
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
            # CDH1 Guideline: Use the strong strength of evidence for canonical splice sites.
            if self.transcript.gene.name == 'CDH1':
                self.criterion = 'CDH1'
                return Strength.Strong
            splice = Splicing(self.vcfrecord, self.transcript, self.genome_version)
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
        :return: string
        """
        return self.functional_region[0]

    @property
    def func_desc(self):
        """
        Description for functional region
        :return:
        """
        if self.consequence in ['splice-5', 'splice-3']:
            splice = Splicing(self.vcfrecord, self.transcript, self.genome_version)
            return splice.func_desc
        else:
            return self.functional_region[1]

    @property
    def functional_region(self):
        """
        PM1 mutational hot spot and/or critical and well-established functional domain.
        """
        if self.transcript.gene.name == 'CDH1':
            is_func = self.get_pHGVS_termination <= 836
            if is_func:
                desc = '<a href="https://www.ncbi.nlm.nih.gov/pubmed/30311375">CDH1 gene-specific criteria</a>: ' \
                       'Truncations in NMD-resistant zone located upstream the most 3′ well-characterized ' \
                       'pathogenic variant c.2506G>T (p.Glu836Ter). '
            else:
                desc = '<a href="https://www.ncbi.nlm.nih.gov/pubmed/30311375">CDH1 gene-specific criteria</a>: ' \
                       'Truncations in NMD-resistant zone located downstream the most 3′ well-characterized ' \
                       'pathogenic variant c.2506G>T (p.Glu836Ter).'
            return is_func, desc

        chrom = self.chrom if 'chr' not in self.chrom else self.chrom.replace('chr', '')
        if self.transcript.strand == '+':
            start = int(self.pos) - 1
            end = int(self.transcript.cds_position.chrom_stop)
        else:
            start = int(self.transcript.cds_position.chrom_start)
            end = int(self.pos)

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

        # variant_score = 0
        # for pos in range(start, end + 1):
        #     if chrom + ":" + str(pos) in pathogenic_dict['score']:
        #         variant_score += pathogenic_dict['score'][chrom + ":" + str(pos)]
        # return in_domain or in_hotspot or in_curated_region or variant_score >= 3
        # return in_domain or in_hotspot or in_curated_region

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
            return Strength.Unset

        gene_name = self.transcript.gene.name
        if gene_name == 'MYH7':
            if self.strength_raw.value >= 2:
                return Strength.Moderate
            else:
                return self.strength_raw
        elif gene_name in pvs1_levels:
            if pvs1_levels[gene_name] == 'L0':
                return self.strength_raw
            elif pvs1_levels[gene_name] == 'L1':
                return self.strength_raw.downgrade(1)
            elif pvs1_levels[gene_name] == 'L2':
                return self.strength_raw.downgrade(2)
            elif pvs1_levels[gene_name] == 'L3':
                return Strength.Unmet
            else:
                return Strength.Unset
        else:
            return Strength.Unset

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
        c_pos, g_pos, interval = self.transcript.closest_potential_start_codon(self.genome)
        chrom = self.chrom if 'chr' not in self.chrom else self.chrom.replace('chr', '')
        if not c_pos:
            return Strength.VeryStrong, 'na', 0, 'IC0'
        variant_score = 0
        for start, end in interval:
            for pos in range(start, end):
                if chrom + ":" + str(pos) in self.pathogenic_dict['count']:
                    variant_score += self.pathogenic_dict['count'][chrom + ":" + str(pos)]
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
        new_stop_codon = self.get_pHGVS_termination
        cds_sizes = [i for i in self.transcript.cds_sizes if i > 0]
        if self.transcript.gene.name == 'GJB2':  # Hearing Loss Guidelines
            return True
        elif len(cds_sizes) <= 1:
            return False
        else:
            nmd_cutoff = sum(cds_sizes[:-1]) - min(50, cds_sizes[-2])
            return new_stop_codon * 3 <= nmd_cutoff

    @property
    def get_pHGVS_termination(self):
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
                # if int(match1.group(1)) / (self.transcript.cds_length/3) > 0.5:
                termination = int(match1.group(1)) + int(match1.group(3))
                # else:
                #    termination = int((self.transcript.cds_length/3)/2)
            elif match2:
                termination = int(match2.group(1))
            else:
                termination = -1

        elif '*' in self.pHGVS or 'X' in self.pHGVS or 'Ter' in self.pHGVS:
            pattern = re.compile(r'p\.\D+(\d+)(\*|X|Ter)')
            match = pattern.search(self.pHGVS)
            termination = int(match.group(1)) if match else -1
        else:
            termination = -1
        return termination

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
        start = int(self.pos) - 5
        end = int(self.pos) + 5
        in_exon_lof = contained_in_bed(self.exon_lof_popmax, self.chrom, start, end)
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
