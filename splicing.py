#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# created: 2019/6/27 17:54

import itertools
import re

from .pyhgvs.models import Transcript
from .maxentpy import maxent
from .maxentpy.maxent import load_matrix5, load_matrix3
from .utils import contained_in_bed
from .read_data import genome_hg19, transcripts_hg19, domain_hg19, hotspot_hg19, curated_region_hg19
from .read_data import genome_hg38, transcripts_hg38, domain_hg38, hotspot_hg38, curated_region_hg38


matrix5 = load_matrix5()
matrix3 = load_matrix3()


class Splicing:
    """
    splice class
    """
    donor_threshold = 3
    acceptor_threshold = 3
    percent_threshold = 0.7

    def __init__(self, vcfrecord, transcript, genome_version):
        self.chrom = vcfrecord.chrom
        self.offset = int(vcfrecord.pos)
        self.ref = vcfrecord.ref
        self.alt = vcfrecord.alt
        self.transcript = transcript
        
        if genome_version in ['hg19', 'GRCh37']:
            self.genome_version = 'hg19'
            self.vep_assembly = 'GRCh37'
            self.genome = genome_hg19
            self.transcripts = transcripts_hg19
            self.domain = domain_hg19
            self.hotspot = hotspot_hg19
            self.curated_region = curated_region_hg19
        else:
            self.genome_version = 'hg38'
            self.vep_assembly = 'GRCh38'
            self.genome = genome_hg38
            self.transcripts = transcripts_hg38
            self.domain = domain_hg38
            self.hotspot = hotspot_hg38
            self.curated_region = curated_region_hg38

        self.type = 'NA'
        self.index = 'NA'
        self.refseq = ''
        self.altseq = ''
        self.__parse()
        self.maxentscore_ref = -1.00
        self.maxentscore_alt = -1.00
        self.maxent_foldchange = 1.00
        self.__calculate_maxentscore()

    @staticmethod
    def get_transcript(transcript):
        if isinstance(transcript, Transcript):
            transcript = transcript
        elif isinstance(transcript, str):
            transcript = self.transcripts.get(transcript)
        else:
            transcript = None
        return transcript

    @staticmethod
    def format_donor(raw_seq):
        donor_exon = 3
        format_seq = raw_seq[:donor_exon].lower() + \
                     raw_seq[donor_exon:donor_exon + 2].upper() + \
                     raw_seq[donor_exon + 2:].lower()
        return format_seq

    @staticmethod
    def format_acceptor(raw_seq):
        acceptor_intron = 20
        format_seq = raw_seq[:acceptor_intron - 2].lower() + \
                     raw_seq[acceptor_intron - 2:acceptor_intron].upper() + \
                     raw_seq[acceptor_intron:].lower()
        return format_seq

    @staticmethod
    def reverse_complement(seq):
        """Retrun a reverse complementary seq"""
        nt_complement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C',
                         'a': 't', 'c': 'g', 't': 'a', 'g': 'c'}
        reverse_seq = list(reversed(seq))
        rev_comp_seq = [nt_complement[k] for k in reverse_seq]
        return ''.join(rev_comp_seq)

    def __calculate_maxentscore(self):
        """
        --- Calculate the maxentscan socre ---
        When a mutation occurs, if the WT score is above the threshold and
        the score variation (between WT and Mutant) is under -10% for HSF (-30% for MaxEnt)
        we consider that the mutation breaks the splice site.
        In the other case, if the WT score is under the threshold and
        the score variation is above +10% for HSF (+30% for MaxEnt) we consider that
        the mutation creates a new splice site.
        """
        maxentscore_alt = maxentscore_ref = -1.00
        if self.type == 'donor':
            if len(self.refseq) == 9:
                maxentscore_ref = maxent.score5(self.refseq, matrix=matrix5)
            if len(self.altseq) == 9:
                maxentscore_alt = maxent.score5(self.altseq, matrix=matrix5)
        elif self.type == 'acceptor':
            if len(self.refseq) == 23:
                maxentscore_ref = maxent.score3(self.refseq, matrix=matrix3)
            if len(self.altseq) == 23:
                maxentscore_alt = maxent.score3(self.altseq, matrix=matrix3)

        maxent_foldchange = maxentscore_alt / maxentscore_ref

        self.maxentscore_ref = round(maxentscore_ref, 2)
        self.maxentscore_alt = round(maxentscore_alt, 2)
        self.maxent_foldchange = round(maxent_foldchange, 2)

    def __parse(self):
        """
        Get the refseq and altseq around splice sites
        """
        donor_exon = 3
        donor_intron = 6
        acceptor_exon = 3
        acceptor_intron = 20
        chrom = self.chrom if 'chr' in self.chrom else 'chr' + self.chrom
        offset = self.offset
        ref = self.ref.strip('.')
        alt = self.alt.strip('.')
        refseq = ''
        altseq = ''
        refseq_start = None
        refseq_end = None
        splice_type = 'NA'
        index = 'NA'
        transcript = self.transcript
        intron_count = len(transcript.intronlist)

        for i, (intron_start, intron_end) in enumerate(transcript.intronlist):
            for offset_i in range(offset, offset+len(self.ref)):
                to_start = offset_i - intron_start
                to_end = offset_i - intron_end
                altseq_exon_end = altseq_intron_end = ''
                if transcript.strand == '+':
                    if -donor_exon < to_start <= donor_intron:
                        splice_type = 'donor'
                        refseq_start = intron_start - donor_exon
                        refseq_end = intron_start + donor_intron
                        refseq = self.genome[chrom][refseq_start:refseq_end].seq
                        if intron_start - donor_exon < offset - 1:
                            altseq_exon_end = self.genome[chrom][refseq_start:offset - 1].seq
                        if offset + len(ref) - 1 < refseq_end + len(ref) - len(alt):
                            altseq_intron_end = self.genome[chrom][offset + len(ref) - 1:
                                                              refseq_end + len(ref) - len(alt)].seq
                        altseq = altseq_exon_end + alt + altseq_intron_end

                        if to_start > 0:
                            index = 'IVS' + str(i + 1) + '+' + str(to_start)
                        else:
                            index = 'EX' + str(i + 1) + '-' + str(1 - to_start)

                    if -acceptor_intron < to_end <= acceptor_exon:
                        splice_type = 'acceptor'
                        refseq_start = intron_end - acceptor_intron
                        refseq_end = intron_end + acceptor_exon
                        refseq = self.genome[chrom][refseq_start:refseq_end].seq
                        if offset + len(ref) - 1 < refseq_end:
                            altseq_exon_end = self.genome[chrom][offset + len(ref) - 1:refseq_end].seq
                        if refseq_start - len(ref) + len(alt) < offset - 1:
                            altseq_intron_end = self.genome[chrom][refseq_start - len(ref) + len(alt):offset - 1].seq
                        altseq = altseq_intron_end + alt + altseq_exon_end
                        if to_end > 0:
                            index = 'EX' + str(i + 2) + '+' + str(to_end)
                        else:
                            index = 'IVS' + str(i + 1) + '-' + str(1 - to_end)
                else:
                    if -acceptor_exon < to_start <= acceptor_intron:
                        splice_type = 'acceptor'
                        refseq_start = intron_start - acceptor_exon
                        refseq_end = intron_start + acceptor_intron
                        refseq = self.genome[chrom][refseq_start:refseq_end].reverse.complement.seq
                        if refseq_start < offset - 1:
                            altseq_exon_end = self.genome[chrom][refseq_start:offset - 1].seq
                        if offset + len(ref) - 1 < refseq_end + len(ref) - len(alt):
                            altseq_intron_end = self.genome[chrom][offset + len(ref) - 1:
                                                                   refseq_end + len(ref) - len(alt)].seq
                        altseq = altseq_exon_end + alt + altseq_intron_end
                        altseq = self.reverse_complement(altseq)

                        if to_start > 0:
                            index = 'IVS' + str(intron_count - i) + '-' + str(to_start)
                        else:
                            index = 'EX' + str(intron_count - i + 1) + '+' + str(1 - to_start)

                    if -donor_intron < to_end <= donor_exon:
                        splice_type = 'donor'
                        refseq_start = intron_end - donor_intron
                        refseq_end = intron_end + donor_exon
                        refseq = self.genome[chrom][refseq_start:refseq_end].reverse.complement.seq
                        if offset + len(ref) - 1 < refseq_end:
                            altseq_exon_end = self.genome[chrom][offset + len(ref) - 1:refseq_end].seq
                        if refseq_start - len(ref) + len(alt) < offset - 1:
                            altseq_intron_end = self.genome[chrom][refseq_start - len(ref) + len(alt):
                                                                   offset - 1].seq
                        altseq = altseq_intron_end + alt + altseq_exon_end
                        altseq = self.reverse_complement(altseq)

                        if to_end > 0:
                            index = 'EX' + str(intron_count - i) + '-' + str(to_end)
                        else:
                            index = 'IVS' + str(intron_count - i) + '+' + str(1 - to_end)

            # Format upper and lower case for better demonstration
            if splice_type == 'donor':
                refseq = self.format_donor(refseq)
                altseq = self.format_donor(altseq)
            elif splice_type == 'acceptor':
                refseq = self.format_acceptor(refseq)
                altseq = self.format_acceptor(altseq)
        self.type = splice_type
        self.index = index
        self.refseq = refseq
        self.altseq = altseq
        self.refseq_start = refseq_start
        self.refseq_end = refseq_end
        # print(self.transcript.full_name, len(self.refseq), len(self.altseq))

    @property
    def cryptic_splice_site(self):
        """
        Search for cryptic splice site
        1) nearby (+/- 20 nts) strong consensus splice sequence
        2) reconstitutes or disrupts in-frame splicing
        3) undergo NMD or not
        Consensus values go from 0 to 100 for HSF, -20 to +20 for MaxEnt.
        The threshold is defined at 65 for HSF, 3 for MaxEnt.
        This means that every signal with a score above the threshold is considered
        to be a splice site (donor or acceptor).
        Cite: http://www.umd.be/HSF3/technicaltips.html
        """
        refscore = self.maxentscore_ref
        chrom = self.chrom if 'chr' in self.chrom else 'chr' + self.chrom
        search_flank = 50
        list1 = list(range(self.refseq_start-1, self.refseq_start-1-search_flank, -1))
        list2 = list(range(self.refseq_start+1, self.refseq_start+1+search_flank, 1))
        search_region = list(itertools.chain.from_iterable(zip(list1, list2)))

        for pos in search_region:
            if self.type == 'donor':
                splice_context = self.genome[chrom][pos: pos + 9].seq
                alt_index = self.offset - pos - 1
                if 0 < alt_index < 9:
                    splice_context = splice_context[:alt_index] + self.alt + \
                                     splice_context[alt_index + len(self.alt):10-len(self.alt)]
                if self.transcript.strand == '-':
                    splice_context = self.reverse_complement(splice_context)

                splice_context = self.format_donor(splice_context)
                if len(splice_context) == 9:
                    maxentscore = maxent.score5(splice_context, matrix=matrix5)
                else:
                    maxentscore = 0
                if splice_context[3:5] in ['GT', self.refseq[3:5]] and maxentscore > 1 and \
                        (maxentscore >= self.donor_threshold or
                         maxentscore / refscore >= self.percent_threshold):
                    return pos, splice_context, maxentscore

            elif self.type == 'acceptor':
                splice_context = self.genome[chrom][pos: pos + 23].seq
                alt_index = self.offset - pos - 1
                if 0 < alt_index < 23:
                    splice_context = splice_context[:alt_index] + self.alt + \
                                     splice_context[alt_index + len(self.alt):24-len(self.alt)]
                if self.transcript.strand == '-':
                    splice_context = self.reverse_complement(splice_context)

                splice_context = self.format_acceptor(splice_context)
                if len(splice_context) == 23:
                    maxentscore = maxent.score3(splice_context, matrix=matrix3)
                else:
                    maxentscore = 0
                if splice_context[18:20] in ['AG', self.refseq[18:20]] and maxentscore > 1 and \
                        (maxentscore >= self.acceptor_threshold or
                         maxentscore / refscore >= self.percent_threshold):
                    return pos, splice_context, maxentscore
        return 0, '', 0

    @property
    def has_cryptic_splice_site(self):
        if self.type != 'NA' and self.cryptic_splice_site[0] > 0:
            cryptic_exon = self.cryptic_coding_exons[self.skipped_exon_id - 1]
            if cryptic_exon and cryptic_exon.chrom_start < cryptic_exon.chrom_end:
                return True
        return False

    @property
    def is_exon_skipping(self):
        if (not self.has_cryptic_splice_site and
                self.maxent_foldchange < self.percent_threshold and
                self.maxentscore_alt < 3):
            return True
        else:
            return False

    @property
    def skipped_exon_length(self):
        splice_match = re.match(r'IVS(\d+)([+|-])(\d+)', self.index)
        splice_match2 = re.match(r'EX(\d+)([+|-])(\d+)', self.index)
        if splice_match:
            intron_id = int(splice_match.group(1))
            if splice_match.group(2) == '+':
                return self.transcript.cds_sizes[intron_id - 1]
            else:
                return self.transcript.cds_sizes[intron_id]
        elif splice_match2:
            exon_id = int(splice_match2.group(1))
            return self.transcript.cds_sizes[exon_id - 1]
        else:
            return 0

    @property
    def skipped_exon_id(self):
        splice_match = re.match(r'IVS(\d+)([+|-])(\d+)', self.index)
        splice_match2 = re.match(r'EX(\d+)([+|-])(\d+)', self.index)
        if splice_match:
            intron_id = int(splice_match.group(1))
            if splice_match.group(2) == '+':
                return intron_id
            else:
                return intron_id + 1
        if splice_match2:
            exon_id = int(splice_match2.group(1))
            return exon_id
        else:
            return 0

    @property
    def preserves_reading_frame(self):
        if self.is_exon_skipping:
            return self.skipped_exon_length % 3 == 0
        elif self.has_cryptic_splice_site:
            return (self.cryptic_splice_site[0] - self.refseq_start) % 3 == 0
        else:
            return True

    @property
    def cryptic_coding_exons(self):
        cryptic_coding_exons = self.transcript.coding_exons
        if cryptic_coding_exons[self.skipped_exon_id - 1] is None:
            return cryptic_coding_exons

        if (self.transcript.strand == '+' and self.type == 'donor') or \
                (self.transcript.strand == '-' and self.type == 'acceptor'):
            cryptic_coding_exons[self.skipped_exon_id - 1] = \
                cryptic_coding_exons[self.skipped_exon_id - 1]._replace(
                    chrom_end=self.cryptic_splice_site[0] + 3)
        elif self.type == 'acceptor':
            cryptic_coding_exons[self.skipped_exon_id - 1] = \
                cryptic_coding_exons[self.skipped_exon_id - 1]._replace(
                    chrom_start=self.cryptic_splice_site[0] + 20)
        else:
            cryptic_coding_exons[self.skipped_exon_id - 1] = \
                cryptic_coding_exons[self.skipped_exon_id - 1]._replace(
                    chrom_start=self.cryptic_splice_site[0] + 6)
        return cryptic_coding_exons

    @property
    def get_trans_seq_info(self):
        trans_seq = ''
        cds_sizes = []
        for exon in self.cryptic_coding_exons:
            if exon and exon.strand == '+':
                cds_sizes.append(exon.chrom_end - exon.chrom_start)
                trans_seq += self.genome[exon.chrom][exon.chrom_start:exon.chrom_end].seq
            elif exon and exon.strand == '-':
                cds_sizes.append(exon.chrom_end - exon.chrom_start)
                trans_seq += self.genome[exon.chrom][exon.chrom_start:
                                                     exon.chrom_end].reverse.complement.seq
            else:
                cds_sizes.append(0)

        stop_codon = 0
        for pos in range(0, len(trans_seq), 3):
            if trans_seq[pos:pos + 3].upper() in ['TAA', 'TAG', 'TGA']:
                stop_codon = pos
                break

        is_nmd_target = False
        if len(cds_sizes) == 1 or len([i for i in cds_sizes if i > 0]) == 1:
            is_nmd_target = True
        else:
            nmd_cutoff = sum(cds_sizes[:-1]) - min(50, cds_sizes[-2])
            if stop_codon <= nmd_cutoff:
                is_nmd_target = True

        return trans_seq, stop_codon, is_nmd_target

    @property
    def is_undergo_NMD(self):
        if self.preserves_reading_frame:
            return False
        elif self.has_cryptic_splice_site:
            return self.get_trans_seq_info[2]
        elif self.skipped_exon_id >= self.transcript.exon_count - 1:
            return False
        else:
            return True

    @property
    def is_critical_to_protein_func(self):
        return self.is_critical_to_protein_func_detail[0]

    @property
    def func_desc(self):
        return self.is_critical_to_protein_func_detail[1]

    @property
    def is_critical_to_protein_func_detail(self):
        """
        Truncated/altered region is critical to protein function.
        """
        chrom = self.chrom if 'chr' not in self.chrom else self.chrom.replace('chr', '')
        if self.has_cryptic_splice_site and self.preserves_reading_frame:
            if self.transcript.strand == '+':
                if self.type == 'acceptor':
                    start = self.transcript.exonlist[self.skipped_exon_id-1][0]
                    end = self.cryptic_splice_site[0]
                else:
                    start = self.cryptic_splice_site[0]
                    end = self.transcript.exonlist[self.skipped_exon_id-1][1]
            else:
                if self.type == 'acceptor':
                    start = self.transcript.exonlist[self.transcript.exon_count - self.skipped_exon_id][0]
                    end = self.cryptic_splice_site[0]
                else:
                    start = self.cryptic_splice_site[0]
                    end = self.transcript.exonlist[self.transcript.exon_count - self.skipped_exon_id][1]
            if start >= end:
                return False, 'NA'
        elif self.has_cryptic_splice_site and not self.preserves_reading_frame:
            # TODO: new stop codon position
            if self.transcript.strand == '+':
                start = self.cryptic_splice_site[0]
                end = self.transcript.cds_position.chrom_stop
            else:
                start = self.transcript.cds_position.chrom_start
                end = self.cryptic_splice_site[0]
        elif self.is_exon_skipping:
            if self.transcript.strand == '+':
                start = self.transcript.exonlist[self.skipped_exon_id-1][0]
                end = self.transcript.exonlist[self.skipped_exon_id-1][1]
            else:
                start = self.transcript.exonlist[self.transcript.exon_count - self.skipped_exon_id][0]
                end = self.transcript.exonlist[self.transcript.exon_count - self.skipped_exon_id][1]
        else:
            return False, 'NA'

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
    def variant_removes_10_percent_of_protein(self):
        if self.has_cryptic_splice_site:
            if self.transcript.strand == '+':
                start = self.refseq_start
                end = self.cryptic_splice_site[0]
                if (start - end) / self.transcript.cds_length > 0.1:
                    return True
                else:
                    return False
            else:
                start = self.cryptic_splice_site[0]
                end = self.refseq_start
                if (start - end) / self.transcript.cds_length > 0.1:
                    return True
                else:
                    return False
        elif self.skipped_exon_length / self.transcript.cds_length > 0.1:
            return True
        else:
            return False
