#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# created: 2019/2/2 16:28

import sys
from collections import namedtuple


VCFRecord = namedtuple('VCFRecord', ("chrom", "pos", 'ref', 'alt'))


def get_inheritance(genename):
    """
    get OMIM inheritance
    :param genename:
    :return:
    """
    return genename


lof_type = ['frameshift', 'stop_gained', 'splice_acceptor', 'splice_donor', 'start_lost']


def vep_consequence_trans(vep_consequence):
    """
    trans vep consequence name to custom name
    :param vep_consequence:
    :return:
    """
    if 'frameshift' in vep_consequence:
        return 'frameshift'
    elif 'stop_gained' in vep_consequence:
        return 'nonsense'
    elif 'splice_donor' in vep_consequence:
        return 'splice-5'
    elif 'splice_acceptor' in vep_consequence:
        return 'splice-3'
    elif 'start_lost' in vep_consequence:
        return 'init-loss'
    else:
        return vep_consequence


def vep2vcf(vep_name, genome):
    """
    convert vep Uploaded_variation name to VCFRecord format
    :param vep_name:
    :param genome:
    :return: VCFRecord
    """
    chrom, pos, refalt = vep_name.split("_")
    pos = int(pos)
    ref, alt = refalt.split("/")
    if ref == "-" or alt == "-":
        pos = pos - 1
        complement = genome[chrom][pos-1:pos].seq
        ref = (complement + ref).replace('-', '')
        alt = (complement + alt).replace('-', '')
    return VCFRecord(chrom.replace('chr', ''), pos, ref, alt)


def get_transcript(trans_name, transcripts):
    """
    Get a transcript using its name or a gene name
    :param trans_name:
    :param transcripts:
    :return: transcript
    """
    transcript = transcripts.get(trans_name)
    if not transcript:
        transcript = transcripts.get(trans_name.split('.')[0])
    if not transcript:
        transcript = None
    return transcript


def create_two_dim_dict(thedict, key1, key2, val):
    """
    add two dimension dict
    :param thedict:
    :param key1:
    :param key2:
    :param val:
    :return: dict
    """
    if key1 in thedict:
        thedict[key1].update({key2: val})
    else:
        thedict.update({key1: {key2: val}})


def create_bed_dict(bed):
    """
    read bed3 / bed6 / bed12 format as a dict
    :param bed:
    :return: bed dict
    """
    bed_dict = dict()
    try:
        with open(bed)as bed:
            for line in bed:
                record = line.strip().split("\t")
                chrom = record[0]
                if len(record) >= 12:
                    block_count = record[9]
                    block_sizes = record[10].split(",")
                    block_starts = record[11].split(",")
                    for i in range(int(block_count)):
                        start = int(record[1]) + int(block_starts[i])
                        end = start + int(block_sizes[i])
                        key = record[3] + '|' + str(start) + '-' + str(end)
                        create_two_dim_dict(bed_dict, key, "chrom", chrom)
                        create_two_dim_dict(bed_dict, key, "start", int(start))
                        create_two_dim_dict(bed_dict, key, "end", int(end))
                else:
                    start = int(record[1])
                    end = int(record[2])
                    if len(record) > 3:
                        key = record[3]
                    else:
                        key = chrom + ":" + str(start) + "-" + str(end)
                    create_two_dim_dict(bed_dict, key, "chrom", chrom)
                    create_two_dim_dict(bed_dict, key, "start", int(start))
                    create_two_dim_dict(bed_dict, key, "end", int(end))
        return bed_dict
    except Exception as e:
        sys.stderr.write(str(e))


def contained_in_bed(bed_dict, chrom, start, end):
    """
    Determine whether a variant is contained in the bed region.
    :param bed_dict:
    :param chrom:
    :param start:
    :param end:
    :return: Boolean value
    """
    chrom = str(chrom)
    if "chr" not in chrom:
        chrom = "chr" + chrom
    # max(start1, start2) < min(end1, end2)
    for key in bed_dict:
        if bed_dict[key]["chrom"] == chrom and \
                max(bed_dict[key]["start"], start) < min(bed_dict[key]["end"], end):
            return True, key
    return False


def read_morbidmap(file):
    """
    :param file: OMIM Phenotype
    :return: dict
    """
    _morbidmap_dict = {}
    try:
        with open(file) as fh:
            for line in fh:
                record = line.strip().split("\t")
                key = record[1] + '.' + record[2]
                if key not in _morbidmap_dict:
                    _morbidmap_dict[key] = record[3]
    except Exception as err:
        sys.stderr.write(err)
    return _morbidmap_dict


def read_pathogenic_site(file):
    """
    :param file: Pathogenic sites file
    :return: dict
    1 score   'practice_guideline': 4,
    1 score   'reviewed_by_expert_panel': 3,
    1 score   'criteria_provided,_multiple_submitters,_no_conflicts': 2,
    1/2 score 'criteria_provided,_conflicting_interpretations': 1,
    1/2 score 'criteria_provided,_single_submitter': 1,
    1/3 score 'no_interpretation_for_the_single_variant': 0,
    1/3 score 'no_assertion_criteria_provided': 0,
    1/3 score 'no_assertion_provided': 0
    """
    _pathogenic_dict = {}
    _pathogenic_dict['count'] = {}
    _pathogenic_dict['score'] = {}
    try:
        with open(file) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                record = line.strip().split("\t")
                key = record[0] + ':' + record[1]
                if record[6] in ['4', '3', '2']:
                    score = 1
                elif record[6] == '1':
                    score = 1/2
                else:
                    score = 1/3

                if key not in _pathogenic_dict['score']:
                    _pathogenic_dict['score'][key] = score
                else:
                    _pathogenic_dict['score'][key] += score

                if key not in _pathogenic_dict['count']:
                    _pathogenic_dict['count'][key] = 1
                else:
                    _pathogenic_dict['count'][key] += 1
    except Exception as err:
        sys.stderr.write(err)

    return _pathogenic_dict


def read_pvs1_levels(file):
    """
    :param file: PVS1 Levels
    :return: dict
    """
    _pvs1_levels = {}
    try:
        with open(file) as fh:
            for line in fh:
                record = line.strip().split("\t")
                gene = record[0]
                level = record[1]
                if gene not in _pvs1_levels:
                    _pvs1_levels[gene] = level
    except Exception as err:
        sys.stderr.write(err)

    return _pvs1_levels


def read_gene_alias(file):
    """
    :param file: gene alias
    :return: dict
    """
    _gene_alias = {}
    try:
        with open(file) as fh:
            for line in fh:
                record = line.strip().split("\t")
                _gene_alias[record[1]] = record[0]
    except Exception as err:
        sys.stderr.write(err)

    return _gene_alias
