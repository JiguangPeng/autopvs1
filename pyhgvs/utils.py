"""
Helper functions.
"""

from .models import Exon
from .models import Position
from .models import Transcript


def read_refgene(infile):
    """
    Iterate through a refGene file.

    GenePred extension format:
    http://genome.ucsc.edu/FAQ/FAQformat.html#GenePredExt

    Column definitions:
    0. uint undocumented id
    1. string name;                 "Name of gene (usually transcript_id from GTF)"
    2. string chrom;                "Chromosome name"
    3. char[1] strand;              "+ or - for strand"
    4. uint txStart;                "Transcription start position"
    5. uint txEnd;                  "Transcription end position"
    6. uint cdsStart;               "Coding region start"
    7. uint cdsEnd;                 "Coding region end"
    8. uint exonCount;              "Number of exons"
    9. uint[exonCount] exonStarts;  "Exon start positions"
    10. uint[exonCount] exonEnds;   "Exon end positions"
    11. int score;                  "Score"
    12. string name2;               "Alternate name (e.g. gene_id from GTF)"
    13. string cdsStartStat;        "enum('none','unk','incmpl','cmpl')"
    14. string cdsEndStat;          "enum('none','unk','incmpl','cmpl')"
    15. lstring exonFrames;         "Exon frame offsets {0,1,2}"
    """
    for line in infile:
        if line.startswith('#'):
            continue
        row = line.rstrip('\n').split('\t')
        exon_starts = map(int, row[9].strip(',').split(','))
        exon_ends = map(int, row[10].strip(',').split(','))
        exons = list(zip(exon_starts, exon_ends))

        yield {
            'bin': row[0],
            'name': row[1],
            'chrom': row[2],
            'strand': row[3],
            'txStart': int(row[4]),
            'txEnd': int(row[5]),
            'cdsStart': int(row[6]),
            'cdsEnd': int(row[7]),
            'exonCount': int(row[8]),
            'exonStarts': [int(i) for i in row[9].strip(',').split(',')],
            'exonEnds': [int(i) for i in row[10].strip(',').split(',')],
            'score': int(row[11]),
            'name2': row[12],
            'cdsStartStat': row[13],
            'cdsEndStat': row[14],
            'exonFrames': row[15],
            'exons': exons,
        }


def make_transcript(transcript_json):
    """
    Make a Transcript form a JSON object.
    """

    transcript_name = transcript_json['name']
    if '.' in transcript_name:
        name, version = transcript_name.split('.')
    else:
        name, version = transcript_name, None

    exonlist = list(zip(transcript_json['exonStarts'], transcript_json['exonEnds']))

    exons = transcript_json['exons']
    exon_frames = transcript_json['exonFrames'].strip(',').split(',')
    if transcript_json['strand'] != '+':
        exons = reversed(exons)
        exon_frames = list(reversed(exon_frames))

    transcript = Transcript(
        name=name,
        version=int(version) if version is not None else None,
        gene=transcript_json['name2'],
        tx_position=Position(
            transcript_json['chrom'],
            transcript_json['txStart'],
            transcript_json['txEnd'],
            transcript_json['strand'] == '+'),
        cds_position=Position(
            transcript_json['chrom'],
            transcript_json['cdsStart'],
            transcript_json['cdsEnd'],
            transcript_json['strand'] == '+'),
        exon_frames=exon_frames,
        exonlist=exonlist)

    # exons = transcript_json['exons']
    # if not transcript.tx_position.is_forward_strand:
    #     exons = reversed(exons)

    for exon_number, (exon_start, exon_end) in enumerate(exons, 1):
        transcript.exons.append(
            Exon(transcript=transcript,
                 tx_position=Position(
                     transcript_json['chrom'],
                     exon_start,
                     exon_end,
                     transcript_json['strand'] == '+'),
                 exon_number=exon_number))

    return transcript


def read_transcripts(trans_file):
    """
    Read all transcripts in a RefGene file.
    """
    transcripts = {}
    for trans in map(make_transcript, read_refgene(trans_file)):
        transcripts[trans.name] = trans
        transcripts[trans.full_name] = trans
        # transcripts[trans.gene.name] = trans
        transcripts[trans.full_name, trans.tx_position] = trans

    return transcripts


def get_transcript_list(transcripts, chrom, offset):
    """
    Get the transcript list for annotation
    """
    transcript_list = []
    for trans in transcripts:
        transcript = transcripts[trans]
        if chrom == transcript.tx_position.chrom and \
                transcript.tx_position.chrom_start <= offset <= transcript.tx_position.chrom_stop and \
                transcript not in transcript_list:
            transcript_list.append(transcript)
    return transcript_list
