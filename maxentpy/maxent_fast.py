'''
maxent_fast.py
Calculate splice site strength
Modified from MaxEntScan perl scripts developed by Gene Yeo and Christopher
Burge
Yeo G and Burge C. Maximum entropy modeling of short sequence motifs with
applications to RNA splicing signals. Journal of Computational Biology,
2004, 11:377-94.
'''

import sys
import math
from maxentpy._hashseq import hashseq
import msgpack
try:
    from string import maketrans
except ImportError:
    maketrans = str.maketrans
import os.path


__all__ = ['score5', 'score3', 'load_matrix']

dir_path = os.path.dirname(os.path.abspath(__file__))

bgd_5 = {'A': 0.27, 'C': 0.23, 'G': 0.23, 'T': 0.27}
cons1_5 = {'A': 0.004, 'C': 0.0032, 'G': 0.9896, 'T': 0.0032}
cons2_5 = {'A': 0.0034, 'C': 0.0039, 'G': 0.0042, 'T': 0.9884}

bgd_3 = {'A': 0.27, 'C': 0.23, 'G': 0.23, 'T': 0.27}
cons1_3 = {'A': 0.9903, 'C': 0.0032, 'G': 0.0034, 'T': 0.0030}
cons2_3 = {'A': 0.0027, 'C': 0.0037, 'G': 0.9905, 'T': 0.0030}


def score5(fa, matrix=None):
    '''
    Calculate 5' splice site strength
    (exon)XXX|XXXXXX(intron)
              **
    >>> round(score5('cagGTAAGT'), 2)
    10.86
    >>> round(score5('gagGTAAGT'), 2)
    11.08
    >>> round(score5('taaATAAGT'), 2)
    -0.12
    >>> matrix = load_matrix(5)
    >>> round(score5('cagGTAAGT', matrix=matrix), 2)
    10.86
    '''
    # check length of fa
    if len(fa) != 9:
        sys.exit('Wrong length of fa!')
    # check matrix
    if not matrix:
        matrix = load_matrix(5)
    # for key elements
    key = fa[3:5].upper()
    score = cons1_5[key[0]] * cons2_5[key[1]] / (bgd_5[key[0]] * bgd_5[key[1]])
    # for rest elements
    rest = (fa[:3] + fa[5:]).upper()
    rest_score = matrix[rest]
    # final score
    return math.log(score * rest_score, 2)


def score3(fa, matrix=None):
    '''
    Calculate 3' splice site strength
    (intron)XXXXXXXXXXXXXXXXXXXX|XXX(exon)
                              **
    >>> round(score3('ttccaaacgaacttttgtAGgga'), 2)
    2.89
    >>> round(score3('tgtctttttctgtgtggcAGtgg'), 2)
    8.19
    >>> round(score3('ttctctcttcagacttatAGcaa'), 2)
    -0.08
    >>> matrix = load_matrix(3)
    >>> round(score3('ttccaaacgaacttttgtAGgga', matrix=matrix), 2)
    2.89
    '''
    # check length of fa
    if len(fa) != 23:
        sys.exit('Wrong length of fa!')
    # check matrix
    if not matrix:
        matrix = load_matrix(3)
    # for key elements
    key = fa[18:20].upper()
    score = cons1_3[key[0]] * cons2_3[key[1]] / (bgd_3[key[0]] * bgd_3[key[1]])
    # for rest elements
    rest = (fa[:18] + fa[20:]).upper()
    rest_score = 1
    rest_score *= matrix[0][hashseq(rest[:7])]
    rest_score *= matrix[1][hashseq(rest[7:14])]
    rest_score *= matrix[2][hashseq(rest[14:])]
    rest_score *= matrix[3][hashseq(rest[4:11])]
    rest_score *= matrix[4][hashseq(rest[11:18])]
    rest_score /= matrix[5][hashseq(rest[4:7])]
    rest_score /= matrix[6][hashseq(rest[7:11])]
    rest_score /= matrix[7][hashseq(rest[11:14])]
    rest_score /= matrix[8][hashseq(rest[14:18])]
    # final score
    return math.log(score * rest_score, 2)


def load_matrix(d):
    matrix_f = os.path.join(dir_path, 'data/matrix%d.msg' % d)
    with open(matrix_f, 'rb') as f:
        matrix = msgpack.unpackb(f.read(), encoding='utf-8')
    return matrix


if __name__ == '__main__':
    import doctest
    doctest.testmod()
