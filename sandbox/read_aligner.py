#! /usr/bin/env python
"""
Error correct reads based on a counting hash from a diginorm step.
Output sequences will be put in @@@.

% python scripts/error-correct-pass2 <counting.ct> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

from __future__ import print_function

import sys
import screed
import os
import khmer
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_loader

from khmer.khmer_args import build_counting_args
from khmer.khmer_args import add_loadhash_args

###

DEFAULT_COVERAGE = 20
DEFAULT_MAX_ERROR_REGION = 40


def main():
    parser = build_counting_args()
    parser.add_argument("--trusted-cov", dest="trusted_cov", type=int, default=2)
    parser.add_argument("--theta", type=float, default=1.0)
    parser.add_argument("input_table")
    parser.add_argument("input_filenames", nargs="+")
    add_loadhash_args(parser)

    args = parser.parse_args()

    counting_ht = args.input_table
    infiles = args.input_filenames

    print('file with ht: %s' % counting_ht, file=sys.stderr)

    print('loading hashtable', file=sys.stderr)
    ht = khmer.load_counting_hash(counting_ht)
    K = ht.ksize()

    aligner = khmer.new_readaligner(ht, args.trusted_cov, args.theta) # counting hash, trusted kmer coverage cutoff, bits theta (threshold value for terminating unproductive alignemnts)

    ### the filtering loop
    for infile in infiles:
        print('aligning', infile, file=sys.stderr)
        for n, record in enumerate(screed.open(infile)):

            name = record['name']
            seq = record['sequence'].upper()
            print(name, file=sys.stderr)
            print(seq, file=sys.stderr)

            score, graph_alignment, read_alignment, truncated = aligner.align(seq)
            print(score, file=sys.stderr)
            print(graph_alignment, file=sys.stderr)
            print(read_alignment, file=sys.stderr)
            print(truncated, file=sys.stderr)
            print(">{0}\n{1}".format(name, graph_alignment))

if __name__ == '__main__':
    main()
