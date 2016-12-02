#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 13:01:07 2016

Handles pairwise portion of pipeline. Gets pairwise alignments for all
in either alpha or beta subunit, and computes summary statistics to be used
for clustering prior to multiple sequence alignment.

@author: sshank
"""

import os
import itertools as it

from Bio import SeqIO


def make_pairwise_sequence_files(subunit):
    print('Splitting main %s file into multiple pairwise files...' % subunit)
    path = os.path.join('data', 'Ciliate_%s.fasta' % subunit)
    sequences = SeqIO.parse(path, 'fasta')
    for sequence1, sequence2 in it.combinations(sequences, 2):
        id1 = sequence1.id
        id2 = sequence2.id
        pairwise_filename = '%s_%s.fasta' % (id1, id2)
        pairwise_path = os.path.join('data', 'pairwise', pairwise_filename)
        with open(pairwise_path, 'w') as output_file:
            SeqIO.write([sequence1, sequence2], output_file, 'fasta')


def perform_pairwise_alignments():
    print('Aligning all files in data/pairwise...')


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-s', '--subunit', metavar='SUBUNIT', dest='subunit',
                        nargs='?', const='both',
                        help='Subunit to process (alpha or beta)')
    parser.add_argument('-a', '--align', dest='align', nargs='?', const=True,
                        help='Perform pairwise alignments.', default=False)
    args = parser.parse_args()
    
    subunit = args.subunit
    if subunit:
        if subunit == 'alpha' or subunit =='both':
            make_pairwise_sequence_files('alpha')
        if subunit == 'beta' or subunit =='both':
            make_pairwise_sequence_files('beta')
    
    align = args.align
    if align:
        perform_pairwise_alignments()
