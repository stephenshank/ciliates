#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 13:01:07 2016

Handles pairwise portion of pipeline. Gets pairwise alignments for all sequences
in either alpha or beta subunit, and computes summary statistics to be used
for clustering prior to multiple sequence alignment.

@author: sshank
"""

import os
import itertools as it
import subprocess

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO


def get_sequences(subunit):
    path = os.path.join('data', 'Ciliate_%s.fasta' % subunit)
    return list(SeqIO.parse(path, 'fasta'))


def make_sequence_files(subunit):
    print('Splitting main %s file into multiple pairwise files...' % subunit)
    sequences = get_sequences(subunit)
    for sequence1, sequence2 in it.combinations(sequences, 2):
        id1 = sequence1.id
        id2 = sequence2.id
        pairwise_filename = '%s__%s.fasta' % (id1, id2)
        pairwise_path = os.path.join('data', 'pairwise', pairwise_filename)
        with open(pairwise_path, 'w') as output_file:
            SeqIO.write([sequence1, sequence2], output_file, 'fasta')


def perform_alignments():
    all_pairwise_files = os.listdir(os.path.join('data','pairwise'))
    already_aligned = lambda file: 'ALIGNED' in file
    desired_pairwise_files = it.filterfalse(already_aligned, all_pairwise_files)
    for filename in desired_pairwise_files:
        title = filename[:-6]
        output_filename = title + '__ALIGNED.fasta'
        input_path = os.path.join('data', 'pairwise', filename)
        output_path = os.path.join('data', 'pairwise', output_filename)
        alignment_command = 'mafft %s > %s' % (input_path, output_path)
        subprocess.call(alignment_command, shell=True)


def get_metrics(first_sequence, second_sequence):
    first_array = np.array(list(str(first_sequence)), dtype='<U1')
    second_array = np.array(list(str(second_sequence)), dtype='<U1')
    gaps = (first_array == '-') | (second_array == '-')
    not_gaps = ~gaps
    percent_identity = sum((first_array==second_array) & not_gaps)/sum(not_gaps)
    total_length = len(first_array)
    fraction_of_gaps = sum(gaps)/total_length
    return percent_identity, fraction_of_gaps


def get_metrics_from_filename(filename):
    path = os.path.join('data', 'pairwise', filename)
    alignment = AlignIO.read(path, 'fasta')
    first_sequence = alignment[0].seq
    second_sequence = alignment[1].seq
    return get_metrics(first_sequence, second_sequence)


def get_all_metrics(subunit):
    sequence_ids = [sequence.id for sequence in get_sequences(subunit)]
    first_ids = []
    second_ids = []
    all_identity = []
    all_gaps = []
    for id1, id2 in it.combinations(sequence_ids, 2):
        filename = '%s__%s__ALIGNED.fasta' % (id1, id2)
        percent_identity, fraction_of_gaps = get_metrics_from_filename(filename)
        first_ids.append(id1)
        second_ids.append(id2)
        all_identity.append(percent_identity)
        all_gaps.append(fraction_of_gaps)
    metrics_df = pd.DataFrame({
        'id1': first_ids, 'id2': second_ids, 'identity': all_identity, 'gaps': all_gaps
    })
    metrics_df['origin'] = '%s pairwise' % subunit
    return metrics_df


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
            make_sequence_files('alpha')
        if subunit == 'beta' or subunit =='both':
            make_sequence_files('beta')
    
    align = args.align
    if align:
        perform_alignments()
