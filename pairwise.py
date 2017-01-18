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

import pandas as pd
from Bio import SeqIO
from Bio import AlignIO

from tools import get_list
from tools import get_metrics
from tools import remove_badchars


def get_pairwise_filename(record1, record2, aligned=''):
    id1 = remove_badchars(record1.id)
    id2 = remove_badchars(record2.id)
    if aligned:
        aligned = '__ALIGNED'
    filename = '%s__%s%s.fasta' % (id1, id2, aligned)
    return filename


def make_sequence_files():
    print('Splitting protein sequences into pairwise files...')
    sequences = get_list('InitialProtein')
    for sequence1, sequence2 in it.combinations(sequences, 2):
        pairwise_filename = get_pairwise_filename(sequence1, sequence2)
        pairwise_path = os.path.join('data', 'pairwise', pairwise_filename)
        with open(pairwise_path, 'w') as output_file:
            SeqIO.write([sequence1, sequence2], output_file, 'fasta')


def perform_alignments():
    print('Performing all pairwise alignments...')
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


def get_metrics_from_filename(filename):
    path = os.path.join('data', 'pairwise', filename)
    alignment = AlignIO.read(path, 'fasta')
    first_sequence = alignment[0].seq
    second_sequence = alignment[1].seq
    first_length = len(str(first_sequence).replace('-', ''))
    second_length = len(str(second_sequence).replace('-', ''))
    percent_identity, gap_fraction = get_metrics(first_sequence, second_sequence)
    pairwise_alignment_metrics = {
        'percent_identity': percent_identity,
        'gap_fraction': gap_fraction,
        'first_length': first_length,
        'second_length': second_length
    }
    return pairwise_alignment_metrics


def get_all_metrics():
    sequences = get_list('InitialProtein')
    first_ids = []
    second_ids = []
    first_lengths = []
    second_lengths = []
    all_identity = []
    all_gaps = []
    for sequence1, sequence2 in it.combinations(sequences, 2):
        filename = get_pairwise_filename(sequence1, sequence2, aligned=True)
        alignment_metrics = get_metrics_from_filename(filename)
        first_ids.append(sequence1.id)
        second_ids.append(sequence2.id)
        first_lengths.append(alignment_metrics['first_length'])
        second_lengths.append(alignment_metrics['second_length'])
        all_identity.append(alignment_metrics['percent_identity'])
        all_gaps.append(alignment_metrics['gap_fraction'])
    metrics_df = pd.DataFrame({
        'id1': first_ids,
        'id2': second_ids,
        'length1': first_lengths,
        'length2': second_lengths,
        'identity': all_identity,
        'gaps': all_gaps
    })
    metrics_df['status'] = 'okay'
    isoforms = (metrics_df.identity == 1) & (metrics_df.gaps > 0)
    duplicates = (metrics_df.identity == 1) & (metrics_df.gaps == 0)
    metrics_df.loc[isoforms, 'status'] = 'isoform'
    metrics_df.loc[duplicates, 'status'] = 'duplicate'
    return metrics_df


if __name__ == '__main__':
    pairwise_dir = os.path.join('data', 'pairwise')
    if not os.path.exists(pairwise_dir):
        os.makedirs(pairwise_dir)
    make_sequence_files()
    perform_alignments()
