#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 10:29:36 2016

@author: sshank
"""

import os
import subprocess
import itertools as it

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from Bio import SeqIO, AlignIO
import seaborn as sns
import matplotlib.pyplot as plt

import pairwise


round_threshold = lambda threshold: str(round(100*threshold))

def get_aligned_filename(subunit, nclust, method, index):
    return '%s_%d_%s_%d_ALIGNED.fasta' % (subunit, nclust, method, index)


def get_cluster_assignments(subunit, nclust, method='average'):
    metrics = pairwise.get_all_metrics(subunit)
    Z = linkage(metrics.dist, method=method)
    assignments = fcluster(Z, nclust, criterion='maxclust')
    print(list(assignments))
    path = os.path.join('data', 'Ciliate_%s.fasta' % subunit)
    sequences = SeqIO.parse(path, 'fasta')
    sequence_clusters = [[] for i in range(nclust)]
    for sequence, assignment in zip(sequences, assignments):
        sequence_clusters[assignment-1].append(sequence)
    for i, cluster in enumerate(sequence_clusters):
        print('--- CLUSTER %d : %d members ---' % (i+1, len(cluster)))
        print([sequence.id for sequence in cluster])
        filename = '%s_%d_%s_%d.fasta' % (subunit, nclust, method, i)
        path = os.path.join('data', 'split', filename)
        with open(path, 'w') as output_file:
            SeqIO.write(cluster, output_file, 'fasta')
        aligned_filename = get_aligned_filename(subunit, nclust, method, i)
        aligned_path = os.path.join('data', 'split', aligned_filename)
        alignment_command = 'mafft %s > %s' % (path, aligned_path)
        subprocess.call(alignment_command, shell=True)
    

def get_all_metrics(subunit, threshold, index):
    rounded_threshold = round_threshold(threshold)
    filename = '%s__%d__ALIGNED.fasta' % (subunit, index)
    alignment_path = os.path.join('data', 'split', rounded_threshold, filename)
    alignment = AlignIO.read(alignment_path, 'fasta')
    gaps = []
    identities = []
    for record1, record2 in it.combinations(alignment, 2):
        identity, gap = pairwise.get_metrics(record1.seq, record2.seq)
        gaps.append(gap)
        identities.append(identity)
    origin = '%s clade %d' % (subunit, index)
    return pd.DataFrame({'gaps':gaps, 'identity':identities, 'origin':origin})


def single_linkage_clustering(subunit, threshold):
    sequences = [sequence.id for sequence in pairwise.get_sequences(subunit)]
    number_of_sequences = len(sequences)
    metrics = pairwise.get_all_metrics(subunit)
    metrics['above_threshold'] = metrics.identity > threshold
    current_assignment = np.arange(number_of_sequences)
    for _, row in metrics[metrics['above_threshold']].iterrows():
        index1 = sequences.index(row.id1)
        index2 = sequences.index(row.id2)
        assignment1 = current_assignment[index1]
        assignment2 = current_assignment[index2]        
        current_assignment[current_assignment==assignment1] = assignment2
    return current_assignment
    

def cluster_and_align(subunit, threshold):
    rounded_threshold = round_threshold(threshold)
    data_directory = os.path.join('data', 'split', rounded_threshold)
    if not os.path.exists(data_directory):
        os.makedirs(data_directory)
    clusters = single_linkage_clustering(subunit, threshold)
    sequences = list(pairwise.get_sequences(subunit))
    for i, cluster in enumerate(set(clusters)):
        filename = '%s__%d.fasta' % (subunit, i)
        unaligned_path = os.path.join(data_directory, filename)
        indices = np.arange(len(clusters))[clusters==cluster]
        if len(indices) > 1:
            cluster_sequences = []
            for index in indices:
                cluster_sequences.append(sequences[index])
            with open(unaligned_path, 'w') as output_file:
                SeqIO.write(cluster_sequences, output_file, 'fasta')
            aligned_filename = '%s__%d__ALIGNED.fasta' % (subunit, i)
            aligned_path = os.path.join(data_directory, aligned_filename)
            alignment_command = 'mafft %s > %s' % (unaligned_path, aligned_path)
            subprocess.call(alignment_command, shell=True)


def get_stats(subunit, index):
    metrics = get_all_metrics(subunit, index)
    sns.distplot(metrics.identity.dropna())
    plt.show()


def create_replicated_clusters():
    get_cluster_assignments('alpha', 3, 'ward')
    get_cluster_assignments('beta', 2, 'ward')


if __name__ == '__main__':
    cluster_and_align('beta', .6)
