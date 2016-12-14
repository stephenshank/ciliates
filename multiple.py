#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 10:29:36 2016

@author: sshank
"""

import os
import subprocess
import itertools as it
import pickle

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


def directory_from_threshold(threshold):
    rounded_identity = round_threshold(threshold['identity'])
    rounded_gap = round_threshold(threshold['gaps'])
    directory_name = '%s_%s' % (rounded_identity, rounded_gap)
    return directory_name


def get_all_metrics(subunit, threshold, index):
    directory_name = directory_from_threshold(threshold)
    filename = '%s__%d__ALIGNED.fasta' % (subunit, index)
    alignment_path = os.path.join('data', 'split', directory_name, filename)
    alignment = AlignIO.read(alignment_path, 'fasta')
    gaps = []
    identities = []
    for record1, record2 in it.combinations(alignment, 2):
        identity, gap = pairwise.get_metrics(record1.seq, record2.seq)
        gaps.append(gap)
        identities.append(identity)
    origin = '%s clade %d' % (subunit, index)
    return pd.DataFrame({'gaps':gaps, 'identity':identities, 'origin':origin})


def single_linkage_clustering(subunit, threshold, variable='both'):
    sequences = [sequence.id for sequence in pairwise.get_sequences(subunit)]
    number_of_sequences = len(sequences)
    metrics = pairwise.get_all_metrics(subunit)
    if variable == 'both':
        correct_identity = metrics.identity > threshold['identity']
        correct_gaps = metrics.gaps < threshold['gaps']
        metrics['cluster_variable'] = correct_identity & correct_gaps
    elif variable == 'identity':
        metrics['cluster_variable'] = metrics.identity > threshold
    current_assignment = np.arange(number_of_sequences)
    for _, row in metrics[metrics['cluster_variable']].iterrows():
        index1 = sequences.index(row.id1)
        index2 = sequences.index(row.id2)
        assignment1 = current_assignment[index1]
        assignment2 = current_assignment[index2]        
        current_assignment[current_assignment==assignment1] = assignment2
    return current_assignment

    
def cluster_and_align(subunit, threshold):
    directory_name = directory_from_threshold(threshold)
    directory_path = os.path.join('data', 'split', directory_name)
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
    clusters = single_linkage_clustering(subunit, threshold)
    sequences = list(pairwise.get_sequences(subunit))
    for i, cluster in enumerate(set(clusters)):
        filename = '%s__%d.fasta' % (subunit, i)
        unaligned_path = os.path.join(directory_path, filename)
        indices = np.arange(len(clusters))[clusters==cluster]
        if len(indices) > 1:
            cluster_sequences = []
            for index in indices:
                cluster_sequences.append(sequences[index])
            with open(unaligned_path, 'w') as output_file:
                SeqIO.write(cluster_sequences, output_file, 'fasta')
            aligned_filename = '%s__%d__ALIGNED.fasta' % (subunit, i)
            aligned_path = os.path.join(directory_path, aligned_filename)
            alignment_command = 'mafft %s > %s' % (unaligned_path, aligned_path)
            subprocess.call(alignment_command, shell=True)


def get_stats(subunit, index):
    metrics = get_all_metrics(subunit, index)
    sns.distplot(metrics.identity.dropna())
    plt.show()

    
def make_and_plot(subunit, threshold, show=False):
    cluster_and_align(subunit, threshold)
    directory_name = directory_from_threshold(threshold)
    data_directory = os.path.join('data', 'split', directory_name)
    valid_file = lambda file: 'ALIGNED' in file and subunit in file
    files = [file for file in os.listdir(data_directory) if valid_file(file)]
    for file in files:
        path = os.path.join(data_directory, file)
        alignment = AlignIO.read(path, 'fasta')
        desired_sequences = [sequence.id for sequence in alignment if 'Contig' in sequence.id]
        contains_desired_sequence = any(desired_sequences)
        enough_members = len(alignment) > 2
        if contains_desired_sequence and enough_members:
            _, index, _ = file.split('__')
            index = int(index)
            print(subunit, index, ': contains', end='')
            print(', '.join(desired_sequences), len(alignment), 'total members')
            metrics = get_all_metrics(subunit, threshold, index)
            fig, axs = plt.subplots(1, 2, figsize=(10, 5))
            metrics.identity.dropna().plot(kind='hist', ax=axs[0])
            axs[0].set_title('Percent identity')
            axs[0].set_xlim([0, 1])
            metrics.gaps.dropna().plot(kind='hist', ax=axs[1])
            axs[1].set_title('Gap proportion')
            axs[1].set_xlim([0, 1])
            plt.show()


def fasta_renamer(filename, threshold):
    threshold_directory = directory_from_threshold(threshold)
    directory = os.path.join('data', 'split', threshold_directory)
    path = os.path.join(directory, filename)
    names = {}
    alignment = AlignIO.read(path, 'fasta')
    for i, sequence in enumerate(alignment):
        names[i] = sequence.id
        sequence.id = str(i)
    phylip_filename = filename.split('.')[0] + '.phylip'
    phylip_path = os.path.join(directory, phylip_filename)
    with open(phylip_path, 'w') as phylip_file:
        AlignIO.write(alignment, phylip_file, 'phylip')
    pickle_filename = filename.split('.')[0] + '.pkl'
    pickle_path = os.path.join(directory, pickle_filename)
    with open(pickle_path, 'wb') as pickle_file:
        pickle.dump(names, pickle_file)


def fasta_to_phylip(threshold):
    threshold_directory = directory_from_threshold(threshold)
    directory = os.path.join('data', 'split', threshold_directory)
    is_alignment = lambda file: 'ALIGNED' in file and file[-5:] == 'fasta'
    files = [file for file in os.listdir(directory) if is_alignment(file)]
    for filename in files:
        fasta_renamer(filename, threshold)    


if __name__ == '__main__':
    threshold = {'identity':.6, 'gaps':.4}
    cluster_and_align('alpha', threshold)
    cluster_and_align('beta', threshold)
    fasta_to_phylip(threshold)
