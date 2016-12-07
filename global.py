#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 10:29:36 2016

@author: sshank
"""

import os
import subprocess
import itertools as it

import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage, fcluster
import matplotlib.pyplot as plt
from Bio import SeqIO, AlignIO

import pairwise


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
    

def get_all_metrics(subunit, nclust, method, index):
    aligned_filename = get_aligned_filename(subunit, nclust, method, index)
    alignment_path = os.path.join('data', 'split', aligned_filename)
    alignment = AlignIO.read(alignment_path, 'fasta')
    gaps = []
    identities = []
    for record1, record2 in it.combinations(alignment, 2):
        identity, gap = pairwise.get_metrics(record1.seq, record2.seq)
        gaps.append(gap)
        identities.append(identity)
    origin = '%s clade %d' % (subunit, index)
    return pd.DataFrame({'gaps':gaps, 'identity':identities, 'origin':origin})


def create_replicated_clusters():
    get_cluster_assignments('alpha', 3, 'ward')
    get_cluster_assignments('beta', 2, 'ward')


if __name__ == '__main__':
    create_replicated_clusters()
