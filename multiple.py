#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 10:29:36 2016

@author: sshank
"""

import os
import subprocess
import itertools as it
import json
import re

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
import matplotlib.pyplot as plt

import pairwise
from tools import get_list
from tools import get_dictionary
from tools import subunit_dictionary


def alignment_path(cluster, extension='fasta'):
    filename = 'cluster_%d__ALIGNED.%s' % (cluster, extension)
    return os.path.join('data', 'clusters', filename)


def alignment_files(extension='fasta'):
    alignment_directory = os.path.join('data', 'clusters')
    all_files = os.listdir(alignment_directory)
    is_alignment = lambda filename: 'ALIGNED' in filename
    is_format = lambda filename: filename.split('.')[-1] == extension
    is_desired = lambda filename: is_alignment(filename) and is_format(filename)
    return filter(is_desired, all_files)
    

def get_all_metrics(cluster):
    alignment = AlignIO.read(alignment_path(cluster), 'fasta')
    gaps = []
    identities = []
    for record1, record2 in it.combinations(alignment, 2):
        identity, gap = pairwise.get_metrics(record1.seq, record2.seq)
        gaps.append(gap)
        identities.append(identity)
    return pd.DataFrame({'gaps':gaps, 'identity':identities})


def single_linkage_clustering(threshold, variable='both'):
    sequences = get_list('InitialProtein')
    sequence_ids = [sequence.id for sequence in sequences]
    number_of_sequences = len(sequences)
    metrics = pairwise.get_all_metrics()
    if variable == 'both':
        correct_identity = metrics.identity > threshold['identity']
        correct_gaps = metrics.gaps < threshold['gaps']
        metrics['cluster_variable'] = correct_identity & correct_gaps
    elif variable == 'identity':
        metrics['cluster_variable'] = metrics.identity > threshold
    current_assignment = np.arange(number_of_sequences)
    for _, row in metrics[metrics['cluster_variable']].iterrows():
        index1 = sequence_ids.index(row.id1)
        index2 = sequence_ids.index(row.id2)
        assignment1 = current_assignment[index1]
        assignment2 = current_assignment[index2]        
        current_assignment[current_assignment==assignment1] = assignment2
    return current_assignment

    
def cluster_and_align(threshold):
    directory_path = os.path.join('data', 'clusters')
    clusters = single_linkage_clustering(threshold)
    sequences = get_list('InitialProtein')
    count = 0
    for cluster in set(clusters):
        indices = np.arange(len(clusters))[clusters==cluster]
        if len(indices) > 1:
            filename = 'cluster_%d.fasta' % count
            unaligned_path = os.path.join(directory_path, filename)
            cluster_sequences = []
            for index in indices:
                cluster_sequences.append(sequences[index])
            with open(unaligned_path, 'w') as output_file:
                SeqIO.write(cluster_sequences, output_file, 'fasta')
            aligned_filename = 'cluster_%d__ALIGNED.fasta' % count
            aligned_path = os.path.join(directory_path, aligned_filename)
            alignment_command = 'mafft %s > %s' % (unaligned_path, aligned_path)
            subprocess.call(alignment_command, shell=True)
            count += 1

    
def plot(cluster):
    metrics = get_all_metrics(cluster)
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    metrics.identity.dropna().plot(kind='hist', ax=axs[0])
    axs[0].set_title('Percent identity')
    axs[0].set_xlim([0, 1])
    metrics.gaps.dropna().plot(kind='hist', ax=axs[1])
    axs[1].set_title('Gap proportion')
    axs[1].set_xlim([0, 1])
    plt.show()


def cluster_info():
    alpha_ids = list(subunit_dictionary('alpha').values())
    beta_ids = list(subunit_dictionary('beta').values())
    i = 0
    total_size = 0
    while os.path.exists(alignment_path(i)):
        alignment = AlignIO.read(alignment_path(i), 'fasta')
        cluster_size = len(alignment)
        print('--- CLUSTER %2d ---' % i)
        print('Number of sequences: %d' % cluster_size)
        contains_alphas = False
        contains_betas = False
        for record in alignment:
            if record.id in alpha_ids:
                contains_alphas = True
            if record.id in beta_ids:
                contains_betas = True
        print('Contains alpha sequences: %s' % contains_alphas)
        print('Contains beta sequences: %s' % contains_betas)
        total_size += cluster_size
        i += 1
    print('Total size:', total_size)


def fasta_renamer(filename):
    directory = os.path.join('data', 'clusters')
    path = os.path.join(directory, filename)
    names = {}
    alignment = AlignIO.read(path, 'fasta')
    for i, record in enumerate(alignment):
        names[i] = record.id
        record.id = str(i)
    phylip_filename = filename.split('.')[0] + '.phylip'
    phylip_path = os.path.join(directory, phylip_filename)
    with open(phylip_path, 'w') as phylip_file:
        AlignIO.write(alignment, phylip_file, 'phylip')
    json_filename = filename.split('.')[0] + '.json'
    json_path = os.path.join(directory, json_filename)
    with open(json_path, 'w') as json_file:
        json.dump(names, json_file)


def fasta_to_phylip():
    directory = os.path.join('data', 'clusters')
    is_alignment = lambda file: 'ALIGNED' in file and file[-5:] == 'fasta'
    files = filter(is_alignment, os.listdir(directory))
    for filename in files:
        fasta_renamer(filename)


def run_prottest():
    prottest_loc = '$HOME/Software/prottest3/dist/prottest-3.4.2.jar'
    args = '-all-matrices -all-distributions'
    for file in alignment_files():
        path = os.path.join('data', 'clusters', file)
        alignment = AlignIO.read(path, 'fasta')
        if len(alignment) < 3:
            continue
        cluster = int(file.replace('.', '_').split('_')[1])
        output_filename = 'prottest_%d.txt' % cluster
        output_path = os.path.join('data', 'clusters', output_filename)
        info = (prottest_loc, alignment_path(cluster), output_path, args)
        command = 'java -jar %s -i %s -o %s %s' % info
        subprocess.call(command, shell=True)


def run_phyml():
    prottest_output_dir = os.path.join('data', 'clusters')
    prottest_regex = 'prottest_\d+\.txt'
    prottest_filenames = filter(
        lambda file: re.match(prottest_regex, file),
        os.listdir(prottest_output_dir)
    )
    clusters = []
    for prottest_filename in prottest_filenames:
        prottest_path = os.path.join(prottest_output_dir, prottest_filename)
        with open(prottest_path, 'r') as prottest_file:
            line = prottest_file.readlines()[421]
        model_list = line.split()[-1].split('+')
        cluster = int(prottest_filename.split('_')[-1].split('.')[0])
        clusters.append(cluster)
        model = model_list[0]
        gamma_supported = len(model_list) > 1
        gamma_arg = '-c 4 ' if gamma_supported else ''
        info = (alignment_path(cluster, extension='phylip'), model, gamma_arg)
        command = 'phyml3 -i %s -d aa -m %s %s-b 100' % info
        subprocess.call(command, shell=True)
    json_filename = 'clusters.json'
    json_path = os.path.join(prottest_output_dir, json_filename)
    with open(json_path, 'w') as json_file:
        json.dump(clusters, json_file)


def write_codon_alignments():
    directory = os.path.join('data', 'clusters')
    clusters_path = os.path.join(directory, 'clusters.json')
    with open(clusters_path, 'r') as json_file:
        clusters = json.load(json_file)
    all_nucleotide_records = get_dictionary('InitialNucleotide')
    for cluster in clusters:
        nucleotide_records = []
        amino_acid_records = []
        id_map_filename = 'cluster_%d__ALIGNED.json' % cluster
        id_map_path = os.path.join(directory, id_map_filename)
        with open(id_map_path, 'r') as json_file:
            id_map = json.load(json_file)
        phylip_filename = 'cluster_%d__ALIGNED.phylip' % cluster
        phylip_path = os.path.join(directory, phylip_filename)
        phylip_records = SeqIO.parse(phylip_path, 'phylip')
        for record in phylip_records:
            amino_acid_records.append(record)
            nucleotide_id = id_map[record.id]
            nucleotide_record = all_nucleotide_records[nucleotide_id]
            nucleotide_record.id = record.id
            nucleotide_record.description = ''
            nucleotide_records.append(nucleotide_record)
        nucleotide_filename = 'cluster_%d__nucleotide.fasta' % cluster
        nucleotide_path = os.path.join(directory, nucleotide_filename)
        SeqIO.write(nucleotide_records, nucleotide_path, 'fasta')
        protein_filename = 'cluster_%d_NumID__ALIGNED.fasta' % cluster
        protein_path = os.path.join(directory, protein_filename)
        SeqIO.write(amino_acid_records, protein_path, 'fasta')
        codon_filename = 'cluster_%d_CODON__ALIGNED.fasta' % cluster
        codon_path = os.path.join(directory, codon_filename)
        info = (protein_path, nucleotide_path, codon_path)
        command = 'pal2nal %s %s -codontable 6 -output fasta > %s' % info
        subprocess.call(command, shell=True)


if __name__ == '__main__':
    threshold = {'identity':.6, 'gaps':.4}
    clusters_dir = os.path.join('data', 'clusters')
    if not os.path.exists(clusters_dir):
        os.makedirs(clusters_dir)
    cluster_and_align(threshold)
    fasta_to_phylip()
    cluster_info()
    run_prottest()
    run_phyml()
    write_codon_alignments()
