#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 11:42:29 2017

@author: sshank
"""

import os
import json

import numpy as np
from Bio import SeqIO


data_dir = 'data'
initial_dir = 'initial'

filenames = {
    'InitialNucleotide': 'nucleotide.fasta',
    'InitialProtein': 'protein.fasta',
}

directories = {
    'InitialNucleotide': os.path.join(data_dir, initial_dir),
    'InitialProtein': os.path.join(data_dir, initial_dir),
}

remove_badchars = lambda string: string.replace('|', '_').replace('.', '_')


def get_path(key):
    directory = directories[key]
    filename = filenames[key]
    return os.path.join(directory, filename)

    
def get_generator(key):
    path = get_path(key)
    return SeqIO.parse(path, 'fasta')


def get_list(key):
    return list(get_generator(key))


def get_dictionary(key):
    return SeqIO.to_dict(get_generator(key))

    
def get_metrics(first_sequence, second_sequence):
    first_array = np.array(list(str(first_sequence)), dtype='<U1')
    second_array = np.array(list(str(second_sequence)), dtype='<U1')
    gaps = (first_array == '-') | (second_array == '-')
    not_gaps = ~gaps
    percent_identity = sum((first_array==second_array) & not_gaps)/sum(not_gaps)
    total_length = len(first_array)
    fraction_of_gaps = sum(gaps)/total_length
    return percent_identity, fraction_of_gaps


def check_sequence_identity(first_sequence, second_sequence):
    first_string = str(first_sequence).upper().replace('\n', '')
    second_string = str(second_sequence).upper().replace('\n', '')
    if len(first_string) != len(second_string):
        return False
    compare_characters = lambda i: first_string[i] == second_string[i]
    return all([compare_characters(i) for i in range(len(first_sequence))])


def subunit_dictionary(subunit):
    path = os.path.join('data', 'initial', '%s.json' % subunit)
    with open(path, 'r') as file:
        subunit_dict = json.load(file)
    return subunit_dict
