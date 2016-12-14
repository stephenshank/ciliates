#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 15:07:50 2016

@author: sshank
"""

import os
import subprocess
import pickle

import ete3


def write_newick(subunit, index):
    directory = os.path.join('data', 'split')
    pickle_filename = '%s__%d__ALIGNED.pkl' % (subunit, index)
    tree_filename = '%s__%d__ALIGNED.phylip_phyml_tree.txt' % (subunit, index)
    pickle_path = os.path.join(directory, pickle_filename)
    tree_path = os.path.join(directory, tree_filename)
    with open(pickle_path, 'rb') as pickle_file:
        names = pickle.load(pickle_file)
    tree = ete3.Tree(tree_path)
    for leaf in tree.get_leaves():
        leaf.name = names[int(leaf.name)]
    new_tree_filename = '%s__%d.new' % (subunit, index)
    new_tree_path = os.path.join('data', 'trees', new_tree_filename)
    tree.write(outfile=new_tree_path)
    

if __name__ == '__main__':
    subprocess.call(['bash', 'prottest.sh'])
    subprocess.call(['bash', 'phyml.sh'])
    write_newick('alpha', 12)
    write_newick('alpha', 18)
    write_newick('alpha', 21)
    write_newick('alpha', 25)
    write_newick('beta', 19)
