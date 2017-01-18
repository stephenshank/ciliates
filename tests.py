#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 11:22:53 2016

@author: sshank
"""

import unittest

from tools import get_list
from tools import get_dictionary
from tools import subunit_dictionary


class TestDataIntegrity(unittest.TestCase):
    
    def setUp(self):
        self.nucleotide_records = get_list('InitialNucleotide')
        self.nucleotide_dictionary = get_dictionary('InitialNucleotide')
        self.protein_dictionary = get_dictionary('InitialProtein')

    def test_translation(self):
        for nucleotide_record in self.nucleotide_records:
            translated_sequence = nucleotide_record.seq.translate(table=6)
            protein_sequence = self.protein_dictionary[nucleotide_record.id].seq
            self.assertEqual(translated_sequence, protein_sequence)
    
    def test_subunits(self):
        for subunit in ['alpha', 'beta']:
            sub_dict = subunit_dictionary(subunit)
            for value in sub_dict.values():
                self.assertTrue(value in self.nucleotide_dictionary.keys())
                self.assertTrue(value in self.protein_dictionary.keys())
            

if __name__ == '__main__':
    unittest.main()
