#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 11:22:53 2016

@author: sshank
"""

import unittest

from pairwise import get_sequences

class TestStringMethods(unittest.TestCase):

    def test_get_sequences(self):
        get_sequences('alpha')
        get_sequences('beta')

if __name__ == '__main__':
    unittest.main()
