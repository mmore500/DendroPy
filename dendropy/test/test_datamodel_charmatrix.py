#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Tests character sequence map.
"""

import unittest
import dendropy
from dendropy.datamodel import charmatrixmodel
from dendropy.test.support import dendropytest

class TaxonCharacterMatrixBasicCrud(dendropytest.ExtendedTestCase):

    def get_taxon_namespace(self, ntax):
        taxon_namespace = dendropy.TaxonNamespace()
        for i in range(ntax):
            label = "T{}".format(i)
            t = taxon_namespace.require_taxon(label=label)
        return taxon_namespace

    def test_setitem_by_taxon(self):
        tns = self.get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        self.assertEqual(len(char_matrix), 0)
        seqs = [
                "abcd",
                [1,2,3,4,],
                ["a", "b", "c", "d",]
                ]
        assert len(seqs) == len(tns)
        for idx, taxon in enumerate(tns):
            self.assertFalse(taxon in char_matrix)
            self.assertNotIn(taxon, char_matrix)
            char_matrix[taxon] = seqs[idx]
        self.assertEqual(len(char_matrix._taxon_sequence_map), len(tns))
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        for idx, taxon in enumerate(tns):
            self.assertTrue(taxon in char_matrix)
            self.assertIn(taxon, char_matrix)
            self.assertTrue(isinstance(char_matrix[taxon], charmatrixmodel.CharacterSequence))
            self.assertEqual(len(char_matrix[taxon]), len(seqs[idx]))
            for c1, c2 in zip(char_matrix[taxon], seqs[idx]):
                self.assertEqual(c1, c2)

    def test_setitem_by_taxon_idx(self):
        tns = self.get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        self.assertEqual(len(char_matrix), 0)
        seqs = [
                "abcd",
                [1,2,3,4,],
                ["a", "b", "c", "d",]
                ]
        assert len(seqs) == len(tns)
        for idx, taxon in enumerate(tns):
            self.assertFalse(taxon in char_matrix)
            self.assertNotIn(taxon, char_matrix)
            char_matrix[idx] = seqs[idx]
        self.assertEqual(len(char_matrix._taxon_sequence_map), len(tns))
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        for idx, taxon in enumerate(tns):
            self.assertTrue(taxon in char_matrix)
            self.assertIn(taxon, char_matrix)
            self.assertTrue(isinstance(char_matrix[taxon], charmatrixmodel.CharacterSequence))
            self.assertEqual(len(char_matrix[taxon]), len(seqs[idx]))
            for c1, c2 in zip(char_matrix[taxon], seqs[idx]):
                self.assertEqual(c1, c2)

    def test_setitem_by_taxon_label(self):
        tns = self.get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        self.assertEqual(len(char_matrix), 0)
        seqs = [
                "abcd",
                [1,2,3,4,],
                ["a", "b", "c", "d",]
                ]
        assert len(seqs) == len(tns)
        for idx, taxon in enumerate(tns):
            self.assertFalse(taxon in char_matrix)
            self.assertNotIn(taxon, char_matrix)
            char_matrix[taxon.label] = seqs[idx]
        self.assertEqual(len(char_matrix._taxon_sequence_map), len(tns))
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        for idx, taxon in enumerate(tns):
            self.assertTrue(taxon in char_matrix)
            self.assertIn(taxon, char_matrix)
            self.assertTrue(isinstance(char_matrix[taxon], charmatrixmodel.CharacterSequence))
            self.assertEqual(len(char_matrix[taxon]), len(seqs[idx]))
            for c1, c2 in zip(char_matrix[taxon], seqs[idx]):
                self.assertEqual(c1, c2)

    def test_setitem_by_taxon_not_in_namespace(self):
        tns = self.get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix()
        t = tns[0]
        seq = ["a", "b"]
        with self.assertRaises(ValueError):
            char_matrix[t] = seq
        char_matrix.taxon_namespace.add_taxon(t)
        char_matrix[t] = seq
        self.assertEqual(len(char_matrix), 1)
        self.assertIn(t, char_matrix)
        self.assertEqual(len(char_matrix[t]), len(seq))
        self.assertTrue(isinstance(char_matrix[t], charmatrixmodel.CharacterSequence))
        for c1, c2 in zip(char_matrix[t], seq):
            self.assertEqual(c1, c2)

    def test_setitem_by_idx_not_in_namespace(self):
        tns = self.get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix()
        with self.assertRaises(IndexError):
            char_matrix[len(tns)] = []

    def test_setitem_by_idx_not_in_namespace(self):
        tns = self.get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix()
        with self.assertRaises(KeyError):
            char_matrix[tns[0].label] = []

    def test_multi_setitem(self):
        tns = self.get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        self.assertEqual(len(char_matrix), 0)
        seqs = [
                "abcd",
                [1,2,3,4,],
                ["a", "b", "c", "d",]
                ]
        t = tns[0]
        for seq in seqs:
            char_matrix[t] = seq
        for taxon in tns:
            if taxon is t:
                self.assertIn(taxon, char_matrix)
            else:
                self.assertNotIn(taxon, char_matrix)
        seq = seqs[-1]
        self.assertEqual(len(char_matrix), 1)
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        self.assertEqual(len(char_matrix[0]), len(seq))
        self.assertTrue(isinstance(char_matrix[0], charmatrixmodel.CharacterSequence))
        for c1, c2 in zip(char_matrix[0], seq):
            self.assertEqual(c1, c2)
        for c1, c2 in zip(char_matrix[0], seqs[1]):
            self.assertNotEqual(c1, c2)

    def test_delitem(self):
        tns = self.get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        self.assertEqual(len(char_matrix), 0)
        seqs = [
                "abcd",
                [1,2,3,4,],
                ["a", "b", "c", "d",]
                ]
        assert len(seqs) == len(tns)
        for idx, taxon in enumerate(tns):
            self.assertFalse(taxon in char_matrix)
            self.assertNotIn(taxon, char_matrix)
            char_matrix[taxon] = seqs[idx]
        self.assertEqual(len(char_matrix._taxon_sequence_map), len(tns))
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        for idx, taxon in enumerate(tns):
            self.assertTrue(taxon in char_matrix)
            self.assertIn(taxon, char_matrix)
            del char_matrix[taxon]
            self.assertFalse(taxon in char_matrix)
            self.assertNotIn(taxon, char_matrix)
        self.assertEqual(len(char_matrix._taxon_sequence_map), 0)
        self.assertEqual(len(char_matrix), 0)

    def test_clear(self):
        tns = self.get_taxon_namespace(3)
        char_matrix = charmatrixmodel.CharacterMatrix(taxon_namespace=tns)
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        self.assertEqual(len(char_matrix), 0)
        seqs = [
                "abcd",
                [1,2,3,4,],
                ["a", "b", "c", "d",]
                ]
        assert len(seqs) == len(tns)
        for idx, taxon in enumerate(tns):
            self.assertFalse(taxon in char_matrix)
            self.assertNotIn(taxon, char_matrix)
            char_matrix[taxon] = seqs[idx]
        self.assertEqual(len(char_matrix._taxon_sequence_map), len(tns))
        self.assertEqual(len(char_matrix), len(char_matrix._taxon_sequence_map))
        char_matrix.clear()
        self.assertEqual(len(char_matrix._taxon_sequence_map), 0)
        self.assertEqual(len(char_matrix), 0)
        for idx, taxon in enumerate(tns):
            self.assertFalse(taxon in char_matrix)
            self.assertNotIn(taxon, char_matrix)

if __name__ == "__main__":
    unittest.main()
