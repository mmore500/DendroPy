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
Test cloning of dataobjects.
"""

import unittest
import copy

from dendropy.test.support import pathmap
from dendropy.test.support import datatest
import dendropy


class TestTreeCloning(datatest.AnnotatedDataObjectVerificationTestCase):

    ## Trees need special attention
    # def assertDistinctButEqualTree(self, tree1, tree2, **kwargs):
    #     otaxa = tree1.taxon_set
    #     ts = dendropy.TaxonSet()
    #     tree1.reindex_taxa(ts, clear=True)
    #     tree2.reindex_taxa(ts)
    #     self.assertIs(tree1.taxon_set, tree2.taxon_set)
    #     self.assertIsNot(tree1.taxon_set, otaxa)
    #     self.assertDistinctButEqual(tree1.taxon_set, otaxa, **kwargs)
    #     treesplit.encode_splits(tree1)
    #     treesplit.encode_splits(tree2)
    #     rfdist = treecalc.robinson_foulds_distance(tree1, tree2)
    #     self.assertAlmostEqual(rfdist, 0)

    def testDeepCopy(self):
        s = pathmap.tree_source_stream("pythonidae.annotated.nexml")
        d = dendropy.DataSet.get_from_stream(s, "nexml")
        tree1 = d.tree_lists[0][0]
        tree2 = copy.deepcopy(tree1)
        self.assertDistinctButEqualTree(tree1, tree2, distinct_taxa=False)

if __name__ == "__main__":
    unittest.main()
