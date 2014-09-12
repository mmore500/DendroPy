#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
DEPRECATED IN DENDROPY 4: USE `dendropy.simulate.charsim`.
"""

import dendropy
from dendropy.simulate import charsim
from dendropy.utility import deprecate

def generate_hky_dataset(seq_len,
                         tree_model,
                         mutation_rate=1.0,
                         kappa=1.0,
                         base_freqs=[0.25, 0.25, 0.25, 0.25],
                         root_states=None,
                         dataset=None,
                         rng=None):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.seqsim.generate_hky_dataset()' function has been replaced with 'dendropy.simulate.charsim.hky85_char_matrix()'.",
            old_construct="from dendropy import seqsim\ndataset = seqsim.generate_hky_dataset(...)",
            new_construct="from dendropy.simulate import charsim\nchar_matrix = discrete.hky85_char_matrix(...)")
    if dataset is None:
        dataset = dendropy.DataSet()
    char_matrix = dataset.new_char_matrix(char_matrix_type="dna", taxon_namespace=tree_model.taxon_namespace)
    return charsim.hky85_char_matrix(
            seq_len=seq_len,
            tree_model=tree_model,
            mutation_rate=mutation_rate,
            kappa=kappa,
            base_freqs=base_freqs,
            root_states=root_states,
            char_matrix=char_matrix,
            rng=rng)

def generate_hky_characters(seq_len,
                            tree_model,
                            mutation_rate=1.0,
                            kappa=1.0,
                            base_freqs=[0.25, 0.25, 0.25, 0.25],
                            root_states=None,
                            char_matrix=None,
                            rng=None):
    pass

def generate_dataset(seq_len,
                     tree_model,
                     seq_model,
                     mutation_rate=1.0,
                     root_states=None,
                     dataset=None,
                     rng=None):
    pass

def generate_char_matrix(seq_len,
                        tree_model,
                        seq_model,
                        mutation_rate=1.0,
                        root_states=None,
                        char_matrix=None,
                        rng=None):

    pass

class SeqEvolver(object):
    def __init__(self,
     seq_model=None,
     mutation_rate=None,
     seq_attr='sequences',
     seq_model_attr="seq_model",
     edge_length_attr="length",
     edge_rate_attr="mutation_rate",
     seq_label_attr='taxon'):
        pass

