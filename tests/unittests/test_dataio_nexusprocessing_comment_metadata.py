#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2026 Jeet Sukumaran, Mark T. Holder, and Matthew Andres Moreno.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Moreno, M. A., Holder, M. T., & Sukumaran, J. (2024). DendroPy 5: a
##     mature Python library for phylogenetic computing. Journal of Open
##     Source Software, 9(101), 6943, https://doi.org/10.21105/joss.06943
##
##############################################################################

"""
Tests for comment metadata extraction (``dendropy.dataio.nexusprocessing``),
including support for a callable ``extract_comment_metadata`` argument and
the BEAST2-version-specific parser added to address main fork issue #145
(https://github.com/jeetsukumaran/DendroPy/issues/145): DendroPy's original
comment metadata parser mis-parses nested list-valued ("vector")
annotations, e.g. ``history_all={{57,0.08,C,T},{134,0.079,A,G}}``.

Also covers ``parse_comment_metadata_beast2_v2_7_8_lark``, a
grammar-driven reimplementation of the same BEAST2 v2.7.8 comment
metadata parser built with Lark's Standalone Mode, via
``Beast2V2_7_8LarkCommentMetadataParsingTestCase``, which reruns the
hand-rolled parser's own test cases against it to confirm parity.
"""

import sys
import os
import unittest
import dendropy
from dendropy.dataio import nexusprocessing

sys.path.insert(0, os.path.dirname(__file__))
from support import dendropytest

# Verbatim (whitespace preserved) from a node comment in a real BEAST/
# TreeAnnotator-format MCC tree (``HA_discrete_MCC.tre``, bundled as an
# example with the FigTree v1.4.4 source distribution).
REAL_MCC_TREE_EXAMPLE_COMMENT = (
        '&rate_range={4.938776751387227E-4,0.036916549293719556},'
        'height_median=9.000000000000002,'
        'length=0.6825227857160625,'
        'state="D",'
        'state.prob=1.0,'
        'rate=0.007334968720519001'
        )

# The nested-list annotation from main fork issue #145.
ISSUE_145_COMMENT = "&history_all={{57,0.08,C,T},{134,0.079,A,G},{4,0.07,C,T}}"


class DendroPyV5_0_0CommentMetadataParsingTestCase(dendropytest.ExtendedTestCase):
    """
    ``parse_comment_metadata_dendropy_v5_0_0`` is DendroPy's original
    (pre-issue-145-fix) comment metadata parser, factored out unchanged
    so that it remains available (and is used as the default parser, for
    backward compatibility) despite its nested-list limitation.
    """

    def test_simple_scalar_values(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0(
                "&rate=0.5,label=foo")
        self.assertEqual(d, {"rate": "0.5", "label": "foo"})

    def test_quoted_and_boolean_values(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0(
                '&name="hello world",flag=true,off=FALSE')
        self.assertEqual(d, {"name": "hello world", "flag": True, "off": False})

    def test_single_level_vector(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0("&x={1,2,3}")
        self.assertEqual(d, {"x": ["1", "2", "3"]})

    def test_nhx_format(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0(
                "&&NHX:S=human:E=1.1.1.1")
        self.assertEqual(d, {"S": "human", "E": "1.1.1.1"})

    def test_nhx_format_with_explicit_prefix(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0(
                "&&NHX:S=human")
        d2 = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0("&&S=human")
        self.assertEqual(d, d2)

    def test_field_value_types_scalar(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0(
                "&age=5", field_value_types={"age": float})
        self.assertEqual(d, {"age": 5.0})

    def test_field_value_types_vector(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0(
                "&x={1,2,3}", field_value_types={"x": int})
        self.assertEqual(d, {"x": [1, 2, 3]})

    def test_unrecognized_comment_returns_empty(self):
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0("just a comment")
        self.assertEqual(d, {})

    def test_issue_145_nested_list_is_mis_parsed(self):
        # Regression-locks the originally-reported bug: the outer vector
        # value is truncated at the first (inner) closing brace.
        d = nexusprocessing.parse_comment_metadata_dendropy_v5_0_0(ISSUE_145_COMMENT)
        self.assertEqual(d["history_all"], ["{57", "0.08", "C", "T"])


class ParseCommentMetadataToAnnotationsBackwardCompatTestCase(dendropytest.ExtendedTestCase):
    """
    ``parse_comment_metadata_to_annotations`` is retained, with its
    original signature and (Annotation-set-returning) behavior, as a
    backward-compatibility wrapper around
    ``parse_comment_metadata_dendropy_v5_0_0`` +
    ``comment_metadata_to_annotations``.
    """

    def _as_dict(self, annotations):
        return {a.name: a.value for a in annotations}

    def test_basic(self):
        annotations = nexusprocessing.parse_comment_metadata_to_annotations("&rate=0.5")
        self.assertEqual(self._as_dict(annotations), {"rate": "0.5"})

    def test_field_name_map(self):
        annotations = nexusprocessing.parse_comment_metadata_to_annotations(
                "&rate=0.5", field_name_map={"rate": "substitution_rate"})
        self.assertEqual(self._as_dict(annotations), {"substitution_rate": "0.5"})

    def test_field_value_types_on_vector(self):
        annotations = nexusprocessing.parse_comment_metadata_to_annotations(
                "&x={1,2,3}", field_value_types={"x": int})
        self.assertEqual(self._as_dict(annotations), {"x": [1, 2, 3]})

    def test_accumulates_into_existing_set(self):
        existing = nexusprocessing.parse_comment_metadata_to_annotations("&a=1")
        combined = nexusprocessing.parse_comment_metadata_to_annotations(
                "&b=2", annotations=existing)
        self.assertIs(combined, existing)
        self.assertEqual(self._as_dict(combined), {"a": "1", "b": "2"})


class Beast2V2_7_8CommentMetadataParsingTestCase(dendropytest.ExtendedTestCase):
    """
    ``parse_comment_metadata_beast2_v2_7_8`` reproduces the comment
    metadata parsing behavior of BEAST2 v2.7.8's ``TreeParser``, based
    on a direct re-implementation of its ANTLR grammar
    (``NewickParser.g4``/``NewickLexer.g4``); unlike the other parsers
    in this module, it correctly (and non-destructively) handles
    arbitrarily-nested list-valued annotations.

    All test bodies are written against ``self.PARSE_FN`` rather than
    calling ``parse_comment_metadata_beast2_v2_7_8`` directly, so that
    ``Beast2V2_7_8LarkCommentMetadataParsingTestCase`` below can reuse
    them unchanged to confirm that
    ``parse_comment_metadata_beast2_v2_7_8_lark`` (a grammar-driven
    reimplementation built with Lark's Standalone Mode) behaves
    identically.
    """

    PARSE_FN = staticmethod(nexusprocessing.parse_comment_metadata_beast2_v2_7_8)

    def test_numbers_and_strings(self):
        d = self.PARSE_FN('&rate=0.0123,label="hello",tag=bareword')
        self.assertEqual(d, {"rate": 0.0123, "label": "hello", "tag": "bareword"})

    def test_single_quoted_value(self):
        d = self.PARSE_FN("&label='hello'")
        self.assertEqual(d, {"label": "hello"})

    def test_negative_and_scientific_notation_numbers(self):
        d = self.PARSE_FN("&a=-1.5,b=4.938776751387227E-4")
        self.assertEqual(d["a"], -1.5)
        self.assertAlmostEqual(d["b"], 4.938776751387227E-4)

    def test_all_numeric_vector(self):
        d = self.PARSE_FN("&hpd={1.1,2.2,3.3}")
        self.assertEqual(d, {"hpd": [1.1, 2.2, 3.3]})

    def test_mixed_vector_falls_back_to_raw_text(self):
        d = self.PARSE_FN('&x={1,"a"}')
        self.assertEqual(d["x"], ["1", '"a"'])

    def test_string_only_vector_falls_back_to_raw_text(self):
        d = self.PARSE_FN("&x={C,T,A,G}")
        self.assertEqual(d["x"], ["C", "T", "A", "G"])

    def test_nested_vector_resolves_issue_145(self):
        d = self.PARSE_FN(ISSUE_145_COMMENT)
        self.assertEqual(
                d["history_all"],
                ["{57,0.08,C,T}", "{134,0.079,A,G}", "{4,0.07,C,T}"])

    def test_doubly_nested_vector(self):
        d = self.PARSE_FN("&x={{{1,2},{3,4}},{5,6}}")
        self.assertEqual(d["x"], ["{{1,2},{3,4}}", "{5,6}"])

    def test_whitespace_is_tolerated(self):
        d = self.PARSE_FN("& rate = 0.5 , hpd = { 1.1 , 2.2 } ")
        self.assertEqual(d, {"rate": 0.5, "hpd": [1.1, 2.2]})

    def test_empty_comment_returns_empty(self):
        self.assertEqual(self.PARSE_FN("&"), {})

    def test_unrecognized_comment_returns_empty(self):
        d = self.PARSE_FN("not a metadata comment")
        self.assertEqual(d, {})

    def test_malformed_missing_value_raises(self):
        with self.assertRaises(ValueError):
            self.PARSE_FN("&rate=")

    def test_malformed_unbalanced_vector_raises(self):
        with self.assertRaises(ValueError):
            self.PARSE_FN("&x={1,2")

    def test_malformed_trailing_content_raises(self):
        with self.assertRaises(ValueError):
            self.PARSE_FN("&rate=0.5,,")

    def test_numeric_only_key_raises(self):
        # per the BEAST2 grammar, an attribute key must lex as ASTRING,
        # not as a number
        with self.assertRaises(ValueError):
            self.PARSE_FN("&123=5")

    def test_real_mcc_tree_example_comment(self):
        # exercised against realistic (non-nested) MCC-tree metadata
        d = self.PARSE_FN(REAL_MCC_TREE_EXAMPLE_COMMENT)
        self.assertEqual(d["state"], "D")
        self.assertEqual(d["state.prob"], 1.0)
        self.assertAlmostEqual(d["rate"], 0.007334968720519001)
        self.assertEqual(
                d["rate_range"],
                [4.938776751387227E-4, 0.036916549293719556])


class Beast2V2_7_8LarkCommentMetadataParsingTestCase(Beast2V2_7_8CommentMetadataParsingTestCase):
    """
    Reruns every ``Beast2V2_7_8CommentMetadataParsingTestCase`` test
    against ``parse_comment_metadata_beast2_v2_7_8_lark`` -- a
    grammar-driven reimplementation of the same BEAST2 v2.7.8 comment
    metadata format, generated (via Lark's Standalone Mode) from the
    grammar at ``dev/grammars/beast2_v2_7_8_comment_metadata.lark`` --
    to confirm it agrees with the hand-rolled recursive-descent parser
    on every case, error handling included.
    """

    PARSE_FN = staticmethod(nexusprocessing.parse_comment_metadata_beast2_v2_7_8_lark)


class CommentMetadataToAnnotationsTestCase(dendropytest.ExtendedTestCase):

    def _as_dict(self, annotations):
        return {a.name: a.value for a in annotations}

    def test_basic(self):
        annotations = nexusprocessing.comment_metadata_to_annotations(
                {"rate": 0.5, "label": "x"})
        self.assertEqual(self._as_dict(annotations), {"rate": 0.5, "label": "x"})

    def test_field_name_map(self):
        annotations = nexusprocessing.comment_metadata_to_annotations(
                {"rate": 0.5}, field_name_map={"rate": "substitution_rate"})
        self.assertEqual(self._as_dict(annotations), {"substitution_rate": 0.5})

    def test_field_value_types_scalar_and_list(self):
        annotations = nexusprocessing.comment_metadata_to_annotations(
                {"rate": 5, "hpd": [1, 2, 3]},
                field_value_types={"rate": float, "hpd": float})
        self.assertEqual(
                self._as_dict(annotations),
                {"rate": 5.0, "hpd": [1.0, 2.0, 3.0]})

    def test_empty_dict_yields_no_annotations(self):
        annotations = nexusprocessing.comment_metadata_to_annotations({})
        self.assertEqual(len(annotations), 0)


class ExtractCommentMetadataCallableIntegrationTestCase(dendropytest.ExtendedTestCase):
    """
    End-to-end tests confirming that a callable ``extract_comment_metadata``
    is honored by the tree readers (for both NEWICK and NEXUS schemas),
    and that using the BEAST2-style parser resolves main fork issue #145.
    """

    NEWICK_STR = (
            "((A[&rate=0.5,hpd={1.1,2.2}]:1,"
            "B[&history_all={{57,0.08,C,T},{134,0.079,A,G}}]:1):1,C:1);"
            )

    def _annotations_by_taxon_label(self, tree):
        result = {}
        for nd in tree:
            if nd.taxon is not None:
                result[nd.taxon.label] = nd.annotations.values_as_dict()
        return result

    def test_default_bool_true_matches_dendropy_v5_0_0(self):
        tree = dendropy.Tree.get(data=self.NEWICK_STR, schema="newick")
        result = self._annotations_by_taxon_label(tree)
        self.assertEqual(result["A"]["rate"], "0.5")
        self.assertEqual(result["B"]["history_all"], ["{57", "0.08", "C", "T"])

    def test_extract_comment_metadata_false(self):
        tree = dendropy.Tree.get(
                data=self.NEWICK_STR, schema="newick",
                extract_comment_metadata=False)
        for nd in tree:
            if nd.taxon is not None and nd.taxon.label == "B":
                self.assertEqual(len(nd.comments), 1)
                self.assertIn("history_all", nd.comments[0])
                self.assertEqual(len(nd.annotations), 0)

    def test_extract_comment_metadata_callable_beast2_resolves_issue_145(self):
        tree = dendropy.Tree.get(
                data=self.NEWICK_STR, schema="newick",
                extract_comment_metadata=nexusprocessing.parse_comment_metadata_beast2_v2_7_8)
        result = self._annotations_by_taxon_label(tree)
        self.assertEqual(
                result["B"]["history_all"],
                ["{57,0.08,C,T}", "{134,0.079,A,G}"])
        self.assertEqual(result["A"]["rate"], 0.5)
        self.assertEqual(result["A"]["hpd"], [1.1, 2.2])

    def test_extract_comment_metadata_callable_beast2_lark_matches(self):
        tree = dendropy.Tree.get(
                data=self.NEWICK_STR, schema="newick",
                extract_comment_metadata=nexusprocessing.parse_comment_metadata_beast2_v2_7_8_lark)
        result = self._annotations_by_taxon_label(tree)
        self.assertEqual(
                result["B"]["history_all"],
                ["{57,0.08,C,T}", "{134,0.079,A,G}"])
        self.assertEqual(result["A"]["rate"], 0.5)
        self.assertEqual(result["A"]["hpd"], [1.1, 2.2])

    def test_extract_comment_metadata_callable_via_nexus_schema(self):
        nexus_str = "#NEXUS\nBegin trees;\n  tree t1 = " + self.NEWICK_STR + "\nEnd;\n"
        tree = dendropy.Tree.get(
                data=nexus_str, schema="nexus",
                extract_comment_metadata=nexusprocessing.parse_comment_metadata_beast2_v2_7_8)
        result = self._annotations_by_taxon_label(tree)
        self.assertEqual(
                result["B"]["history_all"],
                ["{57,0.08,C,T}", "{134,0.079,A,G}"])

    def test_custom_user_callable(self):
        def my_parser(comment):
            return {"raw": comment}
        tree = dendropy.Tree.get(
                data=self.NEWICK_STR, schema="newick",
                extract_comment_metadata=my_parser)
        result = self._annotations_by_taxon_label(tree)
        self.assertEqual(result["A"]["raw"], "&rate=0.5,hpd={1.1,2.2}")

    def test_custom_callable_returning_empty_dict_falls_back_to_comments(self):
        # when metadata extraction yields nothing, the raw comment is
        # kept (as a plain comment) rather than silently dropped
        def no_op_parser(comment):
            return {}
        tree = dendropy.Tree.get(
                data=self.NEWICK_STR, schema="newick",
                extract_comment_metadata=no_op_parser)
        for nd in tree:
            if nd.taxon is not None and nd.taxon.label == "A":
                self.assertEqual(len(nd.annotations), 0)
                self.assertEqual(len(nd.comments), 1)


if __name__ == "__main__":
    unittest.main()
