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
##############################################################################

"""
Standalone (non-pytest/unittest) fuzz-testing harness for the comment
metadata parsers in ``dendropy.dataio.nexusprocessing``
(``parse_comment_metadata_dendropy_v5_0_0``,
``parse_comment_metadata_figtree_v1_4_4``, and
``parse_comment_metadata_beast2_v2_7_8``), added in support of resolving
main fork issue #145 (mis-parsing of nested list-valued comment
metadata).

This is deliberately kept outside of the pytest/unittest suite: it is not
a fixed-input regression check but a randomized, best-effort search for
crashes and misbehavior across the (very large) space of NEWICK/NEXUS
comment metadata strings, run on demand by developers rather than on
every CI invocation. It uses only the Python standard library -- no new
dependency (e.g. ``hypothesis``) is introduced.

Usage::

    python3 dev/fuzz_comment_metadata.py [--iterations N] [--seed N]
                                          [--max-depth N] [--verbose]

Exits with a non-zero status (and prints the offending input) on the
first unexpected failure found; otherwise prints a short summary.
"""

import argparse
import os
import random
import re
import sys
import traceback

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from dendropy.dataio import nexusprocessing

UNQUOTED_CHARS = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789|#*%/.-+_&:"
MUTATION_CHARS = UNQUOTED_CHARS + '{}=,"\' \t'
_NUMERIC_LOOKING_PATTERN = re.compile(r"-?[0-9.]+([eE]-?[0-9]+)?")


def random_unquoted_token(rng, min_len=1, max_len=8):
    n = rng.randint(min_len, max_len)
    return "".join(rng.choice(UNQUOTED_CHARS) for _ in range(n))


_KEY_FIRST_CHARS = UNQUOTED_CHARS.replace("&", "")


def random_key(rng):
    # BEAST2's grammar requires an attribute key to lex as a (non-numeric)
    # ASTRING token; avoid generating purely-numeric-looking keys here so
    # that "well-formed" fuzz cases are actually well-formed. (Rejection
    # of numeric-only keys is itself covered by the unit test suite.)
    #
    # Also avoid a leading '&': all three parsers treat a comment's
    # leading "&&" (vs. a single "&") as a distinct ("NHX-style") prefix
    # marker; a fuzzer-generated first key starting with '&' would, once
    # appended to the comment's own leading '&', silently collide with
    # that convention and change what is being parsed -- an inherent
    # ambiguity of the comment format itself (predating these parsers),
    # not something meaningful to fuzz for well-formed-input coverage.
    while True:
        key = rng.choice(_KEY_FIRST_CHARS) + random_unquoted_token(rng, 0, 9)
        if not _NUMERIC_LOOKING_PATTERN.fullmatch(key):
            return key


def random_number(rng):
    style = rng.randrange(3)
    if style == 0:
        return str(rng.randint(-100000, 100000))
    elif style == 1:
        return "{}.{}".format(rng.randint(0, 10000), rng.randint(0, 999999))
    else:
        mantissa = round(rng.uniform(-9.999999, 9.999999), 6)
        exponent = rng.randint(0, 30)
        sign = rng.choice(("", "-"))
        return "{}E{}{}".format(mantissa, sign, exponent)


def random_quoted_string(rng):
    quote = rng.choice(('"', "'"))
    body = random_unquoted_token(rng, 0, 10).replace(quote, "")
    return quote + body + quote


def random_value(rng, depth, max_depth):
    if depth < max_depth and rng.random() < 0.25:
        n = rng.randint(1, 4)
        elements = [random_value(rng, depth + 1, max_depth) for _ in range(n)]
        return "{" + ",".join(elements) + "}"
    choice = rng.random()
    if choice < 0.5:
        return random_number(rng)
    elif choice < 0.8:
        return random_quoted_string(rng)
    else:
        return random_unquoted_token(rng, 1, 8)


def random_comment(rng, min_attribs=1, max_attribs=6, max_depth=3):
    n = rng.randint(min_attribs, max_attribs)
    parts = ["{}={}".format(random_key(rng), random_value(rng, 0, max_depth))
             for _ in range(n)]
    return "&" + ",".join(parts)


def mutate(rng, text, n_mutations=None):
    chars = list(text)
    if n_mutations is None:
        n_mutations = rng.randint(1, 4)
    for _ in range(n_mutations):
        if not chars:
            op = "insert"
        else:
            op = rng.choice(("insert", "delete", "replace"))
        if op == "insert":
            pos = rng.randint(0, len(chars))
            chars.insert(pos, rng.choice(MUTATION_CHARS))
        elif op == "delete":
            pos = rng.randrange(len(chars))
            del chars[pos]
        else:
            pos = rng.randrange(len(chars))
            chars[pos] = rng.choice(MUTATION_CHARS)
    return "".join(chars)


class FuzzFailure(AssertionError):
    pass


def check_well_formed_case(comment):
    """
    ``comment`` is generated to be grammatically well-formed (per the
    BEAST2 v2.7.8 grammar that all three parsers are, at minimum, a
    superset-tolerant of). None of the three parsers should raise for
    such input, and results should be deterministic.
    """
    try:
        beast2_result = nexusprocessing.parse_comment_metadata_beast2_v2_7_8(comment)
    except nexusprocessing.Beast2CommentMetadataError:
        raise FuzzFailure(
                "beast2_v2_7_8 parser rejected well-formed input: {!r}".format(comment))
    beast2_result_again = nexusprocessing.parse_comment_metadata_beast2_v2_7_8(comment)
    if beast2_result != beast2_result_again:
        raise FuzzFailure(
                "beast2_v2_7_8 parser is non-deterministic for: {!r}\n"
                "  first:  {!r}\n  second: {!r}".format(
                    comment, beast2_result, beast2_result_again))
    for parse_fn in (nexusprocessing.parse_comment_metadata_dendropy_v5_0_0,
            nexusprocessing.parse_comment_metadata_figtree_v1_4_4):
        try:
            parse_fn(comment)
        except ValueError:
            # dendropy_v5_0_0 and figtree_v1_4_4 are regex-driven and
            # generally permissive, but are documented to be able to
            # raise ValueError (e.g. figtree_v1_4_4 on a badly-formatted
            # attribute); that is acceptable, just not any *other*
            # exception type.
            pass


def check_mutated_case(comment):
    """
    ``comment`` is a randomly-mutated (and so, generally, *not*
    well-formed) string. Each parser is allowed to raise only the
    exception type(s) it documents for malformed input; anything else
    (``IndexError``, ``AttributeError``, ``RecursionError``, an infinite
    loop, etc.) is a bug.
    """
    for parse_fn, allowed_exceptions in (
            (nexusprocessing.parse_comment_metadata_dendropy_v5_0_0, (ValueError,)),
            (nexusprocessing.parse_comment_metadata_figtree_v1_4_4, (ValueError,)),
            (nexusprocessing.parse_comment_metadata_beast2_v2_7_8,
                (nexusprocessing.Beast2CommentMetadataError,)),
            ):
        try:
            parse_fn(comment)
        except allowed_exceptions:
            pass
        except Exception as e:
            raise FuzzFailure(
                    "{} raised unexpected {} on mutated input {!r}: {}".format(
                        parse_fn.__name__, type(e).__name__, comment, e))


def check_end_to_end_tree_parsing(rng, max_depth):
    """
    Embeds fuzzed comments into a small randomly-shaped NEWICK tree and
    round-trips it through ``dendropy.Tree.get(..., extract_comment_metadata=...)``,
    exercising the callable-``extract_comment_metadata`` integration path
    (not just the parser functions in isolation).
    """
    import dendropy

    n_leaves = rng.randint(2, 6)
    labels = ["T{}".format(i) for i in range(n_leaves)]

    def build(labels):
        if len(labels) == 1:
            comment = random_comment(rng, max_depth=max_depth)
            return "{}[{}]".format(labels[0], comment)
        split = rng.randint(1, len(labels) - 1)
        left = build(labels[:split])
        right = build(labels[split:])
        comment = random_comment(rng, max_depth=max_depth)
        return "({},{})[{}]".format(left, right, comment)

    newick_str = build(labels) + ";"
    for extract_comment_metadata, allowed_exceptions in (
            (True, ()),
            (False, ()),
            (nexusprocessing.parse_comment_metadata_beast2_v2_7_8,
                (nexusprocessing.Beast2CommentMetadataError,)),
            # figtree_v1_4_4 is regex-driven and, per its JEBL origin,
            # legitimately raises ValueError on certain malformed inputs
            # that only arise as an artifact of its own nested-bracket
            # mishandling (e.g. an empty quoted-string attribute label
            # left over from a mis-parsed nested vector); this is
            # documented, faithful behavior, not a bug -- see the
            # figtree_v1_4_4 unit tests and module-level comments in
            # nexusprocessing.py.
            (nexusprocessing.parse_comment_metadata_figtree_v1_4_4, (ValueError,)),
            ):
        try:
            tree = dendropy.Tree.get(
                    data=newick_str,
                    schema="newick",
                    extract_comment_metadata=extract_comment_metadata)
        except allowed_exceptions:
            continue
        except Exception as e:
            raise FuzzFailure(
                    "dendropy.Tree.get(extract_comment_metadata={!r}) raised "
                    "unexpected {} on {!r}: {}".format(
                        extract_comment_metadata, type(e).__name__, newick_str, e))
        if len(tree.leaf_nodes()) != n_leaves:
            raise FuzzFailure(
                    "tree parsed from {!r} has {} leaves, expected {}".format(
                        newick_str, len(tree.leaf_nodes()), n_leaves))


def probe_deep_nesting(max_depth_to_try, step):
    """
    Not a pass/fail check: reports the nesting depth (if any, up to
    ``max_depth_to_try``) at which the recursive-descent BEAST2 parser
    hits Python's recursion limit. This is informational (deeply nested
    metadata is not expected in practice), surfaced so a human can decide
    whether an explicit depth guard is warranted.
    """
    depth = step
    while depth <= max_depth_to_try:
        comment = "&x=" + ("{" * depth) + "1" + ("}" * depth)
        try:
            nexusprocessing.parse_comment_metadata_beast2_v2_7_8(comment)
        except RecursionError:
            return depth
        except nexusprocessing.Beast2CommentMetadataError:
            pass
        depth += step
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iterations", type=int, default=2000,
            help="number of fuzz iterations per check category (default: 2000)")
    parser.add_argument("--seed", type=int, default=None,
            help="PRNG seed, for reproducing a specific run (default: random)")
    parser.add_argument("--max-depth", type=int, default=4,
            help="maximum nesting depth for generated vector values (default: 4)")
    parser.add_argument("--skip-deep-nesting-probe", action="store_true",
            help="skip the (slower) recursion-limit probe")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    seed = args.seed if args.seed is not None else random.SystemRandom().randrange(2**32)
    rng = random.Random(seed)
    print("dev/fuzz_comment_metadata.py: seed={}, iterations={}, max_depth={}".format(
            seed, args.iterations, args.max_depth))

    n_checked = 0
    try:
        for i in range(args.iterations):
            comment = random_comment(rng, max_depth=args.max_depth)
            if args.verbose:
                print("[well-formed {}] {}".format(i, comment))
            check_well_formed_case(comment)
            n_checked += 1

        for i in range(args.iterations):
            base = random_comment(rng, max_depth=args.max_depth)
            mutated = mutate(rng, base)
            if args.verbose:
                print("[mutated {}] {}".format(i, mutated))
            check_mutated_case(mutated)
            n_checked += 1

        for i in range(args.iterations // 10 or 1):
            if args.verbose:
                print("[end-to-end tree {}]".format(i))
            check_end_to_end_tree_parsing(rng, max_depth=min(args.max_depth, 3))
            n_checked += 1
    except FuzzFailure as e:
        print("\nFUZZ FAILURE after {} checks (seed={}):\n{}".format(
                n_checked, seed, e))
        traceback.print_exc()
        sys.exit(1)

    if not args.skip_deep_nesting_probe:
        recursion_limit_depth = probe_deep_nesting(max_depth_to_try=3000, step=100)
        if recursion_limit_depth is not None:
            print("NOTE: beast2_v2_7_8 parser hits Python's recursion limit "
                    "(RecursionError) at ~{} levels of vector nesting.".format(
                        recursion_limit_depth))
        else:
            print("NOTE: no RecursionError observed for vector nesting up to 3000 levels.")

    print("OK: {} checks passed (seed={})".format(n_checked, seed))


if __name__ == "__main__":
    main()
