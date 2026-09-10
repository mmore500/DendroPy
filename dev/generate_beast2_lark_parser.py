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
Regenerates the checked-in standalone parser used by
``parse_comment_metadata_beast2_v2_7_8_lark()``
(``dendropy.dataio.nexusprocessing``) from its source grammar at
``dev/grammars/beast2_v2_7_8_comment_metadata.lark``.

This is a dev-time-only tool. It requires the third-party ``lark``
package (``pip install lark``); DendroPy itself has no runtime
dependency on ``lark`` -- the generated output at
``src/dendropy/dataio/_beast2_v2_7_8_lark_standalone.py`` is a
self-contained parser module with no import of ``lark``, produced via
Lark's "Standalone Mode" (``python -m lark.tools.standalone``).

Note that the *generated* module (unlike the rest of DendroPy, and
unlike ``lark`` itself, both BSD/MIT-licensed) is subject to a separate
Mozilla Public License, v. 2.0, per the notice embedded in that tool's
own source header; see the top of the generated file, and item 5 of
"NOTICES.rst", for details.

The ``basic`` (non-contextual) lexer is required, rather than the
LALR default of ``contextual``: it is what makes a purely-numeric token
lex as NUMBER (rather than ASTRING) regardless of grammar position,
matching the context-free lexing of BEAST2's own ANTLR grammar and so
correctly rejecting numeric-only attribute keys (e.g. "&123=5").

Usage::

    python3 dev/generate_beast2_lark_parser.py
"""

import os
import subprocess
import sys

DEV_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(DEV_DIR)
GRAMMAR_PATH = os.path.join(DEV_DIR, "grammars", "beast2_v2_7_8_comment_metadata.lark")
OUTPUT_PATH = os.path.join(
        REPO_ROOT, "src", "dendropy", "dataio", "_beast2_v2_7_8_lark_standalone.py")

def main():
    try:
        import lark  # noqa: F401
    except ImportError:
        sys.exit(
                "This tool requires the third-party 'lark' package"
                " (dev-time only; not a DendroPy runtime dependency).\n"
                "Install it with: pip install lark")
    result = subprocess.run(
            [sys.executable, "-m", "lark.tools.standalone",
                "-l", "basic", GRAMMAR_PATH],
            stdout=subprocess.PIPE,
            check=True)
    with open(OUTPUT_PATH, "wb") as dest:
        dest.write(result.stdout)
    print("Wrote: {}".format(OUTPUT_PATH))

if __name__ == "__main__":
    main()
