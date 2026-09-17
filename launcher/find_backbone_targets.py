#!/usr/bin/env python3
"""
Scan Snakemake rule files and print, for every rule whose OUTPUT depends on the
backbone/binning wildcards (window_mode / window_size / step), the glob patterns
that match its produced files, plus the log directory to remove.

A rule is considered "backbone-dependent" iff at least one of its declared
`output:` paths contains {window_mode}, {window_size} or {step}. This is a
structural fact about the pipeline (these wildcards only exist once the
backbone/bin grid has been fixed), not a guess based on directory naming.

Usage: find_backbone_targets.py <rules_dir_or_flat_dir> [--flat]
Prints two blocks to stdout, machine-parseable by the caller (bash):
    OUTPUT_GLOB<TAB><pattern>
    LOG_DIR<TAB><rule_name>
"""
import re
import sys
import glob
import os

BACKBONE_WILDCARDS = ("{window_mode}", "{window_size}", "{step}")

# Keywords that can close an `output:` block at the same rule-body indentation
# level. Indentation in this codebase is inconsistent (tabs in some files,
# 4 spaces in others), so we match on the keyword only, anchored to the start
# of a line after arbitrary leading whitespace - not on a fixed indent width.
BLOCK_ENDERS = r"(?:input|params|log|conda|shell|threads|resources|wildcard_constraints|run|priority|benchmark)\s*:"


def extract_block(text, start_keyword_pattern):
    """Return the text between a `<start_keyword>:` line and the next
    BLOCK_ENDERS keyword, tolerant of tabs or spaces as indentation."""
    m = re.search(r"(?:^|\n)[ \t]*" + start_keyword_pattern + r"\s*:", text)
    if not m:
        return None
    rest = text[m.end():]
    end = re.search(r"\n[ \t]*" + BLOCK_ENDERS, rest)
    return rest[: end.start()] if end else rest


def parse_rule_file(path):
    """Return (rule_name, output_globs) or None if the file has no `rule`."""
    text = open(path).read()
    rmatch = re.search(r"(?:^|\n)rule\s+(\w+)\s*:", text)
    if not rmatch:
        return None
    rule_name = rmatch.group(1)

    out_block = extract_block(text, "output")
    if out_block is None:
        return None

    if not any(w in out_block for w in BACKBONE_WILDCARDS):
        return None  # not backbone-dependent

    raw_paths = re.findall(r'"([^"]+)"', out_block)
    globs = sorted(
        {re.sub(r"\{[^}]+\}", "*", p) for p in raw_paths if "{" in p}
    )
    return rule_name, globs


def main():
    if len(sys.argv) < 2:
        print("usage: find_backbone_targets.py <dir> [--flat]", file=sys.stderr)
        sys.exit(2)
    root = sys.argv[1]
    flat = "--flat" in sys.argv[2:]

    if flat:
        smk_files = sorted(glob.glob(os.path.join(root, "*.smk")))
    else:
        smk_files = sorted(glob.glob(os.path.join(root, "**", "*.smk"), recursive=True))

    seen_rules = set()
    for f in smk_files:
        result = parse_rule_file(f)
        if result is None:
            continue
        rule_name, globs = result
        if rule_name in seen_rules:
            continue  # defensive: avoid double-processing if a rule appears twice
        seen_rules.add(rule_name)
        for g in globs:
            print(f"OUTPUT_GLOB\t{g}")
        print(f"LOG_DIR\t{rule_name}")


if __name__ == "__main__":
    main()