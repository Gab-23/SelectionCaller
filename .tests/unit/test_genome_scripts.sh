#!/usr/bin/env bash
# Unit tests for the genome module (reference genome download + aligner indexing).
# No R scripts here - just checks that bwa/bowtie2 correctly index a toy genome.
# Usage: bash .tests/unit/test_genome_scripts.sh
# Requires: bwa, bowtie2

set -euo pipefail

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

FAIL=0
check() {
    local desc="$1"; local got="$2"; local expected="$3"
    if [ "$got" == "$expected" ]; then
        echo "  [OK] $desc"
    else
        echo "  [FAIL] $desc"
        echo "         expected: $expected"
        echo "         got:      $got"
        FAIL=1
    fi
}

echo "=== 1) bwa index (toy 2-chromosome genome) ==="
python3 -c "
import random
random.seed(1)
seq = ''.join(random.choice('ACGT') for _ in range(2000))
with open('$TMP/toy.fa', 'w') as f:
    f.write('>chr1\n')
    for i in range(0, len(seq), 70):
        f.write(seq[i:i+70] + '\n')
"
bwa index "$TMP/toy.fa" 2>"$TMP/bwa.log"
check "bwa produces the .bwt index file" "$([ -s "$TMP/toy.fa.bwt" ] && echo yes)" "yes"

echo "=== 2) bowtie2-build (same toy genome) ==="
bowtie2-build "$TMP/toy.fa" "$TMP/toy" >"$TMP/bt2.log" 2>&1
check "bowtie2-build produces the .1.bt2 index file" "$([ -s "$TMP/toy.1.bt2" ] && echo yes)" "yes"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi
