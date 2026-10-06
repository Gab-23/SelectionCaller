#!/usr/bin/env bash
# Unit test for the assembly-aware GENCODE download (workflow/rules/gencode/get_gencode_annotation.smk):
# the URL built for each gencode_version, checked through a Snakemake dry-run of the real rule
# (no network). hg38 -> primary release; hg19 -> GRCh37_mapping/gencode.v{N}lift37.*
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
fail() { echo "FAIL: $1"; exit 1; }
pass() { echo "PASS: $1"; }

mkdir -p "$TMP/rules"
cp "$REPO/workflow/rules/gencode/get_gencode_annotation.smk" "$TMP/rules/"
printf 'include: "rules/get_gencode_annotation.smk"\n' > "$TMP/Snakefile"

urls() {   # print the wget URLs Snakemake would run for gencode_version $1
    (cd "$TMP" && snakemake "resources/reference_data/gencode/gencode.v$1.annotation.gff3.gz" \
        -s Snakefile -n -p --cores 1 2>&1) | grep -oE 'https://[^ ]+' | sort -u
}

BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"
got=$(urls 45)
exp=$(printf '%s\n' "$BASE/release_45/gencode.v45.annotation.gff3.gz" "$BASE/release_45/gencode.v45.annotation.gtf.gz" | sort)
[ "$got" == "$exp" ] || { echo "$got"; fail "hg38 URLs"; }
pass "hg38 (gencode_version 45): primary release URLs"

got=$(urls 45lift37)
exp=$(printf '%s\n' "$BASE/release_45/GRCh37_mapping/gencode.v45lift37.annotation.gff3.gz" "$BASE/release_45/GRCh37_mapping/gencode.v45lift37.annotation.gtf.gz" | sort)
[ "$got" == "$exp" ] || { echo "$got"; fail "hg19 URLs"; }
pass "hg19 (gencode_version 45lift37): GRCh37_mapping URLs"

got=$(urls 1)
exp=$(printf '%s\n' "$BASE/release_1/gencode.v1.annotation.gff3.gz" "$BASE/release_1/gencode.v1.annotation.gtf.gz" | sort)
[ "$got" == "$exp" ] || { echo "$got"; fail "single-digit version URLs"; }
pass "single-digit version (toy v1) unaffected"

echo "All GENCODE annotation rule tests passed."