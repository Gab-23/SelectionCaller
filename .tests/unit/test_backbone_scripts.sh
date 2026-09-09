#!/usr/bin/env bash
# Unit test "leggeri" per gli script R del modulo backbone/annotazioni.
# Non richiedono l'intera DAG Snakemake: chiamano direttamente gli script con
# input fittizi (dati toy) e verificano l'output atteso, riga per riga.
#
# Uso:
#   bash .tests/unit/test_backbone_scripts.sh
#
# Richiede: Rscript (con dplyr, magrittr, R.utils, GenomicRanges, rtracklayer), bedtools

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

FAIL=0
check() {
    local desc="$1"; local got="$2"; local expected="$3"
    if [ "$got" == "$expected" ]; then
        echo "  [OK] $desc"
    else
        echo "  [FAIL] $desc"
        echo "         atteso:  $expected"
        echo "         ottenuto: $got"
        FAIL=1
    fi
}

echo "=== 1) get_backbone_file.R (modalita' binned) ==="
printf "chr1\t1000\nchr2\t800\n" > "$TMP/chrom.sizes"
Rscript workflow/scripts/get_backbone_file.R \
    --input_path_chrom_sizes="$TMP/chrom.sizes" \
    --input_path_mask=config/non_existing_file \
    --output_path_full="$TMP/backbone.bed" \
    --output_path_mask="$TMP/backbone.mask.bed" \
    --window_mode=binned --window_size=200 --step=0
check "5 bin su chr1 (1000/200)" "$(grep -c '^chr1' "$TMP/backbone.bed")" "5"
check "4 bin su chr2 (800/200)"  "$(grep -c '^chr2' "$TMP/backbone.bed")" "4"

echo "=== 2) get_actual_chromsizes_from_backbone.R (richiede dplyr) ==="
Rscript workflow/scripts/get_actual_chromsizes_from_backbone.R \
    --input_path="$TMP/backbone.bed" --output_path="$TMP/actual.chrom.sizes"
check "chrom size ricostruita da backbone" "$(cat "$TMP/actual.chrom.sizes")" "$(printf 'chr1\t1000\nchr2\t800')"

echo "=== 3) get_TSSs.R (conversione coordinate GFF3 1-based -> BED 0-based) ==="
cat > "$TMP/toy.gff3" << 'EOF'
##gff-version 3
chr1	HAVANA	transcript	100	500	.	+	.	ID=ENST001.1;gene_id=ENSG001.1;gene_name=GENEA;transcript_type=protein_coding;transcript_support_level=1;tag=basic
chr2	HAVANA	transcript	600	900	.	-	.	ID=ENST002.3;gene_id=ENSG002.3;gene_name=GENEB;transcript_type=protein_coding;transcript_support_level=2;tag=basic
EOF
gzip -c "$TMP/toy.gff3" > "$TMP/toy.gff3.gz"
Rscript workflow/scripts/get_TSSs.R --input_path="$TMP/toy.gff3.gz" --output_path="$TMP/tss.bed"
check "TSS strand + (BED 0-based: 99-100)" "$(grep ENST001 "$TMP/tss.bed" | cut -f1-3)" "$(printf 'chr1\t99\t100')"
check "TSS strand - (BED 0-based: 899-900)" "$(grep ENST002 "$TMP/tss.bed" | cut -f1-3)" "$(printf 'chr2\t899\t900')"

echo "=== 4) get_telomeres_proxy.R (deve derivare i cromosomi dall'input, non hard-codare chr1-22) ==="
printf "chr1\t1000\nchrX\t700\n" > "$TMP/chrom_with_X.sizes"
Rscript workflow/scripts/get_telomeres_proxy.R --input_path="$TMP/chrom_with_X.sizes" --output_path="$TMP/telomeres.bed"
check "presenza telomero iniziale chrX" "$(grep -c '^chrX\s0\s1$' "$TMP/telomeres.bed" || true)" "1"
check "nessuna riga fantasma per cromosomi non presenti (es. chr2)" "$(grep -c '^chr2' "$TMP/telomeres.bed" || true)" "0"

echo "=== 5) split_gtf.R (conversione coordinate GFF3 1-based -> BED 0-based su tutte le feature) ==="
cat > "$TMP/toy_split.gff3" << 'EOF'
##gff-version 3
chr1	HAVANA	gene	100	500	.	+	.	ID=ENSG001.1;gene_id=ENSG001.1;gene_type=protein_coding;gene_name=GENEA
EOF
gzip -c "$TMP/toy_split.gff3" > "$TMP/toy_split.gff3.gz"
Rscript workflow/scripts/split_gtf.R \
    --input_path="$TMP/toy_split.gff3.gz" \
    --output_path_protein_coding_5_UTR="$TMP/o1" --output_path_protein_coding_3_UTR="$TMP/o2" \
    --output_path_protein_coding_gene="$TMP/gene.bed" --output_path_protein_coding_exon="$TMP/o3" \
    --output_path_protein_coding_transcript="$TMP/o4" --output_path_lncRNA_gene="$TMP/o5" \
    --output_path_lncRNA_exon="$TMP/o6" --output_path_miRNA_gene="$TMP/o7" --output_path_miRNA_exon="$TMP/o8" \
    --output_path_rRNA_gene="$TMP/o9" --output_path_rRNA_exon="$TMP/o10"
check "gene GFF3 100-500 -> BED 99-500 (401bp)" "$(cut -f1-3 "$TMP/gene.bed")" "$(printf 'chr1\t99\t500')"

echo "=== 6) process_num_overlaps.R / process_dist_to_genes.R (conversione -1 -> NA) ==="
printf "chr1\t0\t200\t1_1\t.\t.\t.\t-1\n" > "$TMP/dist.bed"
Rscript workflow/scripts/process_dist_to_genes.R \
    --input_dist_protein_coding="$TMP/dist.bed" --input_dist_lncRNA="$TMP/dist.bed" \
    --input_dist_miRNA="$TMP/dist.bed" --input_dist_rRNA="$TMP/dist.bed" \
    --output_dist_protein_coding="$TMP/out_pc" --output_dist_lncRNA="$TMP/out_l" \
    --output_dist_miRNA="$TMP/out_m" --output_dist_rRNA="$TMP/out_r"
check "-1 convertito in NA" "$(tail -1 "$TMP/out_pc" | cut -f2)" "NA"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "TUTTI I TEST SONO PASSATI"
    exit 0
else
    echo "ALCUNI TEST SONO FALLITI"
    exit 1
fi
