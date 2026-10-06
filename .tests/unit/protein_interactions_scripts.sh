#!/usr/bin/env bash
# Unit tests for the protein-interactions module (toy data, hand-computed expectations).
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
TMP=$(mktemp -d)
trap '[ -n "${KEEP_TMP:-}" ] || rm -rf "$TMP"' EXIT
fail() { echo "FAIL: $1"; exit 1; }
pass() { echo "PASS: $1"; }

# Toy CORUM: complexes 1 and 2 share the name "X" (must NOT be merged: grouping is by ID);
# complex 3 is a single subunit with whitespace; complex 4 has an empty subunit.
# Toy Interactome INSIDER: a homodimer (A-A) and the same pair in both orders (C-G, G-C).
Rscript -e '
human_proteincomplex <- data.frame(ComplexID = c(1, 2, 3, 4), ComplexName = c("X", "X", "Y", "Z"),
                                   subunits.Gene.name. = c("A;B;C", "A;D", " E ", "B;;F"),
                                   stringsAsFactors = FALSE)
save(human_proteincomplex, file = file.path(commandArgs(TRUE)[1], "corum.RData"))
HS_interfaces_genenames <- data.frame(P1 = factor(c("A", "A", "C", "G", "D")), P2 = factor(c("B", "A", "G", "C", "H")),  # factors, as in the real file
                                      P1_IRES = "[]", P2_IRES = "[]", stringsAsFactors = FALSE)
save(HS_interfaces_genenames, file = file.path(commandArgs(TRUE)[1], "insider.RData"))' "$TMP"

Rscript "$REPO/workflow/scripts/build_gene_interactors.R" --input_path_CORUM="$TMP/corum.RData" \
    --input_path_InteractomeINSIDER="$TMP/insider.RData" --output_path="$TMP/interactors.tsv" 2> "$TMP/err.log" \
    || { cat "$TMP/err.log"; fail "build_gene_interactors.R on toy data"; }

Rscript - "$TMP" <<'RS' || fail "build_gene_interactors.R output checks"
d <- read.delim(file.path(commandArgs(TRUE)[1], "interactors.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(colnames(d), c("gene", "partner", "source")))
key <- function(x) sort(paste(x$gene, x$partner))
corum <- d[d$source == "CORUM", ]
exp_corum <- c("A B","A C","A D","B A","B C","B F","C A","C B","D A","E NA","F B")
stopifnot(identical(key(corum), sort(exp_corum)))          # no B-D / C-D: complexes 1, 2 not merged by name
ppi <- d[d$source == "InteractomeINSIDER", ]
stopifnot(identical(key(ppi), sort(c("A B","B A","C G","G C","D H","H D"))))  # no homodimer, no duplicates
RS
pass "build_gene_interactors.R: grouping by ComplexID, lonely subunit kept, homodimers/duplicates removed"

# Toy annotated backbone (chr, start, end, binID, comma-separated genes)
# Two extra bins carry a gene symbol with an apostrophe (real case: Interactome INSIDER
# "PSG11s'"): with read.table's default quote handling, the rows between them would merge.
printf "chr1\t0\t100\tchr1_0_100\tA,B\nchr1\t100\t200\tchr1_100_200\tC\nchr1\t200\t300\tchr1_200_300\tQ'\nchr2\t0\t100\tchr2_0_100\tE,G,X1\nchr2\t100\t200\tchr2_100_200\tH\nchr2\t200\t300\tchr2_200_300\tR'#x\n" \
    > "$TMP/backbone.bed"
Rscript "$REPO/workflow/scripts/compute_interactor_features.R" --input_path_backbone="$TMP/backbone.bed" \
    --input_path_interactors="$TMP/interactors.tsv" --output_path_features="$TMP/features.tsv" \
    --output_path_gene_sets="$TMP/gene_sets.tsv" 2> "$TMP/err2.log" \
    || { cat "$TMP/err2.log"; fail "compute_interactor_features.R on toy data"; }

Rscript - "$TMP" <<'RS' || fail "compute_interactor_features.R output checks"
P <- commandArgs(TRUE)[1]
f <- read.delim(file.path(P, "features.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(colnames(f), c("binID","complex_density","n_complex_partners_trans","n_PPIs_trans")))
rownames(f) <- f$binID
# Hand-computed:
#  chr1_0_100   {A,B}: CORUM partners {B,C,D}u{A,C,F} minus {A,B} = {C,D,F}; PPIs {B}u{A} minus bin = {}
#  chr1_100_200 {C}:   CORUM {A,B}; PPIs {G}
#  chr2_0_100   {E,G,X1}: E is a lonely subunit (density 1/3), no partners; PPIs of G = {C}
#  chr2_100_200 {H}:   not a subunit; PPIs {D}
exp <- data.frame(complex_density = c(1, 1, 1/3, 0), n_complex_partners_trans = c(3, 2, 0, 0),
                  n_PPIs_trans = c(0, 1, 1, 1),
                  row.names = c("chr1_0_100","chr1_100_200","chr2_0_100","chr2_100_200"))
stopifnot(isTRUE(all.equal(f[rownames(exp), -1], exp, check.attributes = FALSE)))
stopifnot(nrow(f) == 6, all(f[c("chr1_200_300","chr2_200_300"), "complex_density"] == 0))   # apostrophe/hash bins kept
g <- read.delim(file.path(P, "gene_sets.tsv"), stringsAsFactors = FALSE)
s <- function(b, set) sort(g$gene[g$binID == b & g$gene_set == set])
stopifnot(identical(s("chr1_0_100", "partners_trans"), c("C","D","F")),
          identical(s("chr1_100_200", "all_int_trans"), c("A","B","G")),
          identical(s("chr2_0_100", "genes"), c("E","G","X1")),
          length(s("chr1_0_100", "ppis_trans")) == 0,
          !anyNA(g$gene),
          identical(s("chr2_200_300", "genes"), "R'#x"))
RS
pass "compute_interactor_features.R: hand-computed density, trans partners/PPIs and gene sets"

echo "All protein-interactions tests passed."