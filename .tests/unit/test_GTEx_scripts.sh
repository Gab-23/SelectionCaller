#!/usr/bin/env bash
# Unit tests for the GTEx module (toy data, local HTTP server, no internet).
# Download tests run a copy of the real rule file workflow/rules/GTEx/get_GTEx_median_TPM.smk
# with only GTEx_BASE_URL redirected to a local python http.server; test files are served at
# the exact relative paths of the real release table (v8/rna-seq/..., v11/rna-seq/...).
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
TMP=$(mktemp -d)
SERVER_PID=""
trap 'kill $SERVER_PID 2>/dev/null || true; [ -n "${KEEP_TMP:-}" ] || rm -rf "$TMP"' EXIT

fail() { echo "FAIL: $1"; exit 1; }
pass() { echo "PASS: $1"; }

# --- toy payloads ---
mkdir -p "$TMP/payloads"
printf '#1.2\n2\t2\nName\tDescription\tBreast - Mammary Tissue\tLung\nENSG00000000003.14\tTSPAN6\t10.5\t3.2\nENSG00000000005.5\tTNMD\t0.1\t0\n' \
    | gzip -c > "$TMP/payloads/good.gct.gz"
printf '#1.2\r\n2\t1\r\nName\tDescription\tLung\r\nENSG1\tA\t1\r\nENSG2\tB\t2\r\n' \
    | gzip -c > "$TMP/payloads/crlf.gct.gz"
printf '<html>not a gzip</html>\n' > "$TMP/payloads/notgzip.gct.gz"
printf 'chr1\t1\t2\n' | gzip -c > "$TMP/payloads/notgct.gct.gz"

# --- minimal workflow: copy of the real rule, base URL redirected to the local server ---
PORT=$(python3 -c 'import socket; s=socket.socket(); s.bind(("127.0.0.1",0)); print(s.getsockname()[1]); s.close()')
mkdir -p "$TMP/wf/rules/GTEx" "$TMP/www"
RULE="$TMP/wf/rules/GTEx/get_GTEx_median_TPM.smk"
cp "$REPO/workflow/rules/GTEx/get_GTEx_median_TPM.smk" "$RULE"
grep -q '^GTEx_BASE_URL = "https://' "$RULE" || fail "GTEx_BASE_URL not found in the rule file"
sed -i "s|^GTEx_BASE_URL = .*|GTEx_BASE_URL = \"http://127.0.0.1:$PORT\"|" "$RULE"
printf 'include: "rules/GTEx/get_GTEx_median_TPM.smk"\n' > "$TMP/wf/Snakefile"

# Relative path of each release's file, read from the real rule (so the test follows the table)
relpath() { python3 - "$REPO/workflow/rules/GTEx/get_GTEx_median_TPM.smk" "$1" <<'PY'
import re, sys, ast
src = open(sys.argv[1]).read()
table = ast.literal_eval(re.search(r"GTEx_median_TPM_files = (\{.*?\})", src, re.S).group(1))
print(table[sys.argv[2]])
PY
}
V8=$(relpath v8); V11=$(relpath v11)
serve() { mkdir -p "$TMP/www/$(dirname "$2")"; cp "$TMP/payloads/$1" "$TMP/www/$2"; }   # serve <payload> at <relpath>

(cd "$TMP/www" && exec python3 -m http.server "$PORT" --bind 127.0.0.1 > /dev/null 2>&1) &
SERVER_PID=$!
for _ in $(seq 1 50); do curl -s -o /dev/null "http://127.0.0.1:$PORT/" && break; sleep 0.1; done

run() { (cd "$TMP/wf" && snakemake "$1" -s Snakefile --cores 1 > "$TMP/run.out" 2>&1); }
out() { echo "$TMP/wf/resources/GTEx/$1/gene_median_tpm.gct.gz"; }

# 1. valid file (v8): downloaded byte-identical from the v8 path, no .part, URL logged
serve good.gct.gz "$V8"
run resources/GTEx/v8/gene_median_tpm.gct.gz || { cat "$TMP/run.out"; fail "valid GCT download"; }
cmp -s "$(out v8)" "$TMP/payloads/good.gct.gz" || fail "downloaded file differs from source"
[ ! -e "$(out v8).part" ] || fail ".part file left behind after success"
grep -q "$V8" "$TMP/wf/logs/get_GTEx_median_TPM/v8.log" || fail "v8 URL not written to log"
pass "valid GCT is downloaded unchanged from the v8 location"

# 2. CRLF line endings in the header are accepted (v11 location)
serve crlf.gct.gz "$V11"
run resources/GTEx/v11/gene_median_tpm.gct.gz || { cat "$TMP/run.out"; fail "CRLF GCT rejected"; }
grep -q "$V11" "$TMP/wf/logs/get_GTEx_median_TPM/v11.log" || fail "v11 URL not written to log"
pass "CRLF GCT header accepted (v11 location)"

# 3-5. invalid payloads behind the v8 URL: rule fails and leaves neither output nor .part
for case in "notgzip.gct.gz:non-gzip payload" "notgct.gct.gz:gzip but not GCT" "MISSING:HTTP 404"; do
    payload=${case%%:*}; label=${case#*:}
    rm -f "$(out v8)" "$TMP/www/$V8"
    [ "$payload" == "MISSING" ] || serve "$payload" "$V8"
    if run resources/GTEx/v8/gene_median_tpm.gct.gz; then fail "$label accepted"; fi
    [ ! -e "$(out v8)" ] || fail "$label left an output file"
    [ ! -e "$(out v8).part" ] || fail "$label left a .part file"
    pass "$label rejected, nothing left behind"
done

# 6. release with no download location: explicit error naming it
if run resources/GTEx/v99/gene_median_tpm.gct.gz; then fail "unknown release v99 accepted"; fi
grep -q "GTEx release 'v99' has no download location" "$TMP/run.out" || { cat "$TMP/run.out"; fail "unknown-release error message"; }
pass "unknown release raises an explicit error"

# 7. release not matching the wildcard constraint is not resolvable
if run resources/GTEx/latest/gene_median_tpm.gct.gz; then fail "non-'v<N>' release accepted"; fi
pass "release wildcard constrained to v<N>"

# 8. no links in the configs (they hold local paths and parameters only)
for c in "$REPO/config/config_dataset.yaml" "$REPO/.tests/integration/config/config_dataset.yaml"; do
    ! grep -nE "https?://" "$c" || fail "link found in $c"
done
pass "no URLs in config files"


# =========================================================================================
# process_GTEx.R
# =========================================================================================
P="$TMP/process"; mkdir -p "$P"
# Toy GCT: raw GTEx tissue names (with " - " and spaces), one unmapped tissue (Heart),
# a gene not expressed anywhere (G2), a PAR_Y duplicate after version stripping (G3).
printf '#1.2\n5\t7\n' > "$P/toy.gct"
printf 'Name\tDescription\tBreast - Mammary Tissue\tLung\tKidney - Cortex\tBrain - Cortex\tBrain - Cerebellum\tBrain - Hippocampus\tHeart - Left Ventricle\n' >> "$P/toy.gct"
printf 'ENSG01.3\tG1\t10\t2\t4\t1\t7\t3\t100\n'      >> "$P/toy.gct"
printf 'ENSG02.1\tG2\t0\t0\t0\t0\t0\t0\t5\n'         >> "$P/toy.gct"
printf 'ENSG03.5\tG3\t1\t1\t1\t1\t1\t1\t1\n'         >> "$P/toy.gct"
printf 'ENSG03.5_PAR_Y\tG3\t5\t0\t3\t9\t0\t0\t0\n'   >> "$P/toy.gct"
printf 'ENSG04.2\tG4\t3\t8\t1\t2\t6\t20\t0\n'        >> "$P/toy.gct"
gzip -c "$P/toy.gct" > "$P/toy.gct.gz"
printf 'GTEx_tissue\tTCGA_cohort\nBreast...Mammary.Tissue\tBRCA\nLung\tLUAD\nLung\tLUSC\nKidney...Cortex\tKIRC\nBrain...Cortex\tGBMLGG\nBrain...Cerebellum\tGBMLGG\nBrain...Hippocampus\tGBMLGG\n' > "$P/mapping.tsv"

Rscript "$REPO/workflow/scripts/process_GTEx.R" --input_path_gct="$P/toy.gct.gz" \
    --input_path_mapping="$P/mapping.tsv" --output_path="$P/out.tsv" 2> "$P/err.log" \
    || { cat "$P/err.log"; fail "process_GTEx.R on toy data"; }

Rscript - "$P" <<'RS' || fail "process_GTEx.R output checks"
args <- commandArgs(trailingOnly = TRUE); P <- args[1]
out <- read.delim(file.path(P, "out.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(identical(colnames(out), c("gene_id","gene_name","BRCA","GBMLGG","KIRC","LUAD","LUSC")))
stopifnot(identical(out$gene_id, c("ENSG01","ENSG02","ENSG03","ENSG04")))   # versions stripped, PAR_Y merged
# G1 by hand: cohort values BRCA 10, GBMLGG median(1,7,3)=3, KIRC 4, LUAD 2, LUSC 2 (Lung counted twice)
v <- c(10, 3, 4, 2, 2); z <- (v - mean(v)) / sd(v)
stopifnot(isTRUE(all.equal(unname(unlist(out[1, 3:7])), z, tolerance = 1e-6)))
# G2 not expressed anywhere -> SD 0 -> all NA
stopifnot(all(is.na(out[2, 3:7])))
# G3 duplicate IDs collapsed by max: BRCA 5, GBMLGG median(9,1,1)=1, KIRC 3, Lung 1 (min would give 1,0,1,0,0 - not affine-equivalent)
v <- c(5, 1, 3, 1, 1); z <- (v - mean(v)) / sd(v)
stopifnot(isTRUE(all.equal(unname(unlist(out[3, 3:7])), z, tolerance = 1e-6)))
# every gene's z-scores have mean 0 and SD 1 (where defined)
zz <- as.matrix(out[c(1,3,4), 3:7])
stopifnot(all(abs(rowMeans(zz)) < 1e-6), all(abs(apply(zz, 1, sd) - 1) < 1e-6))

# Oracle: PENNE's original v8 code path (06_Feature_preprocess_mrna_protein_levels.R),
# adapted only for file paths and column positions.
suppressMessages(library(dplyr))
info <- read.delim(file.path(P, "mapping.tsv")); colnames(info) <- c("V6p", "TCGA")
mrna.v8 <- read.delim(gzfile(file.path(P, "toy.gct.gz")), skip = 2)
mrna.v8$Name <- unlist(lapply(mrna.v8$Name, function(x) gsub("\\..*","",x)))
gene.info <- mrna.v8[,1:2]; gene.info <- gene.info[!duplicated(gene.info$Name),]
mrna.v8 <- mrna.v8 %>% group_by(Name) %>% summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))
mrna.v8 <- merge(gene.info, mrna.v8, by = "Name"); rownames(mrna.v8) <- mrna.v8$Name
m <- as.data.frame(t(mrna.v8[,-c(1,2)])); m$V8 <- rownames(m)
m <- merge(info, m, by.x = "V6p", by.y = "V8"); m <- m[complete.cases(m$TCGA),]
m <- m %>% group_by(TCGA) %>% summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)))
m <- as.data.frame(m); rownames(m) <- m$TCGA; m <- as.data.frame(t(m[,-1]))
m$Gene.stable.ID <- rownames(m)
m <- merge(mrna.v8[,1:2], m, by.x = "Name", by.y = "Gene.stable.ID")
m.scaled <- t(apply(as.matrix(m[,-c(1,2)]), 1, function(v) as.numeric(scale(v))))
colnames(m.scaled) <- colnames(m)[-c(1,2)]; rownames(m.scaled) <- m$Name
new <- as.matrix(out[, colnames(m.scaled)]); rownames(new) <- out$gene_id
new <- new[rownames(m.scaled), ]
stopifnot(identical(is.na(new), is.na(m.scaled)))
stopifnot(isTRUE(all.equal(unname(new[!is.na(new)]), unname(m.scaled[!is.na(m.scaled)]), tolerance = 1e-6)))
RS
pass "process_GTEx.R: hand-computed values, z-score properties, identical to PENNE code"

# mapped tissue absent from the GCT -> error naming the tissue (no silent drop)
printf 'Ovary\tOV\n' >> "$P/mapping.tsv"
if Rscript "$REPO/workflow/scripts/process_GTEx.R" --input_path_gct="$P/toy.gct.gz" \
    --input_path_mapping="$P/mapping.tsv" --output_path="$P/out2.tsv" 2> "$P/err2.log"; then
    fail "missing mapped tissue accepted"
fi
grep -q "Ovary" "$P/err2.log" || fail "error message does not name the missing tissue"
pass "process_GTEx.R: missing mapped tissue raises an error naming it"

# real mapping table shipped in config/ is well-formed
Rscript -e 'm <- read.delim(commandArgs(TRUE)[1]); stopifnot(identical(colnames(m), c("GTEx_tissue","TCGA_cohort")), nrow(m) == 40, length(unique(m$TCGA_cohort)) == 23, !anyNA(m), !any(grepl("\r", unlist(m))))' \
    "$REPO/config/mapping_GTEx_TCGA" 2>/dev/null || fail "config/mapping_GTEx_TCGA malformed"
pass "config/mapping_GTEx_TCGA: 40 tissue-cohort pairs, 23 cohorts, clean"

# =========================================================================================
# compute_GTEx_features.R
# =========================================================================================
F="$TMP/features"; mkdir -p "$F"
# Gene sets of the toy backbone used in protein_interactions_scripts.sh
printf 'binID\tgene_set\tgene\n' > "$F/gene_sets.tsv"
for row in "chr1_0_100 genes A" "chr1_0_100 genes B" "chr1_0_100 partners_trans C" "chr1_0_100 partners_trans D" \
           "chr1_0_100 partners_trans F" "chr1_0_100 all_int_trans C" "chr1_0_100 all_int_trans D" "chr1_0_100 all_int_trans F" \
           "chr1_100_200 genes C" "chr1_100_200 partners_trans A" "chr1_100_200 partners_trans B" "chr1_100_200 ppis_trans G" \
           "chr1_100_200 all_int_trans A" "chr1_100_200 all_int_trans B" "chr1_100_200 all_int_trans G" \
           "chr2_0_100 genes E" "chr2_0_100 genes G" "chr2_0_100 genes X1" "chr2_0_100 genes H" "chr2_0_100 ppis_trans C" "chr2_0_100 all_int_trans C" \
           "chr2_100_200 genes H" "chr2_100_200 ppis_trans D" "chr2_100_200 all_int_trans D"; do
    echo "$row" | tr ' ' '\t' >> "$F/gene_sets.tsv"
done
# Expression: F carries an outdated GTEx symbol ("Fold") that GENCODE renames to "F" via its
# Ensembl ID; H has no value (NA) and sits in a set with expressed genes (must be ignored, not
# propagate NA); E has no expression data at all.
printf 'gene_id\tgene_name\tBRCA\tLUAD\nENSGA\tA\t1\t0\nENSGB\tB\t2\t0\nENSGC\tC\t3\t0\nENSGD\tD\t-1\t0\nENSGF\tFold\t5\t0\nENSGG\tG\t0.5\t0\nENSGH\tH\tNA\t0\nENSGX\tX1\t10\t0\n' > "$F/expr.tsv"
printf 'chr1\t0\t10\tENSGA.1\tA\tprotein_coding\tgene\nchr1\t20\t30\tENSGF.3\tF\tprotein_coding\tgene\nchr1\t40\t50\tENSGF.3_PAR_Y\tF\tprotein_coding\tgene\n' > "$F/gencode.bed"

Rscript "$REPO/workflow/scripts/compute_GTEx_features.R" --input_path_gene_sets="$F/gene_sets.tsv" \
    --input_path_expression="$F/expr.tsv" --input_path_gencode="$F/gencode.bed" --GTEx_cohort=BRCA \
    --output_path="$F/out.tsv" 2> "$F/err.log" || { cat "$F/err.log"; fail "compute_GTEx_features.R on toy data"; }

Rscript - "$F" <<'RS' || fail "compute_GTEx_features.R output checks"
o <- read.delim(file.path(commandArgs(TRUE)[1], "out.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(colnames(o), c("binID","GTEx_genes","GTEx_partners_trans","GTEx_PPIs_trans","GTEx_all_int_trans")))
rownames(o) <- o$binID
# Hand-computed medians (F = 5 only reachable through the GENCODE ID harmonisation):
exp <- rbind(chr1_0_100   = c(1.5,  3,   NA,  3),
             chr1_100_200 = c(3,    1.5, 0.5, 1),
             chr2_0_100   = c(5.25, NA,  3,   3),
             chr2_100_200 = c(NA,   NA, -1,  -1))
got <- as.matrix(o[rownames(exp), -1])
stopifnot(identical(unname(is.na(got)), unname(is.na(exp))),
          isTRUE(all.equal(unname(got[!is.na(got)]), unname(exp[!is.na(exp)]))))
RS
grep -q "symbol changed by harmonisation: 1" "$F/err.log" || fail "harmonisation count not logged"
pass "compute_GTEx_features.R: hand-computed medians, NA handling, symbol harmonisation via GENCODE ID"

# hg19 / GENCODE lift37: back-mapped IDs carry a mapping-version suffix (ENSGF.3_1); the
# harmonisation must still match them, giving byte-identical output
sed -E 's/(ENSG[A-Z]+\.[0-9]+)/\1_1/' "$F/gencode.bed" > "$F/gencode_lift37.bed"
grep -q "ENSGF.3_1" "$F/gencode_lift37.bed" || fail "lift37 toy annotation not built"
Rscript "$REPO/workflow/scripts/compute_GTEx_features.R" --input_path_gene_sets="$F/gene_sets.tsv" \
    --input_path_expression="$F/expr.tsv" --input_path_gencode="$F/gencode_lift37.bed" --GTEx_cohort=BRCA \
    --output_path="$F/out_lift37.tsv" 2> "$F/err_lift37.log" || { cat "$F/err_lift37.log"; fail "compute_GTEx_features.R on lift37 IDs"; }
cmp -s "$F/out.tsv" "$F/out_lift37.tsv" || fail "lift37 ID suffix changes the output"
pass "compute_GTEx_features.R: GENCODE lift37 IDs (ENSG...N_1) harmonised identically"

if Rscript "$REPO/workflow/scripts/compute_GTEx_features.R" --input_path_gene_sets="$F/gene_sets.tsv" \
    --input_path_expression="$F/expr.tsv" --input_path_gencode="$F/gencode.bed" --GTEx_cohort=OV \
    --output_path="$F/out2.tsv" 2> "$F/err2.log"; then fail "unknown GTEx cohort accepted"; fi
grep -q "'OV' not found" "$F/err2.log" || fail "unknown cohort error message"
pass "compute_GTEx_features.R: unknown cohort raises an error"

# =========================================================================================
# Snakemake chain on toy inputs, real rule files:
# build_gene_interactors -> compute_interactor_features -> process_GTEx -> compute_GTEx_features
# =========================================================================================
W="$TMP/chain"; mkdir -p "$W/workflow/rules/GTEx" "$W/workflow/rules/protein_interactions" "$W/workflow/scripts" "$W/workflow/envs" \
    "$W/resources/GTEx/v8" "$W/resources/backbone" "$W/resources/reference_data/gencode/PC" \
    "$W/resources/protein_interactions/origin" "$W/config"
cp "$REPO"/workflow/rules/GTEx/*.smk "$W/workflow/rules/GTEx/"
cp "$REPO"/workflow/rules/protein_interactions/*.smk "$W/workflow/rules/protein_interactions/"
cp "$REPO/workflow/envs/r_plain.yml" "$W/workflow/envs/"
cp "$REPO"/workflow/scripts/{process_GTEx,build_gene_interactors,compute_interactor_features,compute_GTEx_features}.R "$W/workflow/scripts/"
Rscript -e '
P <- commandArgs(TRUE)[1]
human_proteincomplex <- data.frame(ComplexID = 1, ComplexName = "X", subunits.Gene.name. = "A;C", stringsAsFactors = FALSE)
save(human_proteincomplex, file = file.path(P, "resources/protein_interactions/origin/corum.RData"))
HS_interfaces_genenames <- data.frame(P1 = "A", P2 = "B", stringsAsFactors = FALSE)
save(HS_interfaces_genenames, file = file.path(P, "resources/protein_interactions/origin/insider.RData"))' "$W"
printf '#1.2\n3\t4\nName\tDescription\tBreast - Mammary Tissue\tLung\nENSGA.1\tA\t5\t1\nENSGB.1\tB\t1\t5\nENSGC.1\tC\t2\t3\n' | gzip -c > "$W/resources/GTEx/v8/gene_median_tpm.gct.gz"
printf 'GTEx_tissue\tTCGA_cohort\nBreast...Mammary.Tissue\tBRCA\nLung\tLUAD\n' > "$W/config/mapping_GTEx"
printf 'chr1\t0\t100\tchr1_0_100\tA\nchr1\t100\t200\tchr1_100_200\tB,C\n' > "$W/resources/backbone/hg38_binned_100_0_v45.mask.filtered.sorted.annotated.bed"
printf 'chr1\t0\t10\tENSGA.1\tA\tprotein_coding\tgene\nchr1\t100\t110\tENSGB.1\tB\tprotein_coding\tgene\nchr1\t120\t130\tENSGC.1\tC\tprotein_coding\tgene\n' \
    > "$W/resources/reference_data/gencode/PC/v45_gene.sorted.bed"
cat > "$W/Snakefile" <<'SMK'
mapping_file_GTEx_location = "config/mapping_GTEx"
CORUM_location = "resources/protein_interactions/origin/corum.RData"
InteractomeINSIDER_location = "resources/protein_interactions/origin/insider.RData"
include: "workflow/rules/protein_interactions/build_gene_interactors.smk"
include: "workflow/rules/protein_interactions/compute_interactor_features.smk"
include: "workflow/rules/GTEx/get_GTEx_median_TPM.smk"
include: "workflow/rules/GTEx/process_GTEx.smk"
include: "workflow/rules/GTEx/compute_GTEx_features.smk"
SMK
(cd "$W" && snakemake data/GTEx/features/hg38_binned_100_0_v45_BRCA_v8 data/protein_interactions/hg38_binned_100_0_v45 \
    -s Snakefile --cores 1 > "$W/run.out" 2>&1) || { cat "$W/run.out"; fail "Snakemake chain on toy inputs"; }
Rscript - "$W" <<'RS' || fail "Snakemake chain output checks"
W <- commandArgs(TRUE)[1]
g <- read.delim(file.path(W, "data/GTEx/features/hg38_binned_100_0_v45_BRCA_v8"))
i <- read.delim(file.path(W, "data/protein_interactions/hg38_binned_100_0_v45"))
# BRCA z: cohort values A (5,1), B (1,5), C (2,3) -> z_BRCA = +0.707, -0.707, -0.707
z <- 1 / sqrt(2)
stopifnot(identical(g$binID, c("chr1_0_100", "chr1_100_200")),
          isTRUE(all.equal(g$GTEx_genes, c(z, -z))),
          isTRUE(all.equal(g$GTEx_partners_trans, c(-z, z))),     # A <-> C are complex partners across bins
          isTRUE(all.equal(g$GTEx_PPIs_trans, c(-z, z))),         # A <-> B are PPI partners across bins
          identical(i$n_complex_partners_trans, c(1L, 1L)), identical(i$n_PPIs_trans, c(1L, 1L)),
          isTRUE(all.equal(i$complex_density, c(1, 0.5))))
RS
pass "Snakemake chain (4 rules, real rule files) produces the expected per-bin features"

echo "All GTEx tests passed."