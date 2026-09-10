# Bin-level genomics pipeline — Phase 1: backbone/core annotations module

This package contains the **first phase** of the backtracking work on the original
Snakemake pipeline: the "backbone/core annotations" module (chrom sizes, blacklist,
genomic bin generation, GENCODE annotation, gene distances/overlaps, GC content,
centromeres/telomeres), rebuilt into the standard Snakemake folder structure, fixed
wherever a bug was found and **verified with real tests** (not just by inspection),
and equipped with a 3-tier GitHub Actions CI.

All modules not yet migrated (epigenomes, IntOGen, DepMap, COSMIC, GDC, fragile
sites, ohnologs/paralogs, reference genome, ATAC-seq, Hi-C, RepliSeq) are present in
the Snakefile as commented-out `include:` lines marked with
`# TODO(phase>1 - not yet migrated/tested)`, so they can be re-activated one at a
time in later phases without losing track of what's left.

## Structure

```
workflow/
  snakefile_bin_level          # production Snakefile (phase 1: only the backbone module is active)
  snakefile_helpers/           # helpers.py, helpers.smk (unchanged)
  rules/
    chrom_sizes/
    blacklisted_regions/
    backbone/
    gencode/
    genes/
    GC/
    centromeres_and_telomeres/
    DepMap/
    IntOGen/
  scripts/                     # fixed R scripts + scripts/helpers/
  envs/                        # conda env yml files used by this module
config/
  config.yaml                  # production config (original, values unchanged)
  mapping_tissue_tumor_type, mapping_gender, sequencers_info,   # <- 1 typo fixed here
  BREAST_selected_samples_atacseq, BREAST_samples_to_exclude_atacseq,
  MCF10A_selected_samples_hic, HMEC_selected_samples_repliseq
.tests/
  unit/test_*.sh                    # fast unit tests on individual R scripts (no Snakemake)
  integration/config/config.yaml    # test config using a fake "hg00" assembly
  integration/data/                 # toy fixtures replacing the real downloads
.github/workflows/ci.yml       # 3 jobs: production dry-run, unit tests, toy integration test
```

## Bugs found and fixed (all verified with real tests, not just by inspection)

1. **`get_actual_chromsizes_from_backbone.R`** — guaranteed crash: `group_by`/
   `summarise` used without `library(dplyr)`. *Fix*: namespaced calls
   `dplyr::group_by(...) %>% dplyr::summarise(...)`.

2. **`get_TSSs.R`** — coordinate bug: the GFF3 (1-based) → BED (0-based) conversion
   was only applied on the `-` strand; TSSs on the `+` strand were shifted by
   **+1bp**. *Fix*: the conversion is now applied on both strands.

3. **`split_gtf.R`** — systemic bug: **no** feature (gene/exon/UTR/transcript, for
   PC/lncRNA/miRNA/rRNA) underwent the 1-based→0-based conversion: every BED file
   produced had its start shifted by 1bp. *Fix*: `start <- start - 1` before writing,
   for every sub-table.

4. **`get_telomeres_proxy.R`** — the "telomere-start" rows (position 0) were
   hard-coded to `chr1`...`chr22`: with datasets other than the full hg38 this
   produces ghost rows, and more importantly **chrX/chrY always lost their start
   telomere** (only the end one, generically inherited from chrom_sizes, remained).
   *Fix*: chromosomes are now derived from the input file instead of being
   hard-coded.

5. **`config/sequencers_info`** — typo `HiSeq3000` (no underscore) vs. the query
   `"HiSeq_3000"` in the Snakefile. Combined with a weakness in `get_is_patterned`
   (substring comparison instead of equality, with `all([])` vacuously true), the
   mismatch caused **`is_patterned_hic` to always evaluate to `True`** regardless of
   the actual value in the file. *Fix*: fixed the typo in the config data. **Note**:
   the weakness in `get_is_patterned` (in `helpers.py`) is still present (it works by
   coincidence when the strings match exactly) — worth revisiting when the
   Hi-C/ATAC-seq/RepliSeq module gets migrated.

6. **`get_cytobands.smk` / `clean_chrom_sizes.smk`** — unquoted `grep -v ["_XMY"]`:
   confirmed with a test that the resulting behaviour (excludes alt-contigs
   containing `_` as well as chrX/chrY/chrM) is the intended one (**confirmed with
   the user: chrX/chrY must always be excluded**), but rewritten explicitly and
   readably (`grep -vE '_|X|M|Y'`).

7. **`get_cytobands.smk`** — the declared `.txt.gz` output was not actually gzipped
   (plain text was written under a `.gz` name). *Fix (user's choice, simpler than my
   initial proposal)*: dropped `.gz` from the output filename instead of
   re-compressing — the content was already plain text, so the name now reflects
   reality. As a consequence, `get_centromeres_coords.smk` uses plain `grep` again
   (not `zgrep`).

8. **`generate_genes_on_bins_obj.smk`** — **a DAG bug that had never actually been
   run**: the declared input did not include `_v{gencode_version}` in the filename,
   but `annotate_backbone.smk` always produces that output **with** that suffix →
   guaranteed `MissingInputException`. Discovered via a real Snakemake dry-run
   (never previously performed). *Fix*: added `_v{gencode_version}` to the rule's
   input/output/log (and updated the target in `rule all`).

9. **Fragility identified (not fixed, only documented)**: if a gene category
   (PC/lncRNA/miRNA/rRNA) has **no** feature at all in the GTF for the chosen
   assembly/version, `bedtools closest -b <empty file>` produces rows with fewer
   columns than expected, and `process_dist_to_genes.R`/`process_num_overlaps.R`
   crash on `colnames<-` (length mismatch). With a full GENCODE annotation on a real
   genome this probably never happens, but it would be worth making the R scripts
   more robust (e.g. an explicit column-count check before renaming) at some point —
   left untouched here since it would change real behaviour, not just syntax.

10. **Directory-structure mismatch**: the `include:`/`conda:`/`Rscript
    workflow/scripts/...` statements in the `.smk` files assumed a
    `workflow/{rules,scripts,envs}/` + `config/` layout that was "flat" in the
    uploaded material. Rebuilt from scratch in this package.

## Phase 2 — DepMap module (essential genes)

Added the "essential genes" branch (`filter_CRISPR_*`, `get_essential_genes`,
`find_coords_essential_genes`, `*_essentiality_score`, `get_bin_level_essentiality`),
re-activated in `rule all` (this was the pipeline's original target). Bugs found and
fixed (verified with realistic DepMap test data, including the column-name
"mangling" R performs on headers like `"GENEA (1)"`):

11. **`get_model_average_essentiality_score.R`** — guaranteed crash when **a single
    gene** passes all essentiality filters (far from a rare case with realistic
    cutoffs): `crispr.f[,-1]` decays from a data.frame to a plain vector when only
    one column remains, and `colMeans()` errors out (`'x' must be an array of at
    least two dimensions`). *Fix*: `drop = FALSE` on both subsetting operations.
    Tested with both 1 and 2+ essential genes (non-regression).

12. **`get_essential_genes.smk`** — the log path pointed, by copy-paste mistake, to
    `logs/filter_CRISPR_probabilities/...` instead of `logs/get_essential_genes/...`.
    *Fix*: corrected the path.

13. **`find_coords_essential_genes.smk`** — the rule inside it was named
    `find_coord_essential_genes` (singular), inconsistent with the (plural) file
    name and with the other rules in the module. *Fix*: renamed for consistency
    (verified: no other part of the code referenced it by name).

No bugs found in: `filter_CRISPR_models.R`, `parse_depmap_colnames.R` (R's
column-name mangling of DepMap headers like `"GENE (12345)"`, including ones with
hyphens such as `HLA-DRB1`, is handled correctly — verified with a dedicated test),
`get_essential_genes.R` (including its already well-implemented handling of genes
that are essential in the probabilities file but missing from the scores file — an
explicit warning message, behaviour tested).

## Phase 3 — IntOGen module (cancer genes, TSG/OG)

Added the "distance to cancer driver genes" branch (oncogenes/tumor suppressors),
re-activated in `rule all`. Found **one critical bug**, verified with a real test:

14. 🔴 **`process_IntOGen.R` — locale-dependent gene-role mapping, a silent and
    potentially serious bug**. The original code determined the readable label
    (`OG`/`ambiguous`/`TSG`) by positionally matching `sort(unique(ROLE))` against
    `c("OG", "ambiguous", "TSG")`, assuming the alphabetical order was
    `["Act", "ambiguous", "LoF"]`. **Verified with a real test that in the POSIX/C
    locale (the default on many HPC clusters, Docker containers, and CI runners) the
    actual order is instead `["Act", "LoF", "ambiguous"]`** — result: genuine tumor
    suppressor genes (`LoF`) were labelled `ambiguous`, and genuinely ambiguous genes
    were labelled `TSG`. A systematic, silent swap (no error, no warning) that would
    have corrupted every downstream TSG/OG feature. *Fix*: replaced with an explicit
    name-based mapping (`c(Act="OG", LoF="TSG", ambiguous="ambiguous")`), independent
    of locale/sort, with an explicit check that makes the script fail loudly
    (instead of silently mismapping) should IntOGen ever introduce a ROLE value not
    covered by the mapping.

No other bugs found in the rest of the module (`sort_cancer_genes`,
`split_cancer_genes`, `get_dist_closest_OG/TSG`, `process_dist_to_cancer_genes`) —
verified with an end-to-end test on realistic IntOGen data.

## Phase 4 — epigenomes/Roadmap module

Added the "chromatin-state per bin" branch, which downloads Roadmap Epigenomics
ChromHMM segmentations for a given tissue/sample type, computes for each bin the
number of base pairs covered by each chromatin state, and averages that across all
matching epigenomes. Re-activated in `rule all`. Found **one bug that would silently
corrupt results**, verified with a real test:

15. 🔴 **`aggregate_epigenomes.R` — NA-fill happens after averaging instead of
    before, silently zeroing out real signal**. Different Roadmap epigenomes can
    have different sets of chromatin states present in their segmentation (a state
    that never occurs anywhere in one epigenome simply has no column after
    `compute_state_length`'s `pivot_wider`). After `dplyr::bind_rows`, such a missing
    column becomes `NA` for that epigenome's rows. The original code computed
    `mean()` (default `na.rm = FALSE`) across epigenomes **before** replacing `NA`
    with `0` — meaning that as soon as a *single* epigenome lacked a given state
    column, `mean()` returned `NA` for *every bin* of that state, which was then
    turned into `0` — discarding real, non-zero signal from every other epigenome
    that *did* have that state.

    **Why this is a real bug and not just a design choice**: a state column missing
    from one epigenome means that epigenome genuinely has 0bp of that state
    everywhere (ChromHMM assigns exactly one state to every genomic position, so an
    entirely-absent state really is "0 across the whole genome" for that epigenome,
    not "unknown"). So when averaging across epigenomes, that epigenome's
    contribution for that state/bin should be `0`, exactly like any other bin where
    the state simply doesn't occur — it should count *in* the average, not void it.
    Confirmed with the pipeline's author: with real, whole-genome Roadmap data this
    situation (a state entirely absent from one epigenome but present in others of
    the same tissue) is expected to be rare, which is why it had not surfaced before,
    but the fix is correct regardless of how often it triggers.

    *Fix*: moved the `NA -> 0` fill to right after `bind_rows`, before
    `group_by`/`summarise`/`mean` — a missing state column for one epigenome now
    correctly counts as "0bp of that state" in the average, instead of voiding the
    average entirely. Verified with a toy scenario where one out of two epigenomes
    lacks a state column: the average now correctly comes out to `mean(value, 0)`
    instead of `0`.

Also renamed the rule in `process_epigenomes.smk` from `process_epigenome` (singular)
to `process_epigenomes` (plural), for consistency with the file name and the rest of
the module (verified: not referenced by name anywhere else), and fixed the
corresponding log path (`logs/process_epigenomes/...`, was `logs/process_epigenome/...`).

No bugs found in `compute_state_length.R` (bin × chromatin-state overlap matrix,
verified with a toy ChromHMM segmentation spanning two states over two bins) or in
the network-scraping rules (`get_mapping_EID_tt`, `get_epigenome_paths`,
`get_epigenome_files`) — **note**: these last three could only be reviewed by
inspection, not verified end-to-end against the real external sources (a public
Google Sheet and the WashU Roadmap file server), since that would require live
network access this environment doesn't have. Their column-position assumptions
(`cut -f 2,4,6,15,17,18` on the Google Sheet export) and regex-based HTML scraping
are inherently more fragile than the rest of the pipeline and are worth an occasional
manual sanity-check against the live sources.

## Phase 5 — COSMIC module (coding/non-coding mutation overlaps)

Added the branch computing, per bin, the number of base pairs covered by
COSMIC-catalogued coding and non-coding mutations (filtered by tumor
site/histology and by patient gender), re-activated in `rule all`. Unlike the
previous phases, **no functional bug was found** in this module — but two
non-obvious points are worth documenting explicitly, both verified with tests:

- **Coordinate handling in `get_COSMIC_mutations_coordinates.smk` is correct, but
  relies on an undocumented convention that could easily be broken by a future
  "fix".** COSMIC's own `GENOME_START`/`GENOME_STOP` columns are 1-based inclusive
  (a single-base substitution has `GENOME_START == GENOME_STOP`). This rule copies
  them as-is into a file named `.bed`, **without** the usual 1-based→0-based
  conversion. Verified with a real test that this is nonetheless correct: every
  downstream consumer (`compute_coding_mutations_overlap.R`,
  `compute_noncoding_mutations_overlap.R`, `get_coding_mutations_description.R`)
  reads this file with plain `read.table()`/`read.delim()` and builds a `GRanges`
  directly from these raw numbers — never via `rtracklayer::import.bed()` (which
  would assume real 0-based BED and add 1 to the start, shifting everything by
  1bp). The backbone, by contrast, genuinely is 0-based BED and *is* read via
  `import.bed()`, which correctly restores its 1-based coordinates for the overlap
  computation. The two omissions cancel out correctly today, but this is fragile:
  if anyone ever reads this specific file via `import.bed()`/`bedtools` in a
  coordinate-sensitive way, it will silently be off by 1bp. Documented with an
  explicit comment in the rule so a future maintainer doesn't "fix" the missing
  conversion (which would break things) or switch a consumer to `import.bed()`
  (which would also break things).

- **`filter_classification_COSMIC.R` assumes `COSMIC_parameters` has exactly one
  underscore** (a single-word site + a single-word histology, e.g.
  `breast_carcinoma`, parsed positionally as `parts[1]`/`parts[2]`). Verified with a
  test that this silently gives the wrong site/histology as soon as either one
  contains an underscore itself (e.g. real COSMIC site names like
  `central_nervous_system` or `haematopoietic_and_lymphoid_tissue`). Not fixed,
  since there is no unambiguous general solution without changing the config
  format (e.g. splitting `COSMIC_parameters` into two separate config fields); the
  tumor type currently in production (`BRCA` → `breast_carcinoma`) is unaffected
  since neither component contains an underscore. Documented with a comment in the
  script and a dedicated (expected-to-currently-pass) test that reproduces the
  limitation, so it is caught immediately if a future tumor type triggers it.

## Phase 6 — GDC_CNA module (copy-number amplification/deletion frequency)

Added the branch computing, per bin, the fraction of patients with an amplification
or a deletion (from TCGA copy-number segmentation data downloaded via `gdc-client`),
re-activated in `rule all`. Two small issues found and fixed (verified with tests):

16. **`get_GDC_CNA_data.smk` — `touch {output}` ran *before* the actual
    `gdc-client download` call**, not after. Snakemake's default `set -euo pipefail`
    means a genuinely failing `gdc-client` call would still abort the script (and
    Snakemake's built-in failed-job cleanup would remove the touched file), so this
    likely never caused silent false successes in practice — but the ordering was
    needlessly confusing and offered no safety margin against a `gdc-client` version
    that logs errors without a non-zero exit code (a real possibility with some
    download clients that don't fail hard on partial errors, which this codebase
    cannot verify without live GDC access). *Fix*: moved the `touch` to run only
    after the download command, for a clear "done means done" semantics with no
    downside.

17. **`compute_ampl_overlap.smk`** had a garbled log path
    (`logs/compute_GDcompute_ampl_overlapC_CNA_overlap/...`, clearly a
    find-and-replace mistake) — harmless (Snakemake would just create that
    oddly-named directory) but confusing for anyone debugging a failed run. *Fix*:
    corrected to `logs/compute_ampl_overlap/...`.

No other bugs found: `get_samples_GDC.R` (cohort/gender filtering, dropping
patients with an ambiguous >1 sample/file match), `filter_GDC_CNA_data.R`
(probe-count and `|segment_mean|` thresholding, plus the diagnostic plots), the
amplification/deletion split (`grep -v "-"` / `grep "-"` on the whole line — safe
here specifically because the upstream `cut` already discards every column that
could contain a literal hyphen other than `segment_mean` itself, verified with a
test), and `process_GDC_CNA_overlap.R` (frequency = overlap count / number of
patients, plus the landscape plots) all behaved correctly on realistic toy TCGA-like
data.

**Testing note**: `gdc-client` needs live network access to the GDC API, which this
environment (and GitHub Actions, by default) does not have. The toy integration test
therefore runs the pipeline in two steps: first building the real, filtered manifest
(exercising `get_samples_GDC.R` and the manifest-filtering logic for real), then
marking the `gdc-client` download step as already done and pointing it at
pre-placed toy segmentation files, before letting the rest of the pipeline
(filtering, splitting, overlap, frequency computation, plotting) run for real.

## Phase 7 — fragile_sites, ohnologs_and_paralogs, gene_density, haploinsufficiency

These four small modules had clearly never been executed even once (unlike every
previous phase): every single rule/script had at least one bug that would have
caused an immediate crash or silent data corruption. Re-activated in `rule all`.

18. 🔴 **Systemic bug across 5 rules**: `compute_and_process_ohnologs.smk`,
    `compute_and_process_paralogs.smk`, `process_haploinsufficiency.smk`,
    `compute_density_cancer_genes.smk` and `compute_density_essential_genes.smk` all
    referenced the annotated backbone as
    `resources/backbone/{...}.mask.filtered.sorted.annotated.bed`, but
    `annotate_backbone` always produces it **with** `_v{gencode_version}` in the
    name (the exact same class of bug found and fixed in Phase 1 for
    `generate_genes_on_bins_obj.smk`) — guaranteed `MissingInputException` on all
    five. *Fix*: added `_v{gencode_version}` to each rule's input. Since ohnologs,
    paralogs and haploinsufficiency didn't have `gencode_version` anywhere in their
    own output name (unlike the two gene_density rules, which already did), their
    output filenames now include it too — a correct, not just cosmetic, change: the
    result genuinely depends on which GENCODE version annotated the backbone, so it
    belongs in the name of a target that was never wired into `rule all` before this
    phase (i.e. this isn't "renaming an established output", it's "correctly naming
    it for the first time").

19. **`process_dist_to_fragile_sites.R`** — `colnames(processed_table) <- c("chr",
    "start", "end", "dist.to.closest.FGS")` was missing `"binID"` (4 names for a
    5-column table). R does not error on this — it silently pads with `NA` and
    shifts every subsequent name by one, so `binID` and the real distance column
    both ended up wrongly named, causing a crash further down (`undefined columns
    selected`). *Fix*: added the missing name.

20. ~~🔴 **`compute_and_process_ohnologs.R`** — `mean()`/`sum()` referenced a missing
    column~~ **RETRACTED — this was my own mistake, not a bug in your code.** I
    initially built a toy `OhnologList.RData` fixture for testing that omitted the
    `n_ohnolog_mmpaper` column, without knowing the real object's schema, and
    concluded from the resulting crash that the column was missing from the
    pipeline's logic. After you shared the real `OhnologList.RData`, it turned out
    `ohnologs[["mmpaper"]]` already contains a clean, pre-computed integer
    `n_ohnolog_mmpaper` column — the original script was correct all along. I have
    reverted my incorrect "fix" (which had added a redundant, and potentially
    subtly different, recomputed column) and re-verified the original code against
    a toy fixture matching the real schema: it works correctly with no changes.
    Sorry for the churn — this is a good example of exactly the kind of
    external-data-format assumption I'd flagged as unverifiable without your data,
    and in this one case my guess when building the toy fixture was wrong.

**Verified against real data you provided** (previously flagged as "assumed, not
verifiable without live/licensed access"):
- `OhnologList.RData` / `ParalogList.RData`: confirmed real schemas as above.
  `ParalogList.RData`'s `n_paralogs` is stored as **character**, confirming the
  explicit `as.numeric()` cast in `compute_and_process_paralogs.R` is necessary
  (not just defensive).
- GDC `clinical.tsv`: confirmed `project.project_id`, `cases.case_id`,
  `cases.submitter_id`, `demographic.gender` are all present exactly as
  `get_samples_GDC.R` expects. Also confirmed the file is a wide, denormalized
  export with one row per diagnosis/treatment record (the same patient appears
  many times) — exactly what the script's `unique()` on the four
  identifying/gender columns is designed to collapse.
- GDC sample sheet: real headers use spaces (`"Case ID"`, `"Sample ID"`, `"File
  ID"`, `"File Name"`), not dots. Confirmed with a direct test that R's
  `read.delim()` (default `check.names = TRUE`) automatically mangles these to
  `Case.ID`/`Sample.ID`/`File.ID`/`File.Name` — exactly matching what the script
  references. No bug.
- COSMIC classification file: confirmed `COSMIC_PHENOTYPE_ID`, `PRIMARY_SITE`,
  `PRIMARY_HISTOLOGY` are present (among other columns the script doesn't use).
  Also confirmed real `PRIMARY_HISTOLOGY`/`PRIMARY_SITE` values that **do** contain
  underscores (e.g. `haematopoietic_neoplasm`, `haematopoietic_and_lymphoid_tissue`)
  — concretely validating bug #15's `COSMIC_parameters` limitation as a real,
  not just hypothetical, risk for tumor types other than the current
  `breast_carcinoma`.

All toy fixtures in `.tests/integration/data/` and the unit tests above were
updated to match these confirmed real schemas.

21. 🔴 **`process_gene_density.R`** — two bugs: (a) `TSGs <-
    read.table(input_path_OGs, ...)` read the **OG** file path instead of
    `input_path_TSGs` — the TSG output silently contained OG data under a TSG
    label; (b) `write.table(..., row.names = T, ...)` (every other script in this
    codebase uses `row.names = F`) added an unlabelled leading row-index column,
    misaligning the header against the data for any downstream reader expecting a
    fixed column count. Both confirmed and fixed with a test.

22. 🔴 **`process_haploinsufficiency.R`** — two bugs: (a) same missing-`library(dplyr)`
    issue as Phase 1's `get_actual_chromsizes_from_backbone.R` (unqualified
    `select`/`group_by`/`summarise`, guaranteed crash); (b) the annotated backbone
    bed file has **no header row** (it's raw `bedtools`/awk output), but was read
    with `header = T` — silently consuming whichever bin happened to be first in
    the file as a "header" and losing it from the analysis entirely. Both confirmed
    with a test (the dropped bin was reproducible and specific to file order) and
    fixed (namespaced `dplyr::` calls; `header = F`).

23. **`process_haploinsufficiency.smk`** — the `shell:` block referenced
    `{input.path_haplo}` and `{input.path_backbone}`, but the rule's `input:` block
    declares them as `haplo` and `backbone` (no `path_` prefix) — this would have
    raised an `AttributeError` the moment Snakemake tried to build the command.
    Also had spaces around every `=` in the argument list (`--haplo_thr =
    {wildcards.haplo_thr}`), which would have broken `R.utils::commandArgs`'s
    `--key=value` parsing even if the attribute names had been correct, and was
    missing the `2>{log}` redirect entirely (so the log file, though declared,
    would never actually have been written). *Fix*: corrected the attribute names,
    removed the spaces, added the log redirect.

24. **`process_gene_density.smk`** — `Rscript process_gene_density.R` was missing
    the `workflow/scripts/` path prefix used by every other rule in the codebase —
    would have failed with "cannot open file" from the repository root. *Fix*:
    added the prefix (and, while at it, added the `2>{log}` redirect this rule was
    also missing).

25. **Config gap**: `haplo_thr` (the pLI threshold for calling a gene
    haploinsufficient) did not exist anywhere in `config.yaml` — consistent with
    `process_haploinsufficiency` never having been wired into `rule all` before.
    **Added `haplo_thrs: [0.9]`** (the standard ExAC/gnomAD literature convention
    for "haploinsufficient"). **Please confirm this matches your intended
    threshold** — it was not something migrated from your existing config, since
    none existed.

26. **Observation, not a bug**: `resources/fragile_sites/origin/fragile_sites_li_{assembly_code}.bed`
    has no `config`-driven location and no rule that produces it (unlike every
    other external dataset in this pipeline, which is either downloaded by a rule
    or has a `..._locations` config entry). It's expected to already exist at that
    fixed path. Worth confirming this is intentional (a manually-curated,
    non-redistributable file from the Li et al. fragile-sites publication) rather
    than a rule that's simply missing.

## Phase 8 — genome module (reference genome + aligner indexes)

Added the branch downloading the reference genome FASTA and building `bwa`/`bowtie2`
indexes from it, re-activated in `rule all`. No functional bugs found — the only
change was cosmetic, for consistency with the fix applied to `get_chrom_sizes.smk`
in Phase 1:

27. **`get_reference_genome.smk`** used the bare `{assembly_code}` inside the
    `shell:` block instead of `{wildcards.assembly_code}`. This works today only
    because the pipeline enforces exactly one `assembly_code` (Snakemake resolves
    unqualified names in shell strings against the local/global Python scope too),
    but it's the same fragile pattern flagged for `get_cytobands.smk` in Phase 1.
    *Fix*: made it explicit (`{wildcards.assembly_code}`), matching the sibling
    `index_reference_genome_bowtie2.smk`, which already did this correctly.

`bwa index` and `bowtie2-build` were verified end-to-end with a toy multi-chromosome
FASTA, producing valid index files. This module has no R scripts and needs no
config changes (the reference genome URL is built entirely from `assembly_code`,
like the other UCSC downloads in Phase 1).

## Phase 9 — dataset module (explore_datatypes, combine_dataset; smooth_dataset deprioritized)

Per your instructions: `explore_datatypes` fully migrated and activated;
`combine_dataset` fully reviewed, fixed, and unit-tested but **not yet wired into
`rule all`** (it needs `data/atacseq/...` and `data/repliseq/...`, neither of which
exists yet — activate once ATAC-seq and RepliSeq are migrated); `smooth_dataset`
deprioritized per your explicit note that smoothing isn't currently needed.

28. **`explore_datatypes.smk` did not exist** (the file was completely empty — no
    rule at all) even though `workflow/scripts/explore_datatypes.R` was already a
    complete, working 118-line script. Built the rule from scratch by reading every
    argument the script actually uses, so **all 6 of its outputs** are covered (GC
    landscape, GC distribution, samples barplot, genes barplot, mutations barplot,
    coding-mutation-descriptions barplot) — not just a subset. Verified end-to-end
    with toy data for every input.

29. **`explore_datatypes.R`** — two argument keys didn't match their variable names:
    `input_path_GC <- args$input_path_ampl` (a GC-content input read from a key
    named after amplifications — almost certainly copy-pasted from a CNA-related
    script) and `output_path_num_mutations_barplot <- args$output_path_mutations_barplot`.
    Neither would have crashed anything (the script is internally consistent as
    long as whatever calls it uses the same argument keys), but since I'm writing
    the calling rule for the first time anyway, renamed both for clarity rather
    than perpetuate a confusing key name with no established caller to stay
    consistent with.

30. 🔴 **The rule as originally sketched grouped outputs with mismatched wildcard
    sets** (e.g. a GC-content plot using only `assembly_code`/`window_mode`/etc.
    alongside a samples plot using only the COSMIC/DepMap/GDC parameters) — this is
    rejected outright by Snakemake itself ("Not all output, log and benchmark files
    of rule ... contain the same wildcards"), caught by a real dry-run, not by
    inspection. *Fix*: every output now carries the full combined wildcard set.

31. 🔴 **`resources/GC/...` vs `data/GC/...` — wrong input path (my own mistake
    while writing the rule, caught by a real test run, not a pre-existing bug)**.
    `explore_datatypes.R` reads the GC-content file with `read.table(header=F)` and
    assigns 5 column names (`chr, start, end, binID, GC.content`) — the *raw*
    `get_GC_content` output format. I initially wired the rule to
    `data/GC/...` (the *final*, 2-column `process_GC_content` output) instead,
    which crashed immediately (`'names' attribute [5] must be the same length as
    the vector [2]`). Fixed to point at `resources/GC/...`.

32. 🔴 **`combine_dataset.R` had three bugs that would each have crashed the script,
    all confirmed with a 20-bin toy dataset covering every one of its ~26 inputs**:
    - `output_path_bool_multicollinearity <- args$output_path_bool_multicollinearity`
      was commented out (apparently as a side effect of commenting out the unused
      PCA block above it), but the resulting variable **is** used later
      (`writeLines(..., con = output_path_bool_multicollinearity)`) — guaranteed
      `object not found` crash. Fixed by uncommenting just that one assignment
      (the PCA-plot one stays commented, consistent with the PCA computation
      itself being commented out too).
    - The final `ggsave()` for the correlation heatmap referenced
      `output_path_corr` (never assigned — the real variable is
      `output_path_corr_plot`) — same class of bug, fixed by using the correct
      variable name.
    - The `ggsave()` for the chromatin-state 0-vs-non-0 boxplot passed
      **the output file path** as the plot object (`ggsave(output_path_..., output_path_...)`)
      instead of the actual computed plot (`boxplots_0s_chromstates`) — fixed.

33. 🔴 **`combine_dataset.smk` — wrong input for epigenomes, discovered while
    building the corrected rule, verified with a toy test**. It pointed at
    `resources/epigenomes/checks/processing/...`, an intermediate file that (per
    `process_epigenomes.smk`, Phase 4) contains only a **path string** to the real
    per-epigenome table, not an actual `binID`-keyed table. `read.delim()`-ing it
    with `header=T` produces a 0-row data.frame with one bizarrely-named column
    (confirmed with a direct test) — `merge(..., by="binID")` would fail
    immediately since no such column exists. Fixed to point at
    `data/epigenomes/...` (the real, final output of `aggregate_epigenomes`,
    already verified correct in Phase 4).

34. **`combine_dataset.smk` was also missing the `input_path_haploinsufficiency`
    input entirely**, even though `combine_dataset.R` reads
    `args$input_path_haploinsufficiency` and merges it in like every other
    feature (confirmed: without it, `read.delim(NULL, ...)` would fail). Added,
    pointing at Phase 7's `data/haploinsufficiency/...` output (which also meant
    adding `{haplo_thr}` to `combine_dataset`'s own wildcard set, since it wasn't
    there before).

35. **`combine_dataset.smk`'s ohnologs/paralogs inputs were missing
    `_v{gencode_version}`**, the same systemic issue fixed across 5 other rules in
    Phase 7 — these two outputs gained that wildcard when re-activated, and this
    consumer needed updating to match.

Also added all of `combine_dataset.R`'s outputs that the original rule was missing
entirely (it only declared 4 of the script's 9): `output_path_zero_frac_plot`,
`output_path_boxplots_scaled`, `output_path_boxplots_0s_non0s_muts`,
`output_path_boxplots_0s_non0s_chromstates`, `output_path_bool_multicollinearity`
— and fixed `output_path_NA_plot` → `output_path_NA_frac_plot` (another
key-name mismatch that would have left `args$output_path_NA_frac_plot` as `NULL`).

**`smooth_dataset.smk`** is confirmed to be a genuine empty stub (no input, output,
or shell at all) despite `workflow/scripts/smooth_dataset.R` existing with real
logic — the rule itself was simply never written. Left untouched, per your
instruction that smoothing isn't currently needed; noted clearly in the Snakefile
so this isn't mistaken for "not yet reviewed" like the other TODOs.

## Phase 10 — RepliSeq module (SRA download through replication timing)

The full NGS branch: SRA prefetch/fastq-dump, `trim_galore`, `bwa mem` alignment,
canonical-chromosome filtering, `picard MarkDuplicates`, mapq filtering, blacklist
removal, `bedtools coverage`, per-sample coverage filtering, and finally the
cross-sample (ES/LS/G1) replication-timing computation. Re-activated in `rule all`
and verified end-to-end with a toy paired-end FASTQ dataset (systematic tiling
coverage across a toy genome, not random sampling — an earlier attempt with random
reads happened to leave zero bins covered by all three sample types, which turned
out to be a toy-data-scale artifact, not a pipeline bug, once diagnosed).

This module had the highest concentration of confirmed bugs of any phase so far —
several were caught only by really executing the full chain, not by inspection.

36. 🔴 **`process_repliseq.R` — the final `data/repliseq/...` output could be empty,
    isolated to the mask-overlap step added very recently (per your confirmation)
    and not yet exercised**. `binID`s have their chromosome's `"chr"` prefix
    stripped (`generate_backbone.R`: `binID <- paste0(gsub("chr","",chr), ...)`), so
    splitting `binID` back into `chr`/`pos` (`separate_wider_delim`) gives a bare
    number ("1", not "chr1"). This `chr` column is then used directly as
    `seqnames.field` for `GenomicRanges::makeGRangesFromDataFrame()`, while the
    backbone **mask** file (read straight from disk a few lines earlier) uses real
    chromosome names ("chr1"). `GenomicRanges::findOverlaps()` requires exact
    seqname matches, so this mismatch means the final overlap always returns zero
    hits — confirmed directly (`replication_timing$chr` prints `"1"`,
    `mask$chr` prints `"chr1"`) and by running the whole rule twice, before/after
    the fix, in the real pipeline.

    **Important scope correction** (found while double-checking this with you): the
    3 diagnostic plots (distribution/landscape/raw-vs-smoothed) are generated
    *before* this step, and derive their own `chr`/`pos` directly from `binID` via
    `sub()` — a code path that was never affected by this bug. So the plots always
    showed real data, exactly as you observed; only the final *data file*
    (`output_path_data`, the one `combine_dataset.R` actually reads) was at risk,
    and only once the mask-overlap code was added. *Fix*: re-add the `"chr"` prefix
    to the split-off chromosome column before building the GRanges. Verified with a
    toy test that reproduced the mismatch and confirmed the fix (0 rows before, 4
    rows after, using the exact same script and data).

37. 🔴 **`find_contiguous_subsets.R`** — two compounding bugs: (a) the classic
    `seq(0)` gotcha (`seq(0)` evaluates to `c(1, 0)`, a length-2 *decreasing*
    sequence, not an empty one), and (b) an `if / else if` chain that could never
    add *both* the opening and closing boundary when there was exactly **one**
    break in a chromosome's bin sequence (a single break is simultaneously "the
    first" and "the last" break, but the `else if` meant only one branch could ever
    fire). Confirmed with a systematic test across 0, 1, 2, 3 and 4 breaks: **the
    crash is specific to exactly 1 break** — 0 breaks works (by an accidental
    cancellation of the two bugs), and 2 or more breaks always worked correctly
    even in the original, unfixed code.

    **Why this likely never surfaced in practice**: a real ENCODE blacklist has
    dozens of discrete regions per chromosome, so filtering real bins against it
    almost always produces either 0 breaks (rare) or several (the typical case) —
    landing on *exactly* 1 is a narrow edge case that real, genome-wide blacklist
    data is unlikely to hit, even though it's easy to construct with toy data (as I
    initially did). The fix is a zero-risk safety net either way: verified as a
    pure non-regression for 0/2/3/4-break inputs, in addition to fixing the 1-break
    crash. *Fix*: explicit early return for zero breaks, `seq_len()` instead of
    `seq()` throughout, and independent (not `else if`) checks for the first/last
    break so both boundaries always get added regardless of how many breaks there
    are.

38. **`process_repliseq.R` reads each sample-type's coverage file with
    `read.table(..., header = T)`**, but `filter_coverage_repliseq.R` writes it
    with `col.names = F` (no header at all, confirmed) — the first bin of every
    single ES/LS/G1 file was silently consumed as a fake header and dropped.
    Confirmed with a toy test (2 rows in, 1 row read back). *Fix*: `header = F`.

39. **`process_repliseq.R`**: `pattern <- samples[samples$sample_type == sample_type,
    "SRA_id"]` can be a multi-element vector if more than one replicate SRA run
    exists for the same sample_type (not the case in your current sample sheets,
    which each have exactly one ES/LS/G1 run — confirmed by inspection — so this
    hasn't triggered yet). `grep(pattern = <vector>, ...)` silently uses only the
    first element, with just a warning, dropping every other replicate with no
    visible error. *Fix*: rather than guess how multiple replicates should be
    combined (sum coverage? average CPM? — a scientific decision, not mine to
    make), the script now explicitly checks that exactly one file matches per
    sample_type and **fails with a clear, actionable error** if not, so this would
    surface immediately rather than silently dropping data if you ever add a
    replicate.

40. 🔴 **`compute_library_size_repliseq.smk` — race condition, confirmed by running
    the real pipeline with parallel jobs (`--cores 2`)**. The rule's only declared
    Snakemake input was a generic "check" marker file; the actual
    `*.noblack.sorted.bam` files it reads are found via a `bash` glob *inside the
    shell script*, never declared as real inputs. Those bam files are `temp()` —
    Snakemake deletes them as soon as their *declared* consumers
    (`compute_coverage_repliseq`, `index_noblacklist_bam_repliseq`) are done, which
    can happen (and did, in testing) before this rule gets scheduled, leaving the
    glob matching nothing and the rule crashing with `cut: ... No such file or
    directory`. *Fix*: declared the actual per-sample bam files as a proper input
    (reusing the same helper function the pipeline already has for exactly this
    file set), which tells Snakemake to keep them alive until this rule has run
    too. Also replaced the runtime glob with Snakemake's own resolved file list.

41. ~~`compute_coverage_repliseq.smk` used the non-mask backbone~~ **RETRACTED — this
    was my mistake, not a bug.** I initially flagged the use of the non-mask
    backbone (`....filtered.sorted.bed`, instead of the `....mask.filtered.sorted.bed`
    variant used by every other analysis branch) as an inconsistency and "fixed" it
    to use the mask. **Confirmed with you that this was intentional**: replication-timing
    smoothing (the `loess()` step in `process_repliseq.R`) needs spatial continuity
    across the *whole* genome to work correctly — restricting the coverage
    computation to an arbitrary masked subset of bins would be scientifically wrong
    for this specific step, unlike every other feature in this pipeline (which are
    independent per-bin values with no spatial-continuity requirement). Reverted to
    the original, correct behavior.

42. Minor, cosmetic fixes for consistency with the rest of the codebase: a stray
    log-path name (`get_distribution_counts_repliseq` → `filter_coverage_repliseq`,
    matching the actual rule), two missing `2>{log}` redirects
    (`filter_coverage_repliseq.smk`, `process_repliseq.smk` — the log files were
    declared but nothing was ever written to them), and a rule renamed
    (`compute_library_size` → `compute_library_size_repliseq`) to match its file
    name, consistent with earlier phases.

**Observation, not changed**: `filter_bam_repliseq.smk` explicitly keeps `chrX`/`chrY`
(`chr{1..22} chrX chrY`) when filtering to canonical chromosomes, which looks like
it contradicts the "chrX/chrY always excluded" convention confirmed with you in
Phase 1. In practice this has no effect on final results: the backbone never has
chrX/chrY bins (confirmed, `clean_chrom_sizes` excludes them), so any chrX/chrY
reads kept here simply never overlap anything in the later `bedtools coverage`
step — harmless, just a few reads' worth of wasted disk/CPU. Left as-is since
changing it wouldn't change any result, only worth knowing about.

**Also verified, not a bug**: `filter_bam_repliseq.smk`'s
`$(echo chr{1..22} chrX chrY)` relies on bash brace expansion, which does *not*
work under `/bin/sh`/`dash` — an initial toy-test failure here was a red herring
caused by testing outside of Snakemake. Confirmed directly with a minimal
Snakemake rule that Snakemake always invokes `shell:` blocks via `bash`, so this
works correctly in the real pipeline.

## Phase 11 — ATAC-seq module (alignment through CNVkit-scaled counts, cohort exploration,
##             bin-level DARs; Genrich peak-calling chain not active per your confirmation)

By far the highest bug density of any phase so far — several rules had clearly
never been executed even once, including one rule file that wasn't even valid
Python syntax. Re-activated in `rule all` (alignment → dedup → blacklist removal
→ CNVkit tumor+normal-panel scaling → quantification+PCA+SVA → cohort exploration
→ bin-level DARs). The peak-calling branch (`call_peaks_bysample_atacseq`,
`call_peaks_conditionwise_atacseq`, `process_peaks_atacseq`, `quantification_atacseq`)
is confirmed **not currently active** and was left out of `rule all`, per your
explicit confirmation; its one guaranteed-crash syntax bug was still fixed for
correctness (see #52), since leaving broken-Python in an unrelated file is a
needless landmine even if nothing currently includes it.

Validated end-to-end with a real toy cohort (8 samples: 4 tumor + 4 NAT, 2
studies, matched donors across conditions) run through the actual tools
(`trim_galore`, `bowtie2`, `samtools`, `picard`, `bedtools`, and — installed
specifically for this phase — `CNVkit` and its R dependency `DNAcopy`), not just
isolated script calls: **231/231 steps** complete for the entire pipeline
(Phases 1 through 11) run together from a clean checkout.

43. 🔴 **`call_peaks_conditionwise_atacseq.smk` was not valid Python** — every
    string in its `output:` block was missing its closing quote
    (`pq_vals_path="...pq_vals,` never closes the string), which would raise an
    immediate `SyntaxError` the moment Snakemake tried to `include:` this file —
    even *before* considering whether the rule itself is active. Also had
    unescaped `{print $1}` / `${id}` inside an embedded `awk`/loop in the shell
    block, which Snakemake's own `{...}` formatting would have choked on next.
    *Fix*: closed every string, escaped the braces (`{{print $1}}`, `${{id}}`).

44. 🔴 **`fix_mate_information_atacseq.smk`'s output filename didn't match what 3
    downstream rules expected** (`name_sort_final_bam_atacseq`,
    `compute_insert_size_metrics_atacseq`, `compute_coverage_atacseq` all expected
    no `{SRA_id}` prefix in the filename; the rule produced one). I initially
    "fixed" this by removing the prefix to match those 3 — which promptly broke
    CNVkit (see #46) since it needs the *opposite* convention. Settled by
    restoring the `{SRA_id}` prefix (the original, correct form) and updating the
    3 other consumers to match instead — see #46 for why this direction is the
    one that's actually required.

45. **`explore_cohort_atacseq.smk`** had two independent bugs against the same
    input: the shell referenced `{input.path_library_size}` (singular) while the
    rule declared `path_library_sizes` (plural) — a guaranteed `AttributeError`;
    and even the *path itself* didn't match any real output of
    `get_counts_atacseq.smk` (missing the `wCNAs/` directory component and the
    `{SVs}` wildcard entirely). *Fix*: corrected the attribute name and the path
    (added `{SVs}` to this rule's own wildcard set, required since Snakemake
    needs every output/log to share the same wildcards as the newly-added path).

46. 🔴 **CNVkit requires each bam passed to `cnvkit.py batch` to have a unique
    *filename* (not just a unique path)** — it derives its internal sample ID
    from the basename alone. Removing the `{SRA_id}` prefix (per #44's initial,
    wrong-direction fix) made every sample's `.fixed.bam` share the exact same
    filename, differing only by directory — confirmed with a real run:
    `Duplicate sample ID 'hg00_20.fixed'`, since CNVkit's normal-panel/leave-one-out
    steps read *multiple* samples' bams in a single invocation, unlike the other
    3 consumers which only ever handle one bam at a time. Restoring the prefix
    (see #44) fixes this; also had to propagate the correct convention back
    through `run_CNVkit_atacseq(_normal).smk`, `normal_paths_atacseq` in the
    Snakefile, `get_LOO_normal_bam_paths_atacseq` in `helpers.py`, and the `.cns`
    path in `process_CNVkit_output_atacseq.smk` (CNVkit names its own output
    after the input bam's basename).

47. 🔴 **`process_CNVkit_output_atacseq.smk` depended on a file
    (`resources/atacseq/files/cnvkit/{SRA_id}/{assembly_code}_{mapq_thr}.fixed.cns`)
    that no rule declared as a Snakemake output anywhere** — it's a genuine side
    effect of `cnvkit.py batch`, but Snakemake had no way to build a DAG edge to
    it, confirmed by a real dry-run failing with `MissingInputException`.
    Declaring it as a formal wildcarded output wasn't possible either (CNVkit's
    filename has no room for `window_mode`/`window_size`/`step`, which Snakemake
    requires to match across every output of the *producing* rule — see #48).
    *Fix*: added a `get_CNVkit_check_atacseq` input function (in `helpers.py`)
    that picks whichever of `run_CNVkit_atacseq`'s or
    `run_CNVkit_atacseq_normal`'s "check" marker applies to the given `SRA_id`,
    giving Snakemake a real dependency edge; the `.cns` path itself moved to
    `params:` (Snakemake doesn't require every referenced path to be a tracked
    input, only genuine `input:` entries — this is intentional, not a workaround).

48. **`process_CNVkit_output_atacseq.smk` was also missing the
    `window_mode`/`window_size`/`step` wildcards entirely**, even though CNVkit's
    `--targets` (and therefore its result) genuinely depends on which backbone
    binning was used — the same class of fix as the `_v{gencode_version}` issue
    from Phase 7. *Fix*: added them; `get_CNVkit_scaling_factors_atacseq.smk`'s
    reference to this file's output was updated to match.

49. 🔴 **`get_counts_atacseq.smk` passed two whole directories
    (`resources/atacseq/files/coverage`, `.../cnvkit`) as formal Snakemake
    `input:`** — neither directory is itself declared as anyone's output (they're
    just the accumulated per-sample files from earlier rules), so a real dry-run
    failed with `MissingInputException` on the bare directory path. *Fix*: moved
    both to `params:`, relying on the rule's two genuine "check" file inputs
    (already present) to guarantee every sample's file exists first.

50. 🔴 **Discovered that `workflow/snakefile_helpers/helpers.py` is imported as a
    real Python module** (`from snakefile_helpers.helpers import *`), which gives
    it its own separate global namespace — functions defined there cannot see
    variables computed in the main Snakefile just because they're assigned
    later in file order. This broke my own new `get_CNVkit_check_atacseq`
    function immediately (`NameError`), and, on inspection, turned out to already
    be true of the *pre-existing* `get_LOO_normal_bam_paths_atacseq` (present
    since before this phase) — it references `SRA_normal_ids_atacseq`, which was
    never visible to it either. Since nothing had ever actually invoked this
    function before ATAC-seq was activated, this had been a **latent, undetected
    bug the whole time**. *Fix*: after computing `SRA_tumor_ids_atacseq` /
    `SRA_normal_ids_atacseq` in the main Snakefile, they're now explicitly
    injected into the `helpers` module's namespace
    (`snakefile_helpers.helpers.SRA_tumor_ids_atacseq = ...`), making both
    functions work correctly.

51. 🔴 **Race condition in `run_CNVkit_atacseq.smk`, confirmed with a real
    8-sample run**: the normal-panel bam paths (`normal_paths_atacseq`) were only
    ever a `params` (a plain string built once at parse time), never a real
    Snakemake `input:` — so nothing guaranteed those files existed before this
    rule ran. With real (non-toy-scale) parallelism this surfaced as
    `FileNotFoundError` reading a normal sample's bam that another, still-running
    job hadn't finished writing yet — same class of bug as
    `compute_library_size_repliseq`'s race condition in Phase 10. *Fix*: declared
    the normal-panel bams as a genuine input (a lambda building the file list),
    so Snakemake now waits for all of them before running this rule. (The sibling
    rule, `run_CNVkit_atacseq_normal.smk`, already did this correctly for its own
    leave-one-out panel — only the tumor-side rule had the bug.)

52. 🔴 **`get_DARs.smk`'s shell command never passed any `--output_path_*`
    arguments to `get_DARs.R` at all** — the rule declares 10 named outputs, but
    the shell block stops after the last `--input_path_*` argument. Every
    `output_path_*` variable inside the script would have been `NULL`, and
    `write.table(..., file = NULL)` / `ggsave(NULL, ...)` fail immediately. This
    is why my own earlier manual testing of `get_DARs.R` (Phase "continue ATAC-seq
    with the same rigor" session) never caught it — I always passed the full,
    correct argument list by hand; only running it through the *real* rule
    surfaced the missing arguments. *Fix*: added all 10 missing
    `--output_path_*` arguments to the shell command.

53. **Removed genuinely dead code**: `do_get_counts_analysis.R` computed
    `n.sv <- num.sv(count.f.norm[,-1], mm1, method="leek")`, but `n.sv` is never
    referenced again anywhere in the file — `sva()` is always called with the
    explicit, user-provided `SVs` count instead. Confirmed unused by inspection.
    Incidentally, this exact call is also numerically fragile with few
    features/samples (`Error in eigen(wm) : infinite or missing values in 'x'`,
    hit repeatedly while building a toy cohort small enough to run quickly in
    CI) — removing the dead computation removes that fragility too, with zero
    change in behavior since the result was never used.

54. Two `stopifnot()` safety nets added, mirroring the pattern from Phase 9-10:
    `do_get_counts_analysis.R`'s `all(counts$binID == scaling_factors$binID)` and
    `do_DA_analysis.R`'s `names(libsize) == colnames(counts)` were both bare
    boolean expressions whose result was never checked. As with the earlier
    instances, both `counts`/`scaling_factors` and `library_sizes`/`counts`
    likely stay aligned in practice (they trace back to the same underlying
    `list.files()` ordering), so this is a defensive, zero-cost addition rather
    than a confirmed active bug.

**Retracted, not a bug** (confirmed with you): I initially flagged the "with
CNAs"/"without CNAs" (`out.list`/`out.list.woCNAs`) labeling in
`get_counts_atacseq.R` as inverted, assuming "woCNAs" meant "the analysis that
doesn't use CNA information at all". You clarified the intended reading is "the
signal with the CNA-driven distortion removed" — i.e. `woCNAs` = *corrected*,
exactly matching what the code does. No change made.

**Confirmed with you, not changed**: `run_CNVkit_atacseq(_normal).smk` use the
non-mask backbone (`....filtered.sorted.bed`) as CNVkit's `--targets`, unlike the
mask variant this pipeline uses everywhere else for per-bin analyses. Given the
precedent from Phase 10 (`compute_coverage_repliseq`), this is very plausibly
intentional for the same reason (CNVkit's own copy-number segmentation benefits
from genome-wide continuity, not an arbitrary masked subset) - flagging for
awareness rather than assuming either way.

**CNVkit scaling factor is in linear space, confirmed with you**:
`process_CNVkit_output_atacseq.smk` computes both the raw log2 ratio and
`2^log2ratio` from CNVkit's `.cns` output, but originally kept the *raw log2*
value in the final `.processed` file - while `get_counts_atacseq.R` then
*divides* raw counts by this value to correct for copy-number-driven signal
inflation. Dividing by a log2 ratio (which can be negative, or near zero) would
produce negative or wildly inflated "corrected" counts - mathematically
inconsistent with a multiplicative scaling factor. Confirmed with you that
scaling must be in linear space; fixed to keep `2^log2ratio` instead.

## Phase 12 — Hi-C module (SRA download through eigenvector/compartment calling), and
##             combine_dataset fully unblocked

The last remaining NGS module. Re-activated in `rule all`: alignment (`bwa mem -5SP` +
`samblaster` + `samtools`) → `pairtools` parse/sort/dedup/select → clipping to the
backbone's actual chromosome extents → per-sample `cooler cload` → merge → ICE
balancing → z-score normalization → eigenvector decomposition. Validated end-to-end
with a real toy cohort (2 samples) run through the actual tools (`bwa`,
`samblaster`, `pairtools`, `cooler` — all installed specifically for this phase):
**34/34 steps** for the Hi-C branch alone, and **257/257 steps for the entire
pipeline (Phases 1 through 12, including `combine_dataset`) run together from a
single clean checkout.**

Since both of `combine_dataset`'s remaining blockers (ATAC-seq DARs, RepliSeq) were
completed in Phases 10-11, and Hi-C — per your explicit request — is now
incorporated into it too, **`combine_dataset` is fully active as of this phase.**

55. 🔴🔴 **The single most consequential bug found in this entire project.**
    `process_bins_hic.R` used `merge(x=backbone, y=bins, by=c("chr","start","end"))`
    to attach `binID`s to cooler's own bin table. R's `merge()` **sorts its result
    alphabetically by the `by` columns by default** — and alphabetically,
    `"chr10"` sorts *before* `"chr2"`, even though it comes after it in every real
    assembly's actual chromosome order. Confirmed directly: a toy input with
    chr1/chr2/chr10 in the correct order came back from `merge()` with chr10
    sandwiched between chr1 and chr2. This matters enormously here specifically
    because `bins` **is cooler's own bin table — its row order *is* the contact
    matrix's row/column order** — and `process_hic.py` assumes `bins[i]`
    corresponds to matrix row `i` when assembling the final output. Left
    unfixed, this would have **silently assigned every eigenvector to the wrong
    bin for roughly half the genome**, on any real assembly with 10 or more
    numbered chromosomes (i.e. every real human or mouse genome) — with no error,
    no warning, just quietly wrong compartment calls. *Fix*: replaced `merge()`
    with an explicit, order-preserving `match()`-based lookup that keeps `bins`'
    original row order exactly. Verified with the same chr1/chr2/chr10 toy case:
    output order now matches cooler's input order precisely.

    Worth noting: `combine_dataset.R`, elsewhere in this same codebase, calls
    `merge(..., by="binID", sort = F)` (explicitly disabling the same alphabetical
    sort) and *additionally* re-sorts numerically at the end
    (`dataset[order(as.numeric(coords[,1]), as.numeric(coords[,2])),]`) — clear
    evidence the pipeline's author was well aware of this exact class of bug and
    handled it correctly elsewhere. This reinforces that `process_bins_hic.R`'s
    plain `merge()` was very likely an isolated oversight, not a deliberate choice.

    **You pushed back on this finding twice, appropriately** — first asking
    whether it was really a bug given you'd used this code without noticing
    anything wrong, then sharing an unrelated real-data excerpt to check. That
    excerpt turned out to be a *different* merge (`by="bin"`, a genuine key-based
    join elsewhere in your own downstream code, not part of this pipeline) that
    is **not** affected by this issue - key-based merges correctly re-associate
    each row's data regardless of final row order, unlike the *positional*
    correspondence `process_bins_hic.R`/`process_hic.py` rely on. Whether this bug
    has actually affected your past results depends on whether your real Hi-C
    backbone/chrom-sizes files list chromosomes past chr9 in natural order
    (chr1, chr2, ..., chr9, chr10, ...) — if so, this fix matters for you; if your
    analyses so far only involved fewer than 10 chromosomes, it would not have
    triggered yet, but the fix protects against it going forward regardless.

56. 🔴 **`deduplicate_pairs_hic.smk` ran `pairtools dedup` in the background and
    waited for it via a `kill -0`/`sleep 1` polling loop, instead of running it
    synchronously.** Confirmed directly (a minimal bash reproduction) that this
    pattern **completely hides a failing background job's exit status** — even
    with Snakemake's default `set -euo pipefail`, a background job that fails
    leaves the polling loop (and therefore the whole script) reporting success,
    since `set -e` only ever applies to synchronous commands and nothing in the
    loop ever inspected the background job's actual exit code. **You confirmed
    you'd specifically tested and needed this background+heartbeat pattern**
    (unspecified reason — possibly a job-scheduler timeout/heartbeat
    requirement), so rather than replacing it with a plain synchronous call, the
    fix keeps the exact same background+heartbeat structure and adds a
    `wait $DEDUP_PID` after the polling loop, which captures and correctly
    propagates the real exit code. Verified both directions with a direct test:
    a failing background job is now caught (script exits non-zero), and a
    succeeding one still reports success as before (non-regression).

57. **`index_bam_hic.smk`** had a malformed log path:
    `logs/index_bam_hic/.{assembly_code}_{SRA_id}log` — a stray leading dot (making
    it a hidden file) and a missing `.` before `log` (concatenating directly onto
    the SRA_id with no separator). *Fix*: `logs/index_bam_hic/{assembly_code}_{SRA_id}.log`.

58. **`merge_cooler.smk`'s glob for collecting per-sample `.cool` files
    (`resources/hic/files/cooler/*/*.cool`) did not filter by
    assembly_code/mapq_thr/window_mode/window_size/step** — with the pipeline's
    current single-value config for these parameters this is harmless (there is
    only ever one combination on disk at a time), but if ever run with multiple
    values (every other equivalent parameter in this codebase supports config
    lists), it would silently merge `.cool` files from different, incompatible
    parameter combinations together. *Fix*: the glob now matches only the
    specific combination being requested.

59. **Discovered while doing the initial file inventory**: two Hi-C rule files
    (`balance_cooler.smk`, `merge_cooler.smk`) don't have "hic" anywhere in their
    filename, so an initial `ls *hic*.smk`-style search missed them entirely —
    caught by re-searching for "hic" in file *contents* instead of names, and
    confirmed necessary since `normalize_matrix_hic.smk` depends on a
    `.balanced.cool` file that only `balance_cooler.smk` produces.

60. **`process_hic.py` crashed (`shape mismatch`) whenever fewer than 10 valid
    (non-all-zero-row) bins remained** after filtering — the code unconditionally
    requested the "top 10" eigenvectors regardless of how many were actually
    available. Hit with sparse toy data (a handful of valid bins); real
    genome-scale Hi-C data has far more than 10 valid bins per chromosome, so
    this is a defensive fallback (`min(10, n_valid)`) rather than a change to
    normal behavior.

61. **Added Hi-C to `combine_dataset`, per your explicit request** (it had no Hi-C
    input at all before this phase). Thanks to `combine_dataset.R`'s existing
    design — it auto-collects every `input_path_*` variable in scope and merges
    each one in by `binID` — adding Hi-C required exactly one new line in the
    script (`input_path_hic <- args$input_path_hic`) plus the corresponding new
    input/argument in the `.smk` rule; no other code changes were needed. Verified
    end-to-end: the final feature matrix correctly includes `eigen_1`...`eigen_10`.
    Note (confirmed with you): the separate, manual `if(add_HiC_repliseq ==
    TRUE){...}` merge you showed from your own downstream code is *older, unrelated*
    code, not part of this pipeline — it isn't reviewed or affected by this change.

62. **`combine_dataset.R`'s pairwise `cor.test()` loop (the multicollinearity
    diagnostic) crashed the entire script — before the feature matrix was even
    written — whenever any two columns had zero variance**, confirmed with the
    real toy Snakemake pipeline (several near-constant toy features triggered
    `Error: not enough finite observations`). Since this diagnostic check runs
    *before* `write.table()` for the actual dataset, a single degenerate column
    pair could block the core deliverable entirely. *Fix*: wrapped `cor.test()` in
    `tryCatch` (a degenerate pair is skipped, recorded as `NA`, rather than
    crashing everything), and added `na.rm = TRUE` to the final
    `any(cor.test.result$FDR < FDR_thr)` check (otherwise a single `NA` — from a
    skipped pair — would make the multicollinearity flag itself `NA` instead of
    `TRUE`/`FALSE`). No change in behavior for real, large-scale data, where
    near-zero-variance features are not expected.

## What is NOT (yet) tested in CI

- **`get_GC_content.R`** requires the Bioconductor package `BSgenome.Hsapiens.UCSC.*`
  specific to the real assembly (~800MB for hg38) — not testable with a fake
  assembly. It is validated only for: (a) syntax (`parse()` with no errors), (b) its
  downstream counterpart `process_GC_content.R`, tested with toy data. A real
  end-to-end test of this rule would need a dedicated, slower CI job, with a real
  assembly (e.g. hg38) and BSgenome package caching.
- **`explore_datatypes`** depends on `resources/GC/...`, the *raw* `get_GC_content`
  output — same limitation as above, so this rule is validated only via its
  dedicated unit test (`test_explore_datatypes_scripts.sh`), not end-to-end in the
  toy integration test.
- **`combine_dataset`** is not part of the toy integration test at all: it isn't
  wired into `rule all` yet (blocked on ATAC-seq/RepliSeq, see Phase 9), and is
  instead fully exercised by its own unit test
  (`test_combine_dataset_scripts.sh`), covering all ~27 inputs and 9 outputs with a
  20-bin toy dataset. **Update (Phase 12): `combine_dataset` is now fully active
  in `rule all` and validated end-to-end via the real toy Snakemake pipeline**
  (see Phase 12) — this note is kept for historical context on how it was first
  validated before all its blockers were resolved.
- **All modules are now migrated and tested**, except:
  - `rules/dataset/smooth_dataset.smk` — confirmed to be a genuine empty stub (no
    input/output/shell), deprioritized per your explicit instruction.
  - ATAC-seq's peak-calling/Genrich branch (`call_peaks_bysample_atacseq`,
    `call_peaks_conditionwise_atacseq`, `process_peaks_atacseq`,
    `quantification_atacseq`) — confirmed not currently active and intentionally
    excluded from `rule all`.
  - Hi-C's `run_HiCQC` — requires a manually pre-created, non-Snakemake-managed
    conda environment (`conda activate hic_qc`), which this environment cannot
    replicate; not on the critical path to the final `.eigen` output, so this
    doesn't block anything, but it hasn't been executed here (same category of
    limitation as `get_GC_content`'s BSgenome dependency).
- **`rules/ML/`** is deliberately excluded from tracking (no `include:`, not even as
  a TODO): it concerns the machine-learning models, out of scope for this work.

## Running the tests locally

```bash
# unit tests on individual R scripts (fast, no Snakemake)
bash .tests/unit/test_backbone_scripts.sh
bash .tests/unit/test_depmap_scripts.sh
bash .tests/unit/test_intogen_scripts.sh
bash .tests/unit/test_epigenomes_scripts.sh
bash .tests/unit/test_cosmic_scripts.sh
bash .tests/unit/test_gdc_scripts.sh
bash .tests/unit/test_fragile_sites_scripts.sh
bash .tests/unit/test_ohnologs_paralogs_scripts.sh
bash .tests/unit/test_gene_density_haploinsufficiency_scripts.sh
bash .tests/unit/test_genome_scripts.sh
bash .tests/unit/test_explore_datatypes_scripts.sh
bash .tests/unit/test_combine_dataset_scripts.sh
bash .tests/unit/test_repliseq_scripts.sh
bash .tests/unit/test_atacseq_scripts.sh
bash .tests/unit/test_hic_scripts.sh

# dry-run on the production DAG (requires empty placeholders for the DepMap/IntOGen/COSMIC
# source files, since they are now part of the active rule all - see ci.yml)
mkdir -p resources/DepMap/models resources/DepMap/probabilities resources/DepMap/scores resources/IntOGen/IntOGen-Drivers-20240920
touch resources/DepMap/models/Model.csv resources/DepMap/probabilities/CRISPRGeneDependency.csv resources/DepMap/scores/CRISPRGeneEffect.csv
touch resources/IntOGen/IntOGen-Drivers-20240920/Compendium_Cancer_Genes.tsv
mkdir -p resources/COSMIC/classification/Cosmic_Classification_Tsv_v103_GRCh38 resources/COSMIC/samples/Cosmic_Sample_Tsv_v103_GRCh38 resources/COSMIC/mutations/coding/Cosmic_GenomeScreensMutant_Tsv_v103_GRCh38 resources/COSMIC/mutations/noncoding/Cosmic_NonCodingVariants_Tsv_v103_GRCh38
touch resources/COSMIC/classification/Cosmic_Classification_Tsv_v103_GRCh38/Cosmic_Classification_v103_GRCh38.tsv resources/COSMIC/samples/Cosmic_Sample_Tsv_v103_GRCh38/Cosmic_Sample_v103_GRCh38.tsv resources/COSMIC/mutations/coding/Cosmic_GenomeScreensMutant_Tsv_v103_GRCh38/Cosmic_GenomeScreensMutant_v103_GRCh38.tsv resources/COSMIC/mutations/noncoding/Cosmic_NonCodingVariants_Tsv_v103_GRCh38/Cosmic_NonCodingVariants_v103_GRCh38.tsv
mkdir -p resources/GDC/CNAs/clinical.cohort.2026-01-28 resources/GDC/CNAs/files/dummy
touch resources/GDC/CNAs/clinical.cohort.2026-01-28/clinical.tsv resources/GDC/CNAs/gdc_sample_sheet.2026-01-28.tsv resources/GDC/CNAs/gdc_manifest.2026-01-28.170055.txt resources/GDC/CNAs/files/dummy/dummy.seg.v2.txt
mkdir -p resources/fragile_sites/origin resources/ohnologs_and_paralogs/origin resources/haploinsufficiency
touch resources/fragile_sites/origin/fragile_sites_li_hg38.bed resources/ohnologs_and_paralogs/origin/OhnologList.RData resources/ohnologs_and_paralogs/origin/ParalogList.RData resources/haploinsufficiency/haploinsufficiency.csv
snakemake -n -s workflow/snakefile_bin_level --configfile config/config.yaml

# real execution on toy data (fake "hg00" assembly)
mkdir -p resources && cp -r .tests/integration/data/resources/. ./resources/
snakemake -s workflow/snakefile_bin_level --configfile .tests/integration/config/config.yaml \
  --cores 2 data/distances/centromeres/hg00_binned_200_0  # ...and the other targets, see ci.yml
```

## Next steps (future phases)

All modules covered by this project's scope are now migrated, fixed, and tested
(Phases 1-12). The only remaining, deliberately-out-of-scope items are:
`dataset/smooth_dataset` (an empty stub, deprioritized per your instruction),
ATAC-seq's peak-calling/Genrich branch (confirmed not currently active), and
Hi-C's `run_HiCQC` (needs a manually pre-created conda environment this setup
can't replicate, not on the critical path). If any of these become relevant
later, the same phase-by-phase process (fixes + toy tests + dedicated CI) applies.
