rule process_CNVkit_output_atacseq:
    output:
        "resources/atacseq/files/cnvkit/{SRA_id}/{SRA_id}_{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}.processed"
    input:
        # NOTE (fix): added window_mode/window_size/step wildcards - CNVkit's --targets is the
        # backbone for a specific binning, so different binnings genuinely produce different
        # CNVkit results and need to be tracked as distinct outputs (same reasoning as the
        # gencode_version fixes elsewhere in this codebase). Also added the "check" dependency
        # (see get_CNVkit_check_atacseq in helpers.py) so Snakemake can actually build the DAG:
        # the .cns file itself is a side effect of cnvkit.py batch, never declared as anyone's
        # Snakemake output, so without this the real dry-run failed with MissingInputException.
        check=get_CNVkit_check_atacseq
    params:
        # the .cns path itself is a `params`, not a formal `input`, on purpose: Snakemake
        # requires every `input:` path to have either a producing rule or to already exist as
        # a source file, but this one can't cleanly be declared as a wildcarded Snakemake
        # output (CNVkit's own filename has no room for window_mode/window_size/step - adding
        # them would conflict with Snakemake's same-wildcards-per-rule-output requirement, see
        # run_CNVkit_atacseq.smk). The "check" input above already guarantees cnvkit.py batch
        # has completed - and therefore that this file physically exists - before this rule runs.
        # NOTE (fix 26/10/2026): the .cns is now read from the per-resolution directory written by run_CNVkit_atacseq /
        # run_CNVkit_atacseq_normal (before, all resolutions shared one directory and overwrote each other's files).
        cns="resources/atacseq/files/cnvkit/{SRA_id}/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}/{SRA_id}_{assembly_code}_{mapq_thr}.fixed.cns"
    log:
        "logs/process_CNVkit_output/{SRA_id}_{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}.log"
    shell:
        """
        cat {params.cns} | \
        awk -v OFS='\t' '{{print $1,$2,$3,$5,2^$5,$6,$7,$8}}' |
        tail -n+2 | \
        awk -v OFS='\t' '{{print $1,$2,$3,$5}}' \
        > {output} 2>{log}
        """