rule compute_matrices_hic:
    output:
        "resources/hic/files/cooler/{SRA_id}/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.cool"
    input:
        # NOTE (revised Option B design): cload's --targets is now the raw, pre-filter,
        # sorted backbone (continuous, no gaps, from sort_raw_backbone_hic.smk) instead of a
        # separately-generated `cooler makebins` bin set - reuses the pipeline's own existing
        # bin construction instead of introducing a second, parallel one. cload discards any
        # extra columns from the targets file (confirmed directly - its own bins() only ever
        # has chrom/start/end), so the binID column here doesn't do anything for cload itself,
        # but having it be the SAME file elsewhere in the pipeline uses keeps everything
        # consistent. Blacklisted/gap regions are excluded LATER, explicitly, in
        # process_bins_hic.R.
        raw_sorted_backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.sorted.bed",
        pairs="resources/hic/files/pairs/{SRA_id}/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}.filtered.dedup.sorted.clipped.pairs"        
    log:
        "logs/compute_matrices_hic/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{SRA_id}.log"
    conda:
        "../../envs/hic_processing_env.yml"
    shell:
        """
        cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \
        --no-symmetric-upper \
        --assembly {wildcards.assembly_code} \
        {input.raw_sorted_backbone} \
        {input.pairs} \
        {output} 2>{log}
        """