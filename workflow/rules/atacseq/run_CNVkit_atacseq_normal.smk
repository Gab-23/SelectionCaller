rule run_CNVkit_atacseq_normal:
    output:
        "resources/atacseq/files/cnvkit/checks/{SRA_normal_id}_{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_normal"
    input:
        metadata=selected_samples_location_atacseq,
        normal_bam="resources/atacseq/files/bam/{SRA_normal_id}/{SRA_normal_id}_{assembly_code}_{mapq_thr}.fixed.bam",
        loo_bams=get_LOO_normal_bam_paths_atacseq,
        backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed",
        genome="resources/reference_data/genome/{assembly_code}.fa"
    log:
        "logs/run_CNVkit_normal/{assembly_code}_{window_mode}_{window_size}_{step}_{SRA_normal_id}_{mapq_thr}.log"
    params:
        num_cores_atacseq=num_cores_atacseq,
        path_dir="resources/atacseq/files/cnvkit/",
        loo_paths=lambda wildcards, input: " ".join(input.loo_bams)
    conda:
        "../../envs/cnvkit_env.yml"
    shell:
        """
        touch {output}

        cnvkit.py batch {input.normal_bam} -n {params.loo_paths} \
        --fasta {input.genome} \
        --method wgs \
        --target-avg-size {wildcards.window_size} \
        --targets {input.backbone} \
        --processes {params.num_cores_atacseq} \
        -d {params.path_dir}/{wildcards.SRA_normal_id}
        """