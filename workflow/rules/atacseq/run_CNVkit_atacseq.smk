rule run_CNVkit_atacseq:
    output:
        "resources/atacseq/files/cnvkit/checks/{SRA_tumor_id}_{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}"
    input:
        metadata=selected_samples_location_atacseq,
        tumor_bam="resources/atacseq/files/bam/{SRA_tumor_id}/{SRA_tumor_id}_{assembly_code}_{mapq_thr}.fixed.bam",
        backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed",
        genome="resources/reference_data/genome/{assembly_code}.fa"
    log:
        "logs/run_CNVkit/{assembly_code}_{window_mode}_{window_size}_{step}_{SRA_tumor_id}_{mapq_thr}.log"
    params:
        num_cores_atacseq=num_cores_atacseq,
        path_dir="resources/atacseq/files/cnvkit/",
        normal_paths_atacseq=" ".join(normal_paths_atacseq)
    conda:
        "../../envs/cnvkit_env.yml"
    shell:
        """
        touch {output}

        cnvkit.py batch {input.tumor_bam} -n {params.normal_paths_atacseq} \
        --fasta {input.genome} \
        --method wgs \
        --target-avg-size {wildcards.window_size} \
        --targets {input.backbone} \
        --processes {params.num_cores_atacseq} \
        -d {params.path_dir}/{wildcards.SRA_tumor_id}
        """
