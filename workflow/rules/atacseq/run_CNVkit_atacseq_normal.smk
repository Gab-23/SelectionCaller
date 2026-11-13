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
        # NOTE (fix 26/10/2026): one output directory PER RESOLUTION (see run_CNVkit_atacseq.smk).
        path_dir=lambda wildcards: (f"resources/atacseq/files/cnvkit/{wildcards.SRA_normal_id}/"
                                    f"{wildcards.assembly_code}_{wildcards.window_mode}_{wildcards.window_size}_{wildcards.step}_{wildcards.mapq_thr}"),
        loo_paths=lambda wildcards, input: " ".join(input.loo_bams)
    conda:
        "../../envs/cnvkit_env.yml"
    shell:
        """
        mkdir -p {params.path_dir}
        cnvkit.py batch {input.normal_bam} -n {params.loo_paths} \
        --fasta {input.genome} \
        --method wgs \
        --target-avg-size {wildcards.window_size} \
        --targets {input.backbone} \
        --processes {params.num_cores_atacseq} \
        -d {params.path_dir} 2>{log}

        touch {output}
        """