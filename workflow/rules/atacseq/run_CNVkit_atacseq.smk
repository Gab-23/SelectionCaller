rule run_CNVkit_atacseq:
    output:
        "resources/atacseq/files/cnvkit/checks/{SRA_tumor_id}_{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}"
    input:
        metadata=selected_samples_location_atacseq,
        tumor_bam="resources/atacseq/files/bam/{SRA_tumor_id}/{SRA_tumor_id}_{assembly_code}_{mapq_thr}.fixed.bam",
        # NOTE (fix): the normal-panel bam paths used to only be a `params` (a plain string
        # built once at parse time), never a real Snakemake `input:` - so nothing guaranteed
        # these files existed before this rule ran. With parallel execution (--cores > 1) this
        # is a genuine race condition: confirmed with a real run where CNVkit failed with
        # "FileNotFoundError" reading a normal sample's bam that hadn't finished being written
        # yet by its own (still-running) rule. Declaring them as a real input tells Snakemake
        # to build/wait for all of them first (same class of fix as
        # compute_library_size_repliseq's race condition, Phase 10).
        normal_bams=lambda wildcards: [
            f"resources/atacseq/files/bam/{sid}/{sid}_{wildcards.assembly_code}_{wildcards.mapq_thr}.fixed.bam"
            for sid in SRA_normal_ids_atacseq
        ],
        backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed",
        genome="resources/reference_data/genome/{assembly_code}.fa"
    log:
        "logs/run_CNVkit/{assembly_code}_{window_mode}_{window_size}_{step}_{SRA_tumor_id}_{mapq_thr}.log"
    params:
        num_cores_atacseq=num_cores_atacseq,
        path_dir="resources/atacseq/files/cnvkit/",
        normal_paths_atacseq=lambda wildcards, input: " ".join(input.normal_bams)
    conda:
        "../../envs/cnvkit_env.yml"
    shell:
        """
        cnvkit.py batch {input.tumor_bam} -n {params.normal_paths_atacseq} \
        --fasta {input.genome} \
        --method wgs \
        --target-avg-size {wildcards.window_size} \
        --targets {input.backbone} \
        --processes {params.num_cores_atacseq} \
        -d {params.path_dir}/{wildcards.SRA_tumor_id} 2>{log}

        touch {output}
        """
