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
        # NOTE (fix 26/10/2026): one output directory PER RESOLUTION. Before, every binning wrote CNVkit's intermediate
        # files (*.targetcoverage.cnn, reference, *.cns) into the same per-sample directory with names that do not
        # contain the resolution, so two resolutions running in parallel overwrote each other's files (confirmed:
        # "bins do not match" with coverage files of 141858 vs 28398 rows in the same directory).
        path_dir=lambda wildcards: (f"resources/atacseq/files/cnvkit/{wildcards.SRA_tumor_id}/"
                                    f"{wildcards.assembly_code}_{wildcards.window_mode}_{wildcards.window_size}_{wildcards.step}_{wildcards.mapq_thr}"),
        normal_paths_atacseq=lambda wildcards, input: " ".join(input.normal_bams)
    conda:
        "../../envs/cnvkit_env.yml"
    shell:
        """
        mkdir -p {params.path_dir}
        cnvkit.py batch {input.tumor_bam} -n {params.normal_paths_atacseq} \
        --fasta {input.genome} \
        --method wgs \
        --target-avg-size {wildcards.window_size} \
        --targets {input.backbone} \
        --processes {params.num_cores_atacseq} \
        -d {params.path_dir} 2>{log}

        touch {output}
        """