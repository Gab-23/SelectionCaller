rule compute_pairs_hic:
    output:
        pairs="resources/hic/files/pairs/{SRA_id}/{assembly_code}.pairs",
        stats="resources/hic/files/pairs/{SRA_id}/{assembly_code}.stats",
    input:
        chrom_sizes="resources/reference_data/chrom_sizes/{assembly_code}.clean.chrom.sizes",
        bam="resources/hic/files/bam/{SRA_id}/{assembly_code}.canon.sorted.bam"
    log:
        "logs/compute_pairs_hic/{assembly_code}_{SRA_id}.log"
    params:
        num_cores_hic=num_cores_hic
    conda:
        "../envs/hic_processing_env.yaml"
    shell:
        """
        pairtools parse -o {output.pairs} -c {input.chrom_sizes} \
            --drop-sam --drop-seq --output-stats {output.stats} \
            --assembly {wildcards.assembly_code} \
            --flip \
            --add-columns mapq \
            --walks-policy mask \
            --nproc-in {params.num_cores_hic} \
            --nproc-out {params.num_cores_hic} \
            {input.bam} 2>{log}
        """