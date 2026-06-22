rule sort_filtered_bam_hic:
    output:
        "resources/hic/files/bam/{SRA_id}/{assembly_code}.canon.sorted.bam"
    input:
        "resources/hic/files/bam/{SRA_id}/{assembly_code}.canon.bam"
    log:
        "logs/sort_filtered_bam_hic/{assembly_code}_{SRA_id}.log"
    params:
        num_cores_hic=num_cores_hic
    conda:
        "../../envs/hic_env.yaml"
    shell:
        """
        samtools sort -@ {params.num_cores_hic} -n {input} -o {output} 2>{log}
        """