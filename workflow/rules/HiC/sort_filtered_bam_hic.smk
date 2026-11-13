rule sort_filtered_bam_hic:
    output:
        "resources/hic/files/bam/{SRA_id}/{assembly_code}.canon.sorted.bam"
    input:
        "resources/hic/files/bam/{SRA_id}/{assembly_code}.canon.bam"
    log:
        "logs/sort_filtered_bam_hic/{assembly_code}_{SRA_id}.log"
    params:
        num_cores_hic=num_cores_hic
    threads:
        num_cores_hic
    resources:
        mem_mb=24000
    conda:
        "../../envs/hic_env.yml"
    shell:
        """
        samtools sort -@ {params.num_cores_hic} -n {input} -o {output} 2>{log}
        """