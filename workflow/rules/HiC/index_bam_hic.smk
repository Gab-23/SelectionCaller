rule index_bam_hic:
    output:
        temp("resources/hic/files/bam/{SRA_id}/{assembly_code}.sorted.bam.bai")
    input:
        "resources/hic/files/bam/{SRA_id}/{assembly_code}.sorted.bam"
    log:
        "logs/index_bam_hic/.{assembly_code}_{SRA_id}log"
    params:
        num_cores_hic=num_cores_hic
    conda:
        "../../envs/hic_env.yaml"
    shell:
        """
        samtools index -@ {params.num_cores_hic} {input} 2>{log}
        """

