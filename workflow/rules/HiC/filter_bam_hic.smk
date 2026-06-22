rule filter_bam_hic:
    output:
        "resources/hic/files/bam/{SRA_id}/{assembly_code}.canon.bam"
    input:
        bam="resources/hic/files/bam/{SRA_id}/{assembly_code}.sorted.bam",
        index="resources/hic/files/bam/{SRA_id}/{assembly_code}.sorted.bam.bai"
    log:
        "logs/filter_bam_hic/{assembly_code}_{SRA_id}.log"
    params:
        num_cores_hic=num_cores_hic
    conda:
        "../../envs/hic_env.yml"
    shell:
        """
        samtools view -@ {params.num_cores_hic} -o {output} {input} $(echo chr{{1..22}} chrX chrY) 2>{log}
        """

