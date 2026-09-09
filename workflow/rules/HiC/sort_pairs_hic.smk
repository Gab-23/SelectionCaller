rule sort_pairs_hic:
    output:
        "resources/hic/files/pairs/{SRA_id}/{assembly_code}.sorted.pairs"
    input:
        "resources/hic/files/pairs/{SRA_id}/{assembly_code}.pairs"
    log:
        "logs/sort_pairs_hic/{assembly_code}_{SRA_id}.log"
    params:
        num_cores_hic=num_cores_hic
    conda:
        "../../envs/hic_processing_env.yml"
    shell:
        """
        pairtools sort --nproc {params.num_cores_hic} -o {output} {input} 2>{log}
        """