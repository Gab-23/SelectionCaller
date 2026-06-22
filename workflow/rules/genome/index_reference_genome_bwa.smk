rule index_reference_genome_bwa:
    output:
        "resources/reference_data/genome/{assembly_code}.fa.bwt"
    input:
        "resources/reference_data/genome/{assembly_code}.fa"
    log:
        "logs/index_reference_genome_bwa/{assembly_code}.log"
    conda:
        "../../envs/bwa_samtools.yml"
    shell:
        """
        bwa index {input} 2>{log}
        """