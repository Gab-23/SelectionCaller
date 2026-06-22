rule index_reference_genome_bowtie2:
    output:
        "resources/reference_data/genome/{assembly_code}.1.bt2"
    input:
        "resources/reference_data/genome/{assembly_code}.fa"
    log:
        "logs/index_reference_genome_bowtie2/{assembly_code}.log"
    conda:
        "../../envs/bowtie2_samtools.yml"
    shell:
        """
        bowtie2-build {input} resources/reference_data/genome/{wildcards.assembly_code} 2>{log}
        """