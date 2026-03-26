rule do_bamqc_repliseq:
    output:
        directory("resources/repliseq/files/bamqc/{assembly_code}_{SRA_id}/")
    input:
        "resources/repliseq/files/bam/{SRA_id}/{assembly_code}.canon.sorted.bam"
    log:
        "logs/do_bamqc_repliseq/{assembly_code}_{SRA_id}.log"
    params:
        java_mem_size_repliseq=java_mem_size_repliseq
    conda:
        "../envs/bamqc.yml"
    shell:
        """
        qualimap bamqc \
        --java-mem-size={params.java_mem_size_repliseq} \
        -bam {input} \
        -outdir {output} \
        -outformat HTML \
        -c 2>{log}
        """


