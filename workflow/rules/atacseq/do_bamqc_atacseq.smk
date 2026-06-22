rule do_bamqc_atacseq:
    output:
        directory("resources/atacseq/files/bamqc/{assembly_code}_{SRA_id}/")
    input:
        "resources/atacseq/files/bam/{SRA_id}/{assembly_code}.canon.sorted.bam"
    log:
        "logs/do_bamqc_atacseq/{assembly_code}_{SRA_id}.log"
    params:
        java_mem_size_atacseq=java_mem_size_atacseq
    conda:
        "../../envs/bamqc.yml"
    shell:
        """
        qualimap bamqc \
        --java-mem-size={params.java_mem_size_atacseq} \
        -bam {input} \
        -outdir {output} \
        -outformat HTML \
        -c 2>{log}
        """


