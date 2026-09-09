rule process_CNVkit_output_atacseq:
    output:
        "resources/atacseq/files/cnvkit/{SRA_id}/{SRA_id}_{assembly_code}_{mapq_thr}.processed"
    input:
        "resources/atacseq/files/cnvkit/{SRA_id}/{SRA_id}_{assembly_code}_{mapq_thr}.fixed.cns"
    log:
        "logs/process_CNVkit_output/{SRA_id}_{assembly_code}_{mapq_thr}.log"
    shell:
        """
        cat {input} | \
        awk -v OFS='\t' '{{print $1,$2,$3,$5,2^$5,$6,$7,$8}}' |
        tail -n+2 | \
        awk -v OFS='\t' '{{print $1,$2,$3,$4}}' \
        > {output}
        """