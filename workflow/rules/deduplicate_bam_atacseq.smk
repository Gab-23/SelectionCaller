rule deduplicate_bam_atacseq:
    output:
        bam=temp("resources/atacseq/files/bam/{SRA_id}/{assembly_code}.dedup.bam"),
        metrics="resources/atacseq/files/bam/{SRA_id}/{assembly_code}.metrics.txt"
    input:
        "resources/atacseq/files/bam/{SRA_id}/{assembly_code}.canon.sorted.bam"
    log:
        "logs/deduplicate_bam_atacseq/{assembly_code}_{SRA_id}.log"
    params:
        pixel_distance_atacseq=pixel_distance_atacseq
    conda:
        "../envs/picard.yml"
    shell:
        """
        picard MarkDuplicates \
            I={input} \
            O={output.bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=true \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE={params.pixel_distance_atacseq} \
            ASSUME_SORTED=true
        """



