rule deduplicate_bam_repliseq:
    output:
        bam="resources/repliseq/files/bam/{SRA_id}/{assembly_code}.dedup.bam",
        metrics="resources/repliseq/files/bam/{SRA_id}/{assembly_code}.metrics.txt"
    input:
        "resources/repliseq/files/bam/{SRA_id}/{assembly_code}.canon.sorted.bam"
    log:
        "logs/deduplicate_bam_repliseq/{assembly_code}_{SRA_id}.log"
    params:
        pixel_distance_repliseq=pixel_distance_repliseq
    conda:
        "../envs/picard.yml"
    shell:
        """
        picard MarkDuplicates \
            I={input} \
            O={output.bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=true \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE={params.pixel_distance_repliseq} \
            ASSUME_SORTED=true
        """



