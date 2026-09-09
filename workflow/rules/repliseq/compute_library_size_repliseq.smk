rule compute_library_size:
    output:
        "resources/repliseq/files/library_size/{assembly_code}_{mapq_thr}"
    input:
        "resources/repliseq/checks/{assembly_code}_{mapq_thr}_noblacklist_alignment_sorting"
    log:
        "logs/compute_library_size/{assembly_code}_{mapq_thr}.log"
    conda:
        "../../envs/bwa_samtools.yml"
    shell:
        """
        touch {output}
        printf "SRA_id\tlibrary_size\n" >> {output}

        FILES="resources/repliseq/files/bam/*/{assembly_code}_{mapq_thr}.noblack.sorted.bam"

        for file in $FILES;
        do
        SAMPLE_NAME=$(cut $file -d "/" -f 5)
        LIB_SIZE=$(samtools view -c $file)
        printf "$SAMPLE_NAME\t$LIB_SIZE\n" >> {output} 2>> {log}
        done
        """