rule compute_library_size_repliseq:
    # NOTE (30/10/2026): {repliseq_cell_line} nel nome (uscita, log, controllo in input): le dimensioni di libreria sono quelle
    # dei campioni della linea cellulare (es. HMEC per il seno, T2P per il LUAD). Prima il file era unico per assembly/mapq e
    # una linea sovrascriveva l'altra; e gli allineamenti .noblack.sorted.bam sono temp(), quindi non si poteva ricalcolarlo.
    output:
        "resources/repliseq/files/library_size/{assembly_code}_{mapq_thr}_{repliseq_cell_line}"
    input:
        # NOTE (fix): this used to depend only on the "check" marker file from the
        # all_sort_blacklisted_alignment_repliseq collector rule, while the actual
        # .noblack.sorted.bam files it reads (via a runtime shell glob) were never
        # declared as real Snakemake inputs. Those files are temp() - Snakemake deletes
        # them as soon as their DECLARED consumers are done, which can (and did, in
        # testing with parallel execution) happen before this rule gets a chance to run,
        # causing a race condition where the glob matches nothing. Declaring the actual
        # per-sample bam files as inputs (reusing the same helper function the collector
        # rule itself uses) tells Snakemake to keep them alive until this rule has run too.
        bams=get_sort_noblacklisted_alignment_outputs_repliseq,
        check="resources/repliseq/checks/{assembly_code}_{mapq_thr}_{repliseq_cell_line}_noblacklist_alignment_sorting"
    log:
        "logs/compute_library_size/{assembly_code}_{mapq_thr}_{repliseq_cell_line}.log"
    conda:
        "../../envs/bwa_samtools.yml"
    shell:
        """
        touch {output}
        printf "SRA_id\\tlibrary_size\\n" >> {output}

        for file in {input.bams};
        do
        SAMPLE_NAME=$(cut -d "/" -f 5 <<< "$file")
        LIB_SIZE=$(samtools view -c "$file")
        printf "$SAMPLE_NAME\\t$LIB_SIZE\\n" >> {output}
        done 2>{log}
        """