checkpoint clean_chrom_sizes:
    # Diventato un checkpoint (da rule semplice): la lista dei cromosomi che ne deriva serve
    # a costruire il wildcard {chrom} per la catena Hi-C cis-only (vedi normalize_matrix_hic,
    # process_bins_hic, process_hic, aggregate_hic_eigen) - in una run pulita questo file non
    # esiste ancora al parsing dello Snakefile, quindi la lista va letta DOPO che questa regola
    # ha girato, non prima. Snakemake ricalcola il DAG a valle di un checkpoint una volta che
    # è stato eseguito.
    output:
        "resources/reference_data/chrom_sizes/{assembly_code}.clean.chrom.sizes"
    input:
        "resources/reference_data/chrom_sizes/{assembly_code}.chrom.sizes"
    log:
        "logs/clean_chrom_sizes/{assembly_code}.log"
    shell:
        """
        cat {input} | grep -vE '_|X|M|Y' > {output} 2> {log}
        """