rule aggregate_hic_eigen:
    # Ricompone i .eigen per-cromosoma in un unico file genome-wide, nello STESSO path che
    # combine_dataset.smk si aspettava gia' (nessuna modifica a combine_dataset.smk/.R).
    output:
        "data/hic/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}.eigen"
    input:
        lambda wildcards: expand(
            "data/hic/{{assembly_code}}_{{mapq_thr}}_{{window_mode}}_{{window_size}}_{{step}}_{{hic_cell_line}}_{chrom}.eigen",
            chrom=get_chroms_list(wildcards)
        )
    log:
        "logs/aggregate_hic_eigen/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}.log"
    shell:
        """
        awk 'FNR==1 && NR!=1 {{ next }} {{ print }}' {input} > {output} 2> {log}
        """