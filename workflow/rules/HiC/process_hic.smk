rule process_hic:
    # Autovettori per cromosoma. Il segno di ogni autovettore e' orientato in modo che la
    # correlazione con il contenuto in GC del bin sia positiva (phasing, fix 20/10/2026):
    # senza questo il segno era arbitrario e diverso tra cromosomi.
    output:
        "data/hic/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}_{chrom}.eigen"
    input:
        matrix="resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}_{chrom}.balanced.norm.mat",
        bins="resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}_{chrom}.balanced.norm.processed.bins",
        gc="data/GC/{assembly_code}_{window_mode}_{window_size}_{step}"
    log:
        "logs/process_hic/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}_{chrom}.log"
    conda:
        "../../envs/hic_processing_env.yml"
    shell:
        """
        python workflow/scripts/process_hic.py \
        --input_path_matrix={input.matrix} \
        --input_path_bins={input.bins} \
        --input_path_gc={input.gc} \
        --phasing=gc \
        --output_path_eigen={output} 2>{log}
        """