rule merge_coolers:
    # NOTE (02/11/2026): {hic_cell_line} nel nome (uscita, controllo, log): una matrice unita per linea cellulare / tessuto.
    # NOTE (02/11/2026, fix): si uniscono SOLO le matrici delle corse della linea (file dei campioni della linea, tramite
    # get_matrices_outputs_hic), dichiarate come input. Prima un glob (cooler/*/<parametri>.cool) prendeva le matrici di
    # TUTTE le corse presenti su disco: con piu' linee (es. MCF10A e polmone) le avrebbe unite insieme.
    # Le corse della linea (repliche tecniche e donatori) vengono sommate, come prima.
    output:
        "resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}.cool"
    input:
        coolers=get_matrices_outputs_hic,
        check="resources/hic/checks/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}_matrices"
    log:
        "logs/merge_coolers/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}.log"
    conda:
        "../../envs/hic_processing_env.yml"
    shell:
        """
        echo "matrici unite ({wildcards.hic_cell_line}): {input.coolers}" >> {log}
        cooler merge {output} {input.coolers} 2>>{log}
        """