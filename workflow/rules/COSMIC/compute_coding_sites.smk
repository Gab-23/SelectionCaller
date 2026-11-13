# Esposizione per classe di mutazione: sostituzioni singole possibili per bin (sinonime, missense, troncanti,
# altre codificanti), dalla sequenza di riferimento, sugli stessi trascritti dei conteggi per classe
# (COSMIC_transcript_mode, parametro solo da config definito in compute_coding_mutations_counts.smk).
# ATTENZIONE: con --rerun-triggers mtime un cambio di COSMIC_transcript_mode NON fa ripartire la regola:
# usare --forcerun compute_coding_sites (oltre a compute_coding_mutations_counts).

rule compute_coding_sites:
    output:
        path_sites="data/num_overlaps/mutations/coding_sites/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}"
    input:
        path_genome="resources/reference_data/genome/{assembly_code}.fa",                                  # da get_reference_genome
        path_gtf="resources/reference_data/gencode/gencode.v{gencode_version}.annotation.gtf.gz",          # da get_gencode_annotation
        path_backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.sorted.bed"
    params:
        transcript_mode=str(TRANSCRIPT_MODE[0])
    conda:
        "../../envs/gc_content.yml"
    resources:
        mem_mb=32000
    log:
        "logs/compute_coding_sites/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.log"
    shell:
        """
        Rscript workflow/scripts/compute_coding_sites.R \
        --input_path_genome={input.path_genome} \
        --input_path_gtf={input.path_gtf} \
        --input_path_backbone={input.path_backbone} \
        --transcript_mode={params.transcript_mode} \
        --output_path={output.path_sites} 2>{log}
        """