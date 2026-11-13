# Soglia per i campioni ipermutati: parametro solo da config (NON un wildcard, non entra nel
# nome dei file), come dist_caps. Esattamente un valore: "none" oppure un numero di eventi
# CODIFICANTI per campione (sostituzioni e indel) oltre il quale il campione viene escluso.
# ATTENZIONE: con --rerun-triggers mtime un cambio di soglia NON fa ripartire la regola:
# usare --forcerun compute_coding_mutations_counts.
HYPERMUTATOR_MAX = config.get("COSMIC_hypermutator_max_mutations", ["none"])
if len(HYPERMUTATOR_MAX) != 1:
    raise ValueError("COSMIC_hypermutator_max_mutations is a config-only parameter (not a wildcard): give exactly one value")

# Trascritti per la classe delle mutazioni (colonne mut_n_syn / mis / trunc / other_coding): parametro solo
# da config, come la soglia. "canonical" = Ensembl_canonical del GTF; "all" = classe piu' grave tra tutti.
# Il GTF e' quello di get_gencode_annotation: la versione di GENCODE e' un wildcard ed entra nel nome
# dell'output (come per i file CDS), perche' la classe delle mutazioni dipende dall'annotazione.
# Stessa ATTENZIONE: dopo un cambio usare --forcerun compute_coding_mutations_counts.
TRANSCRIPT_MODE = config.get("COSMIC_transcript_mode", ["canonical"])
if len(TRANSCRIPT_MODE) != 1 or TRANSCRIPT_MODE[0] not in ("canonical", "all"):
    raise ValueError("COSMIC_transcript_mode is a config-only parameter: give exactly one value, 'canonical' or 'all'")

# NOTE (28/10/2026): come si classificano gli EVENTI (parametro solo da config, stessa ATTENZIONE sul --forcerun):
#   "reference" (default): ogni sostituzione singola e' classificata sulla sequenza di riferimento e sui trascritti di
#                          COSMIC_transcript_mode, con lo stesso codice dei siti (helpers/coding_sites_lib.R);
#   "cosmic":              classe dalle righe del file COSMIC (comportamento precedente). Con "canonical" perde gli eventi
#                          dei geni per cui COSMIC non riporta il trascritto canonico Ensembl (8.7% nella coorte BRCA;
#                          tutti quelli di AKT1, CASP8, RUNX1).
EVENT_ANNOTATION = config.get("COSMIC_event_annotation", ["reference"])
if len(EVENT_ANNOTATION) != 1 or EVENT_ANNOTATION[0] not in ("reference", "cosmic"):
    raise ValueError("COSMIC_event_annotation is a config-only parameter: give exactly one value, 'reference' or 'cosmic'")

rule compute_coding_mutations_counts:
    # Conteggi per bin a livello di evento (campione x mutazione), dal file COSMIC filtrato
    # NON collassato: mut_n_events (ricorrenza inclusa), mut_n_samples, mut_n_sites; e, per le sole
    # sostituzioni singole, eventi per classe: mut_n_syn, mut_n_mis, mut_n_trunc, mut_n_other_coding.
    output:
        path_counts="data/num_overlaps/mutations/coding_counts/{assembly_code}_{window_mode}_{window_size}_{step}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}",
        path_samples="results/dataset/COSMIC_mutation_burden/{assembly_code}_{window_mode}_{window_size}_{step}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}.tsv"
    input:
        path_mutations="resources/COSMIC/mutations/coding/{gender_COSMIC}_{COSMIC_parameters}.filtered.bed",
        path_backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.sorted.bed",
        path_gtf="resources/reference_data/gencode/gencode.v{gencode_version}.annotation.gtf.gz",      # da get_gencode_annotation
        path_genome="resources/reference_data/genome/{assembly_code}.fa"                                # da get_reference_genome (riannotazione)
    params:
        hypermutator_max=str(HYPERMUTATOR_MAX[0]),
        transcript_mode=str(TRANSCRIPT_MODE[0]),
        event_annotation=str(EVENT_ANNOTATION[0])
    conda:
        "../../envs/gc_content.yml"          # NOTE (28/10/2026): stesso ambiente di compute_coding_sites (Biostrings); prima r_roadmap.yml
    resources:
        mem_mb=32000                          # genoma + siti di un cromosoma alla volta, come compute_coding_sites
    log:
        "logs/compute_coding_mutation_counts/{assembly_code}_{window_mode}_{window_size}_{step}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}.log"
    shell:
        """
        Rscript workflow/scripts/compute_coding_mutations_counts.R \
        --input_path_mutations={input.path_mutations} \
        --input_path_backbone={input.path_backbone} \
        --input_path_gtf={input.path_gtf} \
        --input_path_genome={input.path_genome} \
        --hypermutator_max={params.hypermutator_max} \
        --transcript_mode={params.transcript_mode} \
        --event_annotation={params.event_annotation} \
        --output_path={output.path_counts} \
        --output_path_samples={output.path_samples} 2>{log}
        """