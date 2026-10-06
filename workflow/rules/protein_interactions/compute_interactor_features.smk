rule compute_interactor_features:
    # Per-bin complex density and number of complex partners / PPIs in trans, plus the
    # per-bin gene sets (genes, partners_trans, ppis_trans, all_int_trans) used by the GTEx
    # features. Depends on the backbone only through its gene annotation.
    output:
        features="data/protein_interactions/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        gene_sets="data/protein_interactions/gene_sets/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.tsv"
    input:
        backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.mask.filtered.sorted.annotated.bed",
        interactors="data/protein_interactions/gene_interactors.tsv"
    log:
        "logs/compute_interactor_features/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.log"
    conda:
        "../../envs/r_plain.yml"
    shell:
        """
        Rscript workflow/scripts/compute_interactor_features.R \
        --input_path_backbone={input.backbone} \
        --input_path_interactors={input.interactors} \
        --output_path_features={output.features} \
        --output_path_gene_sets={output.gene_sets} 2> {log}
        """