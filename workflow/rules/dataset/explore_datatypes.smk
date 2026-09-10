rule explore_datatypes:
    output:
        path_GC_landscape="results/plots/datatypes/GC_landscape/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{tumor_type_GDC}_{gender_GDC}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}_{tumor_type_IntOGen}_{essentiality_cutoff}_{num_assays_cutoff}.pdf",
        path_GC_distribution="results/plots/datatypes/GC_distribution/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{tumor_type_GDC}_{gender_GDC}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}_{tumor_type_IntOGen}_{essentiality_cutoff}_{num_assays_cutoff}.pdf",
        path_samples_barplot="results/plots/datatypes/samples/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{tumor_type_GDC}_{gender_GDC}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}_{tumor_type_IntOGen}_{essentiality_cutoff}_{num_assays_cutoff}.pdf",
        path_genes_barplot="results/plots/datatypes/genes/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{tumor_type_GDC}_{gender_GDC}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}_{tumor_type_IntOGen}_{essentiality_cutoff}_{num_assays_cutoff}.pdf",
        path_num_mutations_barplot="results/plots/datatypes/num_mutations/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{tumor_type_GDC}_{gender_GDC}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}_{tumor_type_IntOGen}_{essentiality_cutoff}_{num_assays_cutoff}.pdf",
        path_coding_mut_descriptions_barplot="results/plots/datatypes/coding_mut_descriptions/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{tumor_type_GDC}_{gender_GDC}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}_{tumor_type_IntOGen}_{essentiality_cutoff}_{num_assays_cutoff}.pdf"
    input:
        # NOTE: this rule file was completely empty in the original codebase (no rule at all),
        # even though workflow/scripts/explore_datatypes.R already existed and is fully
        # functional (verified with toy data covering every input/output below). Built from
        # scratch by reading the script's actual argument usage. Snakemake requires every
        # output/log of a rule to share the exact same wildcard set, so - even though most of
        # these wildcards aren't semantically needed by every individual plot (e.g. GC_landscape
        # doesn't depend on essentiality_cutoff) - all of them appear in every output filename.
        path_GC="resources/GC/{assembly_code}_{window_mode}_{window_size}_{step}",
        path_models_DepMap="resources/DepMap/models/{oncotree_code_DepMap}_{gender_DepMap}",
        path_manifest_GDC="resources/GDC/CNAs/{tumor_type_GDC}_{gender_GDC}_manifest",
        path_samples_COSMIC="resources/COSMIC/samples/{gender_COSMIC}_{COSMIC_parameters}",
        path_OGs="resources/IntOGen/OG/{assembly_code}_v{gencode_version}_{tumor_type_IntOGen}.sorted.bed",
        path_TSGs="resources/IntOGen/TSG/{assembly_code}_v{gencode_version}_{tumor_type_IntOGen}.sorted.bed",
        path_essential="resources/DepMap/essential/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.bed",
        path_coding_mutations="resources/COSMIC/mutations/coding/{gender_COSMIC}_{COSMIC_parameters}_w_description.sorted.bed",
        path_noncoding_mutations="resources/COSMIC/coordinates/noncoding/{gender_COSMIC}_{COSMIC_parameters}.sorted.bed"
    log:
        "logs/explore_datatypes/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{tumor_type_GDC}_{gender_GDC}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}_{tumor_type_IntOGen}_{essentiality_cutoff}_{num_assays_cutoff}.log"
    conda:
        "../../envs/r_plotting.yml"
    shell:
        """
        Rscript workflow/scripts/explore_datatypes.R \
        --input_path_GC={input.path_GC} \
        --input_path_models_DepMap={input.path_models_DepMap} \
        --input_path_manifest_GDC={input.path_manifest_GDC} \
        --input_path_samples_COSMIC={input.path_samples_COSMIC} \
        --input_path_OGs={input.path_OGs} \
        --input_path_TSGs={input.path_TSGs} \
        --input_path_essential={input.path_essential} \
        --input_path_coding_mutations={input.path_coding_mutations} \
        --input_path_noncoding_mutations={input.path_noncoding_mutations} \
        --output_path_GC_landscape={output.path_GC_landscape} \
        --output_path_GC_distribution={output.path_GC_distribution} \
        --output_path_samples_barplot={output.path_samples_barplot} \
        --output_path_genes_barplot={output.path_genes_barplot} \
        --output_path_num_mutations_barplot={output.path_num_mutations_barplot} \
        --output_path_coding_mut_descriptions_barplot={output.path_coding_mut_descriptions_barplot} 2>{log}
        """
