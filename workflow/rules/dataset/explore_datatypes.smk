rule explore_datatypes:
    output:
        path_GC_landscape="results/plots/datatypes/GC_landscape/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}_{tumor_type_IntOGen}_{essentiality_cutoff}_{num_assays_cutoff}.pdf",
        path_GC_distribution="results/plots/datatypes/GC_distribution/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}_{tumor_type_IntOGen}_{essentiality_cutoff}_{num_assays_cutoff}.pdf",
        path_genes_barplot="results/plots/datatypes/genes/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}_{tumor_type_IntOGen}_{essentiality_cutoff}_{num_assays_cutoff}.pdf",
        path_num_mutations_barplot="results/plots/datatypes/num_mutations/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}_{tumor_type_IntOGen}_{essentiality_cutoff}_{num_assays_cutoff}.pdf",
        path_coding_mut_descriptions_barplot="results/plots/datatypes/coding_mut_descriptions/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}_{tumor_type_IntOGen}_{essentiality_cutoff}_{num_assays_cutoff}.pdf"
    input:
        # NOTE: Snakemake requires every output/log of a rule to share the exact same wildcard
        # set, so all wildcards appear in every output filename, even those a given plot does
        # not depend on (e.g. GC_landscape does not depend on essentiality_cutoff).
        path_GC="resources/GC/{assembly_code}_{window_mode}_{window_size}_{step}",
        path_OGs="resources/IntOGen/OG/{assembly_code}_v{gencode_version}_{tumor_type_IntOGen}.sorted.bed",
        path_TSGs="resources/IntOGen/TSG/{assembly_code}_v{gencode_version}_{tumor_type_IntOGen}.sorted.bed",
        path_essential="resources/DepMap/essential/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.bed",
        path_coding_mutations="resources/COSMIC/mutations/coding/{gender_COSMIC}_{COSMIC_parameters}_w_description.sorted.bed",
        path_noncoding_mutations="resources/COSMIC/coordinates/noncoding/{gender_COSMIC}_{COSMIC_parameters}.sorted.bed"
    log:
        "logs/explore_datatypes/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}_{tumor_type_IntOGen}_{essentiality_cutoff}_{num_assays_cutoff}.log"
    conda:
        "../../envs/r_plotting.yml"
    shell:
        """
        Rscript workflow/scripts/explore_datatypes.R \
        --input_path_GC={input.path_GC} \
        --input_path_OGs={input.path_OGs} \
        --input_path_TSGs={input.path_TSGs} \
        --input_path_essential={input.path_essential} \
        --input_path_coding_mutations={input.path_coding_mutations} \
        --input_path_noncoding_mutations={input.path_noncoding_mutations} \
        --output_path_GC_landscape={output.path_GC_landscape} \
        --output_path_GC_distribution={output.path_GC_distribution} \
        --output_path_genes_barplot={output.path_genes_barplot} \
        --output_path_num_mutations_barplot={output.path_num_mutations_barplot} \
        --output_path_coding_mut_descriptions_barplot={output.path_coding_mut_descriptions_barplot} 2>{log}
        """