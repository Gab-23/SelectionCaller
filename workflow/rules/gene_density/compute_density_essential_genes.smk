rule compute_density_essential_genes:
    output:
        "resources/num_overlaps/density/essential/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.bed"
    input:
        backbone_annotated="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.sorted.annotated.bed",
        essential_genes="resources/DepMap/essential/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.sorted.bed",
    log:
        "logs/compute_density_essential_genes/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.log"
    conda:
        "../../envs/bedtools.yml"
    shell:
        """
        bedtools intersect -a {input.backbone_annotated} -b {input.essential_genes} -wa -wb | cut -f 1,2,3,4,5,9 > {output}
        """