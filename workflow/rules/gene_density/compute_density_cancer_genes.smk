rule compute_density_cancer_genes:
    output:
        OGs="resources/num_overlaps/density/OG/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{tumor_type_IntOGen}.bed",
        TSGs="resources/num_overlaps/density/TSG/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{tumor_type_IntOGen}.bed"
    input:
        backbone_annotated="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.mask.filtered.sorted.annotated.bed",
        OGs="resources/IntOGen/OG/{assembly_code}_v{gencode_version}_{tumor_type_IntOGen}.sorted.bed",
        TSGs="resources/IntOGen/TSG/{assembly_code}_v{gencode_version}_{tumor_type_IntOGen}.sorted.bed"
    log:
        "logs/compute_density_cancer_genes/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{tumor_type_IntOGen}.log"
    conda:
        "../../envs/bedtools.yml"
    shell:
        """
        bedtools intersect -a {input.backbone_annotated} -b {input.OGs} -wa -wb | cut -f 1,2,3,4,5,10 > {output.OGs} 2> {log}
        bedtools intersect -a {input.backbone_annotated} -b {input.TSGs} -wa -wb | cut -f 1,2,3,4,5,10 > {output.TSGs} 2> {log}
        """