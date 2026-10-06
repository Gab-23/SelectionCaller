rule compute_GTEx_features:
    # Per-bin median tissue-relative expression of the bin's genes, complex partners, PPIs
    # and all interactors in trans, in the GTEx cohort matched to the tumour type.
    output:
        "data/GTEx/features/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{GTEx_cohort}_{GTEx_release}"
    input:
        gene_sets="data/protein_interactions/gene_sets/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.tsv",
        expression="data/GTEx/{GTEx_release}/tissue_relative_expression.tsv",
        gencode="resources/reference_data/gencode/PC/v{gencode_version}_gene.sorted.bed"
    wildcard_constraints:
        GTEx_release="v[0-9]+",
        GTEx_cohort="[A-Z]+"
    log:
        "logs/compute_GTEx_features/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{GTEx_cohort}_{GTEx_release}.log"
    conda:
        "../../envs/r_plain.yml"
    shell:
        """
        Rscript workflow/scripts/compute_GTEx_features.R \
        --input_path_gene_sets={input.gene_sets} \
        --input_path_expression={input.expression} \
        --input_path_gencode={input.gencode} \
        --GTEx_cohort={wildcards.GTEx_cohort} \
        --output_path={output} 2> {log}
        """