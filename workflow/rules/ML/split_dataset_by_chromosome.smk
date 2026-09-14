rule split_dataset_by_chromosome:
    output:
        "results/ML/{dataset_id}/per_chrom/{chrom}/dataset.tsv"
    input:
        dataset=input_path_dataset
    log:
        "logs/split_dataset_by_chromosome/{dataset_id}_{chrom}.log"
    conda:
        "../../envs/r_roadmap.yml"
    shell:
        """
        Rscript workflow/scripts/split_dataset_by_chromosome.R \
        --input_path={input} \
        --chrom={wildcards.chrom} \
        --output_path={output} 2>{log}
        """