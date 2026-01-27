rule process_dist_to_genes:
    output:
        dist_protein_coding="data/distances/PC/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        dist_lncRNA="data/distances/lncRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        dist_miRNA="data/distances/miRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        dist_rRNA="data/distances/rRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
    input:
        dist_protein_coding="resources/distances/PC/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.bed",
        dist_lncRNA="resources/distances/lncRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.bed",
        dist_miRNA="resources/distances/miRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.bed",
        dist_rRNA="resources/distances/rRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.bed",
    log:
        "logs/process_dist_to_genes/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.log"
    conda:
        "../envs/r_roadmap.yml"
    shell:
        """
        Rscript workflow/scripts/process_dist_to_genes.R \
            --input_dist_protein_coding={input.dist_protein_coding} \
            --input_dist_lncRNA={input.dist_lncRNA} \
            --input_dist_miRNA={input.dist_miRNA} \
            --input_dist_rRNA={input.dist_rRNA} \
            --output_dist_protein_coding={output.dist_protein_coding} \
            --output_dist_lncRNA={output.dist_lncRNA} \
            --output_dist_miRNA={output.dist_miRNA} \
            --output_dist_rRNA={output.dist_rRNA} 2>{log}
        """