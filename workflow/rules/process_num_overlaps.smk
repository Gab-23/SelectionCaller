rule process_num_overlaps:
    output:
        num_overlaps_protein_coding="data/num_overlaps/PC/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        num_overlaps_lncRNA="data/num_overlaps/lncRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        num_overlaps_miRNA="data/num_overlaps/miRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        num_overlaps_rRNA="data/num_overlaps/rRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}"
    input:
        num_overlaps_protein_coding="resources/num_overlaps/PC/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.bed",
        num_overlaps_lncRNA="resources/num_overlaps/lncRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.bed",
        num_overlaps_miRNA="resources/num_overlaps/miRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.bed",
        num_overlaps_rRNA="resources/num_overlaps/rRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.bed"
    log:
        "logs/process_num_overlaps/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.log"
    conda:
        "../envs/r_roadmap.yml"
    shell:
        """
        Rscript workflow/scripts/process_num_overlaps.R \
            --input_num_overlaps_protein_coding={input.num_overlaps_protein_coding} \
            --input_num_overlaps_lncRNA={input.num_overlaps_lncRNA} \
            --input_num_overlaps_miRNA={input.num_overlaps_miRNA} \
            --input_num_overlaps_rRNA={input.num_overlaps_rRNA} \
            --output_num_overlaps_protein_coding={output.num_overlaps_protein_coding} \
            --output_num_overlaps_lncRNA={output.num_overlaps_lncRNA} \
            --output_num_overlaps_miRNA={output.num_overlaps_miRNA} \
            --output_num_overlaps_rRNA={output.num_overlaps_rRNA}
        """