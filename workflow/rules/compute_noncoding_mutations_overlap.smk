rule compute_noncoding_mutations_overlap:
    output:
        "data/num_overlaps/mutations/noncoding/{assembly_code}_{window_mode}_{window_size}_{step}"
    input:
        path_mutation="resources/COSMIC/coordinates/noncoding/coordinates.sorted.bed",
        path_backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed"
    conda:
        "../envs/r_roadmap.yml"
    log:
        "logs/compute_noncoding_mutations_overlap/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    shell:
        """
        Rscript workflow/scripts/compute_noncoding_mutations_overlap.R \
        --input_path_mutation={input.path_mutation} \
        --input_path_backbone={input.path_backbone} \
        --output_path_mutation={output} 2>{log}
        """