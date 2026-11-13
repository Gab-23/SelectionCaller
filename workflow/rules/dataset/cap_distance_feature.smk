# Caps one processed distance table (binID + one distance column) at dist_caps[0] bp: every bin farther
# than the cap from the closest locus gets the cap. bedtools closest -d gives unsigned, same-chromosome
# distances, so the cap is symmetric (upstream = downstream) and chromosome-specific by construction.
# NA (chromosome with no locus) stays NA. Which features go through this rule is decided in
# combine_dataset.smk (capped_distance_features in the config).
# Config-only parameter: the cap is NOT encoded in file names. It is a param, so changing it in the
# config makes Snakemake re-run this rule (rerun-triggers include params by default), but the dataset
# built with the old cap is overwritten, not kept beside the new one.
rule cap_distance_feature:
    output:
        path_capped="data/distances_capped/{dist_feature}/{dist_name}",
        path_plot="results/plots/dataset/distances_capped/{dist_feature}/{dist_name}.pdf"
    input:
        "data/distances/{dist_feature}/{dist_name}"
    params:
        dist_cap=config["dist_caps"][0]
    log:
        "logs/cap_distance_feature/{dist_feature}/{dist_name}.log"
    conda:
        "../../envs/r_plotting.yml"
    shell:
        """
        Rscript workflow/scripts/cap_distance_feature.R \
            --input_path={input} \
            --dist_cap={params.dist_cap} \
            --output_path={output.path_capped} \
            --output_path_plot={output.path_plot} 2>{log}
        """