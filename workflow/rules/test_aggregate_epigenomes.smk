rule test_aggregate_epigenomes:
    output:
	    "unittest/roadmap/input/aggregate_epigenomes_observed_output.tsv"
    input:
        "unittest/roadmap/input/process_epigenomes_observed_output"
    conda:
        "../envs/r_roadmap.yml"
    shell:
        """
        Rscript workflow/scripts/aggregate_epigenomes.R \
        --input_path={input} \
        --output_path={output}
        """