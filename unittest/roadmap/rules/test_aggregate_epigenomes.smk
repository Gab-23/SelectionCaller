rule test_aggregate_epigenomes:
    output:
	    "unittest/roadmap/output/aggregate_epigenomes_dummy_result.tsv"
    input:
        "unittest/roadmap/output/process_epigenomes_dummy_result"
    conda:
        "../envs/r_roadmap.yml"
    shell:
        """
        Rscript workflow/scripts/aggregate_epigenomes.R \
        --input_path={input} \
        --output_path={output}
        """