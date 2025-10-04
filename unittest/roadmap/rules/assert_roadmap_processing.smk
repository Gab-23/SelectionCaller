checkpoint assert_roadmap_processing:
    output:
        "unittest/roadmap/output/assert_roadmap_processing_log"
    input:
        observed_output"unittest/roadmap/output/aggregate_epigenomes_dummy_result.tsv",
        expected_output="unittest/roadmap/output/expected_output"
    shell:
        """
        Rscript assert_roadmap_processing.R \
        --observed_output={input.observed_output} \
        --expected_output={input.expected_output}
        --output_path={output} \ 
        """
    