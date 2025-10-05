checkpoint assert_roadmap_processing:
    output:
        "unittest/roadmap/output/assert_roadmap_processing_log"
    input:
        observed_output="unittest/roadmap/input/aggregate_epigenomes_observed_output.tsv",
        expected_output="unittest/roadmap/input/aggregate_epigenomes_expected_output.tsv"
    shell:
        """
        diff_output=$(diff {input.observed_output} {input.expected_output} | wc -l)
        echo $diff_output > {output}
        rm unittest/roadmap/input/*processed.tsv
        rm unittest/roadmap/input/process_epigenomes_observed_output
        rm unittest/roadmap/input/aggregate_epigenomes_observed_output.tsv
        """
    