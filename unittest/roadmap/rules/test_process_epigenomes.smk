rule test_process_epigenomes:
    output:
	    "unittest/roadmap/output/process_epigenomes_dummy_result"
    input:
        dummy_epigenomes_paths="unittest/roadmap/input/dummy_bed_locations",
        dummy_backbone_path="unittest/roadmap/input/dummy_backbone"
    conda:
        "../envs/r_roadmap.yml"
    shell:
        """
        Rscript workflow/scripts/process_epigenomes.R \
        --epigenomes_paths={input.epigenomes_paths} \
        --backbone_path={input.backbone_path} \
        --output_path={output} \
        --is_test_call=T
        """