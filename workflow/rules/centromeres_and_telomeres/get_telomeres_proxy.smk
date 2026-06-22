rule get_telomeres_proxy:
    output:
        "resources/reference_data/telomeres/{assembly_code}.bed"
    input:
        "resources/reference_data/chrom_sizes/{assembly_code}.clean.chrom.sizes"
    conda:
        "../../envs/r_roadmap.yml"
    log:
        "logs/get_telomeres_proxy/{assembly_code}.log"
    shell:
        """
        Rscript workflow/scripts/get_telomeres_proxy.R \
            --input_path={input} \
            --output_path={output} 2>{log}
        """