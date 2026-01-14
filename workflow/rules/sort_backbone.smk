rule sort_backbone:
    output:
        "resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed"
    input:
        "resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/sort_backbone/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    shell:
        """
	     bedtools sort -i {input} > {output} 2> {log}
        """