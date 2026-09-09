rule smooth_dataset:
    output:
    input:
    log:
        "logs/smooth_dataset/"
    conda:
        "../../envs/r_spdep.yml"
    shell:
        """
        """