rule all_train_models:
    input:
        get_train_model_outputs
    output:
        all_done="results/ML/{dataset_id}/checks/{chrom}_target_{target_col}_train_models"
    shell:
        """
        echo {input} >> {output}
        """