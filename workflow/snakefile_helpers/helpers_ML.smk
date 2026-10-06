rule all_train_models:
    input:
        get_train_model_outputs
    output:
        all_done="results/ML/{dataset_id}/checks/{training_mode}_{chrom}_target_{target_col}_train_models"
    shell:
        """
        echo {input} >> {output}
        """