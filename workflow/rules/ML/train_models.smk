rule train_models:
    output:
        check=""
    input:
        path_dataset="",
        path_is_autocorrelated="",
        path_is_multicollinear="",
        path_models="",
        path_blockings=""
    params:
        selected_model=config["selected_models"][0],
        optimizer=config["optimizers"][0],
        n_iter=config["n_iters"][0],
        n_jobs=config["n_jobss"][0]
    log:
        "logs/"
    conda:
        "../../envs/py_ML.yml"
    shell:
        """
        python3 workflow/scripts/train_models.py \
        --input_path_dataset={input.path_dataset} \
        --input_path_is_autocorrelated={input.path_is_autocorrelated} \
        --input_path_is_multicollinear={input.path_is_multicollinear} \
        --input_path_models={input.path_models} \
        --input_path_blockings={input.path_blockings} \
        --selected_model={params.selected_model} \
        --outer_cv={wildcards.outer_cv} \
        --inner_cv={wildcards.inner_cv} \
        --test_size={wildcards.test_size} \
        --optimizer={params.optimizer} \
        --n_iter={params.n_iter} \
        --n_jobs={params.n_jobs} \
        --output_path_warning={output.path_warning} \
        --output_path_hyperparams={output.path_hyperparams} \
        --output_path_selected_features={output.path_selected_features} \
        --output_path_performances={output.path_performances}              
        """