rule train_models:
    # {training_mode} (per_chrom | genome_wide): in genome_wide {chrom} == "all" and the
    # dataset is the full feature matrix (get_training_dataset bypasses the split).
    output:
        models=[
            "results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/models/{model}_{blocking}/{model}_{blocking}_" + str(i) + ".pkl"
            for i in range(OUTER_CV)
        ],
        warning="results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/logs/{model}_{blocking}_warning.log",
        hyperparams="results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/hyperparams/hyperparams_{model}_{blocking}.tsv",
        selected_features="results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/selfeats/selfeats_file_{model}_{blocking}.tsv",
        performances="results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/performances/perf_file_{model}_{blocking}.tsv"
    input:
        dataset=get_training_dataset,
        is_autocorrelated=input_path_is_autocorrelated,
        is_multicollinear=input_path_is_multicollinear,
        models=input_path_models,
        blockings=input_path_blockings
    params:
        model_dir="results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/models/{model}_{blocking}",
        test_size=TEST_SIZE,
        outer_cv=OUTER_CV,
        inner_cv=INNER_CV,
        optimizer=OPTIMIZER,
        n_iter=N_ITER,
        n_jobs=N_JOBS_ML,
        model_n_jobs=MODEL_N_JOBS_ML,
        cols_to_exclude=get_cols_to_exclude
    log:
        "logs/train_models/{dataset_id}_{training_mode}_{chrom}_{target_col}_{model}_{blocking}.log"
    conda:
        "../../envs/py_ml_env.yml"
    shell:
        """
        mkdir -p {params.model_dir}
        # NOTE: output_path_warning is opened in append mode by train_models.py itself - remove
        # any stale copy from a previous run first, so re-runs don't silently accumulate
        # duplicate warning entries on top of old ones.
        rm -f {output.warning}

        python -m workflow.scripts.train_models \
        --input_path_dataset={input.dataset} \
        --input_path_is_autocorrelated={input.is_autocorrelated} \
        --input_path_is_multicollinear={input.is_multicollinear} \
        --input_path_models={input.models} \
        --input_path_blockings={input.blockings} \
        --target_col={wildcards.target_col} \
        --cols_to_exclude {params.cols_to_exclude:q} \
        --selected_model={wildcards.model} \
        --selected_blocking={wildcards.blocking} \
        --outer_cv={params.outer_cv} \
        --inner_cv={params.inner_cv} \
        --test_size={params.test_size} \
        --optimizer={params.optimizer} \
        --n_iter={params.n_iter} \
        --n_jobs={params.n_jobs} \
        --model_n_jobs={params.model_n_jobs} \
        --output_dir_models={params.model_dir} \
        --output_path_warning={output.warning} \
        --output_path_hyperparams={output.hyperparams} \
        --output_path_selected_features={output.selected_features} \
        --output_path_performances={output.performances} 2>{log}
        """