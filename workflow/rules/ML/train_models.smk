rule train_models:
    # NOTE: MODEL and BLOCKING are genuine Snakemake wildcards (one job per combination, run
    # in parallel by Snakemake) - train_models.py is always invoked with a specific
    # --selected_model/--selected_blocking (never "all"), letting Snakemake itself handle the
    # "all models x all blockings" fan-out declared in rule_all, instead of looping internally
    # in Python as a single monolithic job.
    output:
        models=[
            f"results/ML/{dataset_id}/models/{{model}}_{{blocking}}/{{model}}_{{blocking}}_{i}.pkl"
            for i in range(OUTER_CV)
        ],
        warning=f"results/ML/{dataset_id}/logs/{{model}}_{{blocking}}_warning.log",
        hyperparams=f"results/ML/{dataset_id}/hyperparams/hyperparams_{{model}}_{{blocking}}.tsv",
        selected_features=f"results/ML/{dataset_id}/selfeats/selfeats_file_{{model}}_{{blocking}}.tsv",
        performances=f"results/ML/{dataset_id}/performances/perf_file_{{model}}_{{blocking}}.tsv"
    input:
        dataset=input_path_dataset,
        is_autocorrelated=input_path_is_autocorrelated,
        is_multicollinear=input_path_is_multicollinear,
        models=input_path_models,
        blockings=input_path_blockings
    params:
        model_dir=f"results/ML/{dataset_id}/models/{{model}}_{{blocking}}",
        test_size=TEST_SIZE,
        outer_cv=OUTER_CV,
        inner_cv=INNER_CV,
        optimizer=OPTIMIZER,
        n_iter=N_ITER,
        n_jobs=N_JOBS_ML,
        model_n_jobs=MODEL_N_JOBS_ML
    log:
        f"logs/train_model/{dataset_id}_{{model}}_{{blocking}}.log"
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
