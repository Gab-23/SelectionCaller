rule aggregate_and_visualize_models:
    output:
        performance="",
        selfeats="",
        hyperparam=""
    input:
        check="",
        performance_dir="",
        selfeats_dir="",
        hyperparam_dir=""
    log:
        "logs/aggregate_models/.log"
    shell:
        """
        Rscript workflow/scripts/aggregate_models.R \
            --input_dir_performance={input.performance_dir} \
            --input_dir_selfeats={input.selfeats_dir} \
            --input_dir_hyperparam={input.hyperparam_dir} \
            --output_path_performance={output.performance} \
            --output_path_selfeats={output.selfeats} \
            --output_path_hyperparam={output.hyperparam} 2> {log}
        """