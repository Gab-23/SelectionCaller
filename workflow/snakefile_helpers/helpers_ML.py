import os
import pandas as pd
from snakemake.exceptions import WorkflowError


def check_dataset_exists(dataset_paths):
    """Fail fast, with a clear and actionable message, if the feature matrix this pipeline
    depends on hasn't been generated yet. snakefile_ML is intentionally decoupled from
    snakefile_dataset_bin_level (no include: of combine_dataset.smk), so Snakemake's own DAG-building
    has no way to know the real reason a dataset path is missing - without this check, a
    missing dataset would surface as a generic MissingInputException deep into DAG resolution,
    which is far less helpful than failing immediately, before anything else runs, with a
    message that tells the user exactly what to do.
    """
    missing = [p for p in dataset_paths if not os.path.isfile(p)]
    if missing:
        missing_list = "\n  - ".join(missing)
        raise WorkflowError(
            "snakefile_ML: the following configured dataset(s) do not exist on disk:\n"
            f"  - {missing_list}\n"
            "These are expected to be produced by snakefile_dataset_bin_level's combine_dataset rule. "
            "Run that pipeline first (or fix input_path_datasets in config/config_ML.yaml if "
            "the path is simply wrong) before running snakefile_ML."
        )


def read_model_registry(path):
    """Read the MODELS file (model / handles_NAs / needs_scaling / is_tree_based / is_linear)
    and return it as a DataFrame indexed by model name."""
    if not os.path.isfile(path):
        raise WorkflowError(
            f"snakefile_ML: input_path_models ({path}) does not exist. This file lists the "
            "candidate regressors and their properties (handles_NAs, needs_scaling, "
            "is_tree_based, is_linear) - see MODELS in the repository root for the expected format."
        )
    return pd.read_csv(path, header=0, sep="\t")


def read_blockings(path):
    """Read the BLOCKINGS file (one blocking value per line, e.g. 5/10/15/20, optionally
    "noblocking") and return it as a list of strings."""
    if not os.path.isfile(path):
        raise WorkflowError(
            f"snakefile_ML: input_path_blockings ({path}) does not exist. This file lists the "
            "spatial-blocking sizes to run (one per line); add a literal \"noblocking\" line "
            "yourself if you also want a plain (non-spatial) train/test split - it is not "
            "added automatically."
        )
    with open(path) as f:
        return [line.rstrip() for line in f.readlines() if line.rstrip() != ""]


def resolve_selection(selected, available, kind):
    """Resolve a selected_models/selected_blockings config entry (["all"], or an explicit
    list) against what's actually available in the corresponding registry file, failing
    clearly on anything that doesn't match rather than silently ignoring it."""
    if list(selected) == ["all"]:
        return list(available)

    unknown = [x for x in selected if x not in available]
    if unknown:
        raise WorkflowError(
            f"snakefile_ML: the following selected {kind} are not present in the registry "
            f"file: {unknown}. Available {kind}: {list(available)}."
        )
    return list(selected)
