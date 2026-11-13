import pandas 

def get_SRA_ids(selected_samples_location):
    samples = pandas.read_csv(selected_samples_location, sep = "\t", header = 0)
    SRA_ids = samples["SRA_id"].tolist()
    return SRA_ids

def get_SRA_condition_ids(selected_samples_location, condition):
    samples = pandas.read_csv(selected_samples_location, sep = "\t", header = 0)
    SRA_ids = samples[samples.condition == condition]["SRA_id"].tolist()
    return SRA_ids

def get_conditions(selected_samples_location):
    samples = pandas.read_csv(selected_samples_location, sep = "\t", header = 0)
    conditions = samples["condition"].unique().tolist()
    return conditions

def get_is_patterned(sequencers_info_location, sequencer_names):
    sequencers = pandas.read_csv(sequencers_info_location, sep = "\t", header = 0)
    is_patterned = all(sequencers[[sequencer_name in sequencer_names for sequencer_name in sequencers.sequencer_name.tolist()]]["is_patterned"])
    return is_patterned

def map_tt_acronym(mapping_file_location, tumor_type, cohort_project):
    try:
        mapping = pandas.read_csv(mapping_file_location, header = 0, sep = "\t")
        tissue_type_Roadmap = mapping[mapping.common_name == tumor_type].tissue_type_Roadmap.values[0]
        tumor_type_IntOGen = mapping[mapping.common_name == tumor_type].tumor_type_IntOGen.values[0]
        tumor_type_GDC = mapping[mapping.common_name == tumor_type].tumor_type_GDC.values[0]
        oncotree_code_DepMap = mapping[mapping.common_name == tumor_type].oncotree_code_DepMap.values[0]
        COSMIC_parameters = mapping[mapping.common_name == tumor_type].COSMIC_parameters.values[0]
        GTEx_cohort = mapping[mapping.common_name == tumor_type].GTEx_cohort.values[0]
        if cohort_project == "TCGA":
            cohort_IntOGen = mapping[mapping.common_name == tumor_type].TCGA_cohort_IntOGen.values[0]
            cohort_GDC = mapping[mapping.common_name == tumor_type].cohort_GDC.values[0]
        else:
            pass    
        return [tissue_type_Roadmap, tumor_type_IntOGen, cohort_IntOGen, oncotree_code_DepMap, COSMIC_parameters, tumor_type_GDC, cohort_GDC, GTEx_cohort]
    except FileNotFoundError as err:
        print(f"{mapping_file_location} does not exist")

def map_cell_line(mapping_file_location, tumor_type, column):
    # NOTE (02/11/2026): valore della colonna `column` della mappatura per il tipo tumorale (linea cellulare / tessuto di
    # un dato di riferimento, es. repliseq_cell_line, hic_cell_line); errore se la colonna manca o il valore e' vuoto.
    mapping = pandas.read_csv(mapping_file_location, header = 0, sep = "\t")
    if column not in mapping.columns:
        raise ValueError(f"{mapping_file_location}: manca la colonna {column}")
    values = mapping[mapping.common_name == tumor_type][column].values
    if len(values) != 1 or pandas.isna(values[0]) or str(values[0]).strip() == "":
        raise ValueError(f"{mapping_file_location}: {column} mancante o non unico per {tumor_type}")
    return str(values[0]).strip()

def map_hic_cell_line(mapping_file_location, tumor_type):
    # NOTE (02/11/2026): linea cellulare / tessuto dell'Hi-C (colonna hic_cell_line, es. BRCA -> MCF10A, LUAD -> LUNG)
    return map_cell_line(mapping_file_location, tumor_type, "hic_cell_line")

def map_repliseq_cell_line(mapping_file_location, tumor_type):
    # NOTE (30/10/2026): linea cellulare del Repli-seq per il tipo tumorale (colonna repliseq_cell_line della mappatura,
    # es. BRCA -> MCF10A, LUAD -> T2P). Funzione separata da map_tt_acronym per non cambiarne il valore di ritorno.
    mapping = pandas.read_csv(mapping_file_location, header = 0, sep = "\t")
    if "repliseq_cell_line" not in mapping.columns:
        raise ValueError(f"{mapping_file_location}: manca la colonna repliseq_cell_line")
    values = mapping[mapping.common_name == tumor_type].repliseq_cell_line.values
    if len(values) != 1 or pandas.isna(values[0]) or str(values[0]).strip() == "":
        raise ValueError(f"{mapping_file_location}: repliseq_cell_line mancante o non unico per {tumor_type}")
    return str(values[0]).strip()

def map_gender(mapping_file_location, gender):
    try:
        mapping = pandas.read_csv(mapping_file_location, header = 0, sep = "\t")
        if gender in ["M", "F"]:
            gender_COSMIC = mapping[mapping.common_name == gender].gender_COSMIC.values[0]
            gender_DepMap = mapping[mapping.common_name == gender].gender_DepMap.values[0]
            gender_GDC = mapping[mapping.common_name == gender].gender_GDC.values[0]
        elif gender == "ALL":
            gender_COSMIC = "ALL"
            gender_DepMap = "ALL"
            gender_GDC = "ALL"
        else:
            raise ValueError(f"Gender [{gender}] must be in ['M','F','ALL']")
        return [gender_COSMIC, gender_DepMap, gender_GDC]
    except FileNotFoundError as err:
        print(f"{mapping_file_location} does not exist")

def get_LOO_normal_bam_paths_atacseq(wildcards):
    other_ids = [sid for sid in SRA_normal_ids_atacseq if sid != wildcards.SRA_normal_id]
    return [
        f"resources/atacseq/files/bam/{sid}/{sid}_{wildcards.assembly_code}_{wildcards.mapq_thr}.fixed.bam"
        for sid in other_ids
    ]

def get_CNVkit_check_atacseq(wildcards):
    # NOTE (fix): process_CNVkit_output_atacseq's input .cns file is a real side effect of
    # `cnvkit.py batch` (run by run_CNVkit_atacseq for tumor samples, run_CNVkit_atacseq_normal
    # for normal/NAT samples), but nothing declared it as a Snakemake output anywhere -
    # Snakemake's real dry-run failed with MissingInputException since it had no rule to
    # connect the two. This input function picks whichever of the two "check" marker files
    # actually applies to the given SRA_id, giving Snakemake a real dependency edge to build
    # the DAG on, without having to declare the .cns file itself as a formal Snakemake output
    # (which would need window_mode/window_size/step wildcards that CNVkit's own output
    # filename doesn't have room for, conflicting with Snakemake's requirement that every
    # output of a rule share the same wildcard set).
    if wildcards.SRA_id in SRA_tumor_ids_atacseq:
        return f"resources/atacseq/files/cnvkit/checks/{wildcards.SRA_id}_{wildcards.assembly_code}_{wildcards.window_mode}_{wildcards.window_size}_{wildcards.step}_{wildcards.mapq_thr}"
    elif wildcards.SRA_id in SRA_normal_ids_atacseq:
        return f"resources/atacseq/files/cnvkit/checks/{wildcards.SRA_id}_{wildcards.assembly_code}_{wildcards.window_mode}_{wildcards.window_size}_{wildcards.step}_{wildcards.mapq_thr}_normal"
    else:
        raise ValueError(f"SRA_id {wildcards.SRA_id} is neither a known tumor nor a known NAT/normal sample")