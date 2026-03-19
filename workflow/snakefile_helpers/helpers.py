import pandas 

def get_SRA_ids(selected_samples_location):
    samples = pandas.read_csv(selected_samples_location, sep = "\t", header = 0)
    SRA_ids = samples["SRA_id"].tolist()
    return SRA_ids

def get_conditions(selected_samples_location):
    samples = pandas.read_csv(selected_samples_location, sep = "\t", header = 0)
    conditions = samples["condition"].unique().tolist()
    return conditions

def get_is_patterned(sequencers_info_location, sequencer_names):
    sequencers = pandas.read_csv(sequencers_info_location, sep = "\t", header = 0)
    is_patterned = all(sequencers[[sequencer_name in sequencer_names for sequencer_name in sequencers.sequencer_name.tolist()]]["is_patterned"])
    return is_patterned

def get_epigenome_ids(wildcards, tissue_type_Roadmap):

    mapping_file = checkpoints.get_mapping_EID_tt.get(**wildcards).output[0]
    mapping = pandas.read_csv(mapping_file, sep = "\t", header = 0)
    
    tissue_type = [tissue_type_Roadmap] if tissue_type_Roadmap != "anyTissue" else mapping.tissue_type.tolist()
    sample_type = config["sample_types"] if config["sample_types"][0] != "anySample" else mapping.sample_type.tolist()
    
    epigenome_ids = mapping[((mapping.tissue_type.isin(tissue_type)) & (mapping.sample_type.isin(sample_type)))]["epigenome_id"].tolist()
    
    if len(epigenome_ids) == 0:
        raise ValueError("Empty DataFrame")
    else:
        return epigenome_ids

def map_tt_acronym(mapping_file_location, tumor_type, cohort_project):
    try:
        mapping = pandas.read_csv(mapping_file_location, header = 0, sep = "\t")
        tissue_type_Roadmap = mapping[mapping.common_name == tumor_type].tissue_type_Roadmap.values[0]
        tumor_type_IntOGen = mapping[mapping.common_name == tumor_type].tumor_type_IntOGen.values[0]
        tumor_type_GDC = mapping[mapping.common_name == tumor_type].tumor_type_GDC.values[0]
        oncotree_code_DepMap = mapping[mapping.common_name == tumor_type].oncotree_code_DepMap.values[0]
        COSMIC_parameters = mapping[mapping.common_name == tumor_type].COSMIC_parameters.values[0]
        if cohort_project == "TCGA":
            cohort_IntOGen = mapping[mapping.common_name == tumor_type].TCGA_cohort_IntOGen.values[0]
            cohort_GDC = mapping[mapping.common_name == tumor_type].cohort_GDC.values[0]
        else:
            pass    
        return [tissue_type_Roadmap, tumor_type_IntOGen, cohort_IntOGen, oncotree_code_DepMap, COSMIC_parameters, tumor_type_GDC, cohort_GDC]
    except FileNotFoundError as err:
        print(f"{mapping_file_location} does not exist")

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