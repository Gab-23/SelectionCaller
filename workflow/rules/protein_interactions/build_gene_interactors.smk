rule build_gene_interactors:
    # Gene -> partner table from CORUM (complexes grouped by ComplexID) and Interactome
    # INSIDER (homodimers removed). Gene-level: independent of assembly_code and backbone.
    output:
        "data/protein_interactions/gene_interactors.tsv"
    input:
        CORUM=CORUM_location,
        InteractomeINSIDER=InteractomeINSIDER_location
    log:
        "logs/build_gene_interactors/build_gene_interactors.log"
    conda:
        "../../envs/r_plain.yml"
    shell:
        """
        Rscript workflow/scripts/build_gene_interactors.R \
        --input_path_CORUM={input.CORUM} \
        --input_path_InteractomeINSIDER={input.InteractomeINSIDER} \
        --output_path={output} 2> {log}
        """