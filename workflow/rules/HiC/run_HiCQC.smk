rule run_HiCQC:
    output:
        directory("resources/hic/files/qc/hicqc/{assembly_code}_{SRA_id}/")
    input:
        "resources/hic/files/bam/{SRA_id}/{assembly_code}.canon.sorted.bam"
    log:
        "logs/run_HiCQC/{assembly_code}_{SRA_id}.log"
    params:
        hic_qc_script=hic_qc_script_location 
    shell:
        """
        source ~/miniconda3/etc/profile.d/conda.sh
        conda activate hic_qc

        python {params.hic_qc_script} -b {input} -o {output} 2>{log}
        """
