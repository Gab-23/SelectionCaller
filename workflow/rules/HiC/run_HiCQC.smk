rule run_HiCQC:
    output:
        html="resources/hic/files/qc/hicqc/{SRA_id}/{assembly_code}_qc_report.html"
    input:
        "resources/hic/files/bam/{SRA_id}/{assembly_code}.canon.sorted.bam"
    log:
        "logs/run_HiCQC/{assembly_code}_{SRA_id}.log"
    params:
        prefix="resources/hic/files/qc/hicqc/{SRA_id}/{assembly_code}",
        hic_qc_script=hic_qc_script_location
    shell:
        """
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate hic_qc
        python {params.hic_qc_script} -b {input} -o {params.prefix} 2>{log}
        """