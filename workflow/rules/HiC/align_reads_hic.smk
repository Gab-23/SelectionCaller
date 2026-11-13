rule align_reads_hic:
    # NOTE (02/11/2026): (1) threads dichiarati (= num_cores_hic, gli stessi passati a bwa): prima la regola risultava a 1 core,
    # quindi Snakemake poteva lanciare piu' allineamenti insieme sullo stesso nodo, ognuno con num_cores_hic thread e ~6 GB
    # per l'indice di hg38 -> memoria esaurita e bwa terminato a meta' (intestazione del BAM troncata, vista sul polmone).
    # (2) memoria dichiarata (mem_mb), per lo stesso motivo. (3) gli errori di bwa e samblaster finiscono nel log (prima solo
    # quelli di samtools). Il comando di allineamento e' INVARIATO.
    output:
        "resources/hic/files/bam/{SRA_id}/{assembly_code}.bam"
    input:
        fasta="resources/reference_data/genome/{assembly_code}.fa",
        index="resources/reference_data/genome/{assembly_code}.fa.bwt",
        FWD="resources/hic/files/fastq/{SRA_id}/{SRA_id}_1.fastq.gz",
        REV="resources/hic/files/fastq/{SRA_id}/{SRA_id}_2.fastq.gz"
    log:
        "logs/align_reads_hic/{assembly_code}_{SRA_id}.log"
    params:
        num_cores_hic=num_cores_hic
    threads:
        num_cores_hic
    resources:
        mem_mb=16000
    conda:
        "../../envs/hic_env.yml"
    shell:
        """
        bwa mem -t {params.num_cores_hic} -5SP {input.fasta} {input.FWD} {input.REV} 2>{log}.bwa | \
        samblaster 2>{log}.samblaster | \
        samtools view -S -h -b -F 2316 > {output} 2>{log}
        """