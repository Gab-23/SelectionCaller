rule split_gtf:
    output:
        path_protein_coding_5_UTR="resources/reference_data/gencode/PC/v{gencode_version}_5_UTR.bed",
        path_protein_coding_3_UTR="resources/reference_data/gencode/PC/v{gencode_version}_3_UTR.bed",
        path_protein_coding_gene="resources/reference_data/gencode/PC/v{gencode_version}_gene.bed",
        path_protein_coding_exon="resources/reference_data/gencode/PC/v{gencode_version}_exon.bed",
        path_protein_coding_transcript="resources/reference_data/gencode/PC/v{gencode_version}_transcript.bed",
        path_lncRNA_gene="resources/reference_data/gencode/lncRNA/v{gencode_version}_gene.bed",
        path_lncRNA_exon="resources/reference_data/gencode/lncRNA/v{gencode_version}_exon.bed",
        path_miRNA_gene="resources/reference_data/gencode/miRNA/v{gencode_version}_gene.bed",
        path_miRNA_exon="resources/reference_data/gencode/miRNA/v{gencode_version}_exon.bed",
        path_rRNA_gene="resources/reference_data/gencode/rRNA/v{gencode_version}_gene.bed",
        path_rRNA_exon="resources/reference_data/gencode/rRNA/v{gencode_version}_exon.bed"
    input:
        "resources/reference_data/gencode/gencode.v{gencode_version}.annotation.gff3.gz"
    conda:
        "../../envs/r_roadmap.yml"
    log:
        "logs/split_gtf/v{gencode_version}.log"
    shell:
        """
	    Rscript workflow/scripts/split_gtf.R \
            --input_path={input} \
            --output_path_protein_coding_5_UTR={output.path_protein_coding_5_UTR} \
            --output_path_protein_coding_3_UTR={output.path_protein_coding_3_UTR} \
            --output_path_protein_coding_gene={output.path_protein_coding_gene} \
            --output_path_protein_coding_exon={output.path_protein_coding_exon} \
            --output_path_protein_coding_transcript={output.path_protein_coding_transcript} \
            --output_path_lncRNA_gene={output.path_lncRNA_gene} \
            --output_path_lncRNA_exon={output.path_lncRNA_exon} \
            --output_path_miRNA_gene={output.path_miRNA_gene} \
            --output_path_miRNA_exon={output.path_miRNA_exon} \
            --output_path_rRNA_gene={output.path_rRNA_gene} \
            --output_path_rRNA_exon={output.path_rRNA_exon} 2>{log}
        """