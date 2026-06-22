rule sort_gene_classes:
    output:
        path_protein_coding_5_UTR="resources/reference_data/gencode/PC/v{gencode_version}_5_UTR.sorted.bed",
        path_protein_coding_3_UTR="resources/reference_data/gencode/PC/v{gencode_version}_3_UTR.sorted.bed",
        path_protein_coding_gene="resources/reference_data/gencode/PC/v{gencode_version}_gene.sorted.bed",
        path_protein_coding_exon="resources/reference_data/gencode/PC/v{gencode_version}_exon.sorted.bed",
        path_lncRNA_gene="resources/reference_data/gencode/lncRNA/v{gencode_version}_gene.sorted.bed",
        path_lncRNA_exon="resources/reference_data/gencode/lncRNA/v{gencode_version}_exon.sorted.bed",
        path_miRNA_gene="resources/reference_data/gencode/miRNA/v{gencode_version}_gene.sorted.bed",
        path_miRNA_exon="resources/reference_data/gencode/miRNA/v{gencode_version}_exon.sorted.bed",
        path_rRNA_gene="resources/reference_data/gencode/rRNA/v{gencode_version}_gene.sorted.bed",
        path_rRNA_exon="resources/reference_data/gencode/rRNA/v{gencode_version}_exon.sorted.bed"
    input:
        path_protein_coding_5_UTR="resources/reference_data/gencode/PC/v{gencode_version}_5_UTR.bed",
        path_protein_coding_3_UTR="resources/reference_data/gencode/PC/v{gencode_version}_3_UTR.bed",
        path_protein_coding_gene="resources/reference_data/gencode/PC/v{gencode_version}_gene.bed",
        path_protein_coding_exon="resources/reference_data/gencode/PC/v{gencode_version}_exon.bed",
        path_lncRNA_gene="resources/reference_data/gencode/lncRNA/v{gencode_version}_gene.bed",
        path_lncRNA_exon="resources/reference_data/gencode/lncRNA/v{gencode_version}_exon.bed",
        path_miRNA_gene="resources/reference_data/gencode/miRNA/v{gencode_version}_gene.bed",
        path_miRNA_exon="resources/reference_data/gencode/miRNA/v{gencode_version}_exon.bed",
        path_rRNA_gene="resources/reference_data/gencode/rRNA/v{gencode_version}_gene.bed",
        path_rRNA_exon="resources/reference_data/gencode/rRNA/v{gencode_version}_exon.bed"
    conda:
        "../../envs/bedtools.yml"
    log:
        "logs/sort_gene_classes/v{gencode_version}.log"
    shell:
        """
        bedtools sort -i {input.path_protein_coding_5_UTR} > {output.path_protein_coding_5_UTR} 2>> {log}
        bedtools sort -i {input.path_protein_coding_3_UTR} > {output.path_protein_coding_3_UTR} 2>> {log}
        bedtools sort -i {input.path_protein_coding_gene} > {output.path_protein_coding_gene} 2>> {log}
        bedtools sort -i {input.path_protein_coding_exon} > {output.path_protein_coding_exon} 2>> {log}
        bedtools sort -i {input.path_lncRNA_gene} > {output.path_lncRNA_gene} 2>> {log}
        bedtools sort -i {input.path_lncRNA_exon} > {output.path_lncRNA_exon} 2>> {log}
        bedtools sort -i {input.path_miRNA_gene} > {output.path_miRNA_gene} 2>> {log}
        bedtools sort -i {input.path_miRNA_exon} > {output.path_miRNA_exon} 2>> {log}
        bedtools sort -i {input.path_rRNA_gene} > {output.path_rRNA_gene} 2>> {log}
        bedtools sort -i {input.path_rRNA_exon} > {output.path_rRNA_exon} 2>> {log}
        """