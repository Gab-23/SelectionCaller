rule get_dist_num_overlaps_genes:
    output:
        num_overlaps_protein_coding="resources/num_overlaps/PC/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.bed",
        dist_protein_coding="resources/distances/PC/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.bed",
        num_overlaps_lncRNA="resources/num_overlaps/lncRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.bed",
        dist_lncRNA="resources/distances/lncRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.bed",
        num_overlaps_miRNA="resources/num_overlaps/miRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.bed",
        dist_miRNA="resources/distances/miRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.bed",
        num_overlaps_rRNA="resources/num_overlaps/rRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.bed",
        dist_rRNA="resources/distances/rRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.bed",
    input:
        path_protein_coding_gene="resources/reference_data/gencode/PC/v{gencode_version}_gene.sorted.bed",
        path_lncRNA_gene="resources/reference_data/gencode/lncRNA/v{gencode_version}_gene.sorted.bed",
        path_miRNA_gene="resources/reference_data/gencode/miRNA/v{gencode_version}_gene.sorted.bed",
        path_rRNA_gene="resources/reference_data/gencode/rRNA/v{gencode_version}_gene.sorted.bed",
        path_backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed"
    log:
        "logs/get_dist_num_overlaps_genes/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.log"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools intersect -a {input.path_backbone} -b {input.path_protein_coding_gene} -c > {output.num_overlaps_protein_coding} 2>> {log}
        bedtools closest -a {input.path_backbone} -b {input.path_protein_coding_gene} -t first -d | cut -f 1,2,3,4,8,9,11,12 > {output.dist_protein_coding} 2>> {log}
        bedtools intersect -a {input.path_backbone} -b {input.path_lncRNA_gene} -c > {output.num_overlaps_lncRNA} 2>> {log}
        bedtools closest -a {input.path_backbone} -b {input.path_lncRNA_gene} -t first -d | cut -f 1,2,3,4,8,9,11,12 > {output.dist_lncRNA} 2>> {log}
        bedtools intersect -a {input.path_backbone} -b {input.path_miRNA_gene} -c > {output.num_overlaps_miRNA} 2>> {log}
        bedtools closest -a {input.path_backbone} -b {input.path_miRNA_gene} -t first -d | cut -f 1,2,3,4,8,9,11,12 > {output.dist_miRNA} 2>> {log}
        bedtools intersect -a {input.path_backbone} -b {input.path_rRNA_gene} -c > {output.num_overlaps_rRNA} 2>> {log}
        bedtools closest -a {input.path_backbone} -b {input.path_rRNA_gene} -t first -d | cut -f 1,2,3,4,8,9,11,12 > {output.dist_rRNA} 2>> {log}
        """