options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                                  <- args$input_path
output_path_protein_coding_5_UTR            <- args$output_path_protein_coding_5_UTR
output_path_protein_coding_3_UTR            <- args$output_path_protein_coding_3_UTR
output_path_protein_coding_gene             <- args$output_path_protein_coding_gene
output_path_protein_coding_exon             <- args$output_path_protein_coding_exon
output_path_protein_coding_transcript       <- args$output_path_protein_coding_transcript
output_path_lncRNA_gene                     <- args$output_path_lncRNA_gene
output_path_lncRNA_exon                     <- args$output_path_lncRNA_exon
output_path_miRNA_gene                      <- args$output_path_miRNA_gene
output_path_miRNA_exon                      <- args$output_path_miRNA_exon
output_path_rRNA_gene                       <- args$output_path_rRNA_gene
output_path_rRNA_exon                       <- args$output_path_rRNA_exon

gtf <- rtracklayer::import.gff3(input_path)

protein_coding_5_UTR <- as.data.frame(gtf[(GenomicRanges::mcols(gtf)$gene_type == "protein_coding") & (GenomicRanges::mcols(gtf)$type == "five_prime_UTR")])
protein_coding_3_UTR <- as.data.frame(gtf[(GenomicRanges::mcols(gtf)$gene_type == "protein_coding") & (GenomicRanges::mcols(gtf)$type == "three_prime_UTR")])
protein_coding_gene <- as.data.frame(gtf[(GenomicRanges::mcols(gtf)$gene_type == "protein_coding") & (GenomicRanges::mcols(gtf)$type == "gene")])
protein_coding_exon <- as.data.frame(gtf[(GenomicRanges::mcols(gtf)$gene_type == "protein_coding") & (GenomicRanges::mcols(gtf)$type == "exon")])
protein_coding_transcript <- as.data.frame(gtf[(GenomicRanges::mcols(gtf)$transcript_type == "protein_coding") & (GenomicRanges::mcols(gtf)$type == "transcript")])
protein_coding_transcript$ID <- gsub("\\.[0-9]*", "", protein_coding_transcript$ID)

lncRNA_gene <- as.data.frame(gtf[(GenomicRanges::mcols(gtf)$gene_type == "lncRNA") & (GenomicRanges::mcols(gtf)$type == "gene")])
lncRNA_exon <- as.data.frame(gtf[(GenomicRanges::mcols(gtf)$gene_type == "lncRNA") & (GenomicRanges::mcols(gtf)$type == "exon")])

miRNA_gene <- as.data.frame(gtf[(GenomicRanges::mcols(gtf)$gene_type == "miRNA") & (GenomicRanges::mcols(gtf)$type == "gene")])
miRNA_exon <- as.data.frame(gtf[(GenomicRanges::mcols(gtf)$gene_type == "miRNA") & (GenomicRanges::mcols(gtf)$type == "exon")])

rRNA_gene <- as.data.frame(gtf[(GenomicRanges::mcols(gtf)$gene_type == "rRNA") & (GenomicRanges::mcols(gtf)$type == "gene")])
rRNA_exon <- as.data.frame(gtf[(GenomicRanges::mcols(gtf)$gene_type == "rRNA") & (GenomicRanges::mcols(gtf)$type == "exon")])

output_paths <- list(output_path_protein_coding_5_UTR,
                     output_path_protein_coding_3_UTR,
                     output_path_protein_coding_gene,
                     output_path_protein_coding_exon,
                     output_path_protein_coding_transcript,
                     output_path_lncRNA_gene,
                     output_path_lncRNA_exon,
                     output_path_miRNA_gene,
                     output_path_miRNA_exon,
                     output_path_rRNA_gene,
                     output_path_rRNA_exon)

sub_gtfs <- list(protein_coding_5_UTR,
                 protein_coding_3_UTR,
                 protein_coding_gene,
                 protein_coding_exon,
                 protein_coding_transcript,
                 lncRNA_gene,
                 lncRNA_exon,
                 miRNA_gene,
                 miRNA_exon,
                 rRNA_gene,
                 rRNA_exon)

for (idx in seq(length(output_paths))){

    output_path <- output_paths[[idx]]
    sub_gtf <- sub_gtfs[[idx]]
    sub_gtf$start <- sub_gtf$start - 1

    write.table(x = sub_gtf[,c("seqnames","start","end", "ID", "gene_name", "gene_type","type")],                                                 # write output
            file = output_path, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = F)

}

















