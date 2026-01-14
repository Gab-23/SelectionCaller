options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path             <- args$input_path
output_path            <- args$output_path

gencode_gr <- rtracklayer::import.gff3(input_path)
GenomicRanges::mcols(gencode_gr) <- GenomicRanges::mcols(gencode_gr)[c("ID", "gene_id", "gene_name", "transcript_support_level", "tag", "type", "transcript_type")]

gencode_gr$tag <- sapply(gencode_gr$tag, function(x) paste(x, collapse = "; "))
gencode_gr$ID <- gsub("\\.[0-9]*", "", gencode_gr$ID)
names(gencode_gr) <- gencode_gr$ID

pc_transcripts_gr <- gencode_gr[gencode_gr$type == "transcript" & gencode_gr$transcript_type == "protein_coding" & GenomicRanges::seqnames(gencode_gr) != "chrM"]

output_dataset <- cbind(as.data.frame(pc_transcripts_gr@seqnames),                # pc_transcripts_tss_gr
                        as.data.frame(pc_transcripts_gr@ranges),
                        as.data.frame(pc_transcripts_gr@strand))[,c(1,2,3,5,6)]   # [,c(1:3,5)] here we need to enforce end = start  + 1

colnames(output_dataset) <- c("chr", "start", "end", "names", "strand")
output_dataset$start <- ifelse(output_dataset$strand == "-", (output_dataset$end - 1), output_dataset$start)
output_dataset$end <- output_dataset$start + 1
output_dataset$strand <- NULL

write.table(output_dataset,
            file=output_path,
            row.names=F,
            col.names=F,
            quote=F,
            sep="\t")