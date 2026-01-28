options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_mutation       <-  args$input_path_mutation
input_path_backbone       <-  args$input_path_backbone
output_path_mutation      <-  args$output_path_mutation

`%>%` <- magrittr::`%>%`

mutation <- read.table(file = input_path_mutation, sep = "\t", header = F)
mutation.gr <- GenomicRanges::GRanges(seqnames = mutation$V1, 
                                   ranges=IRanges::IRanges(start=mutation$V2, 
                                                           end=mutation$V3),
                                   mut_ID=mutation$V4)

backbone  <- rtracklayer::import.bed(input_path_backbone)

hits <- GenomicRanges::findOverlaps(query = backbone, subject = mutation.gr)      # compute overlaps

query <- backbone[S4Vectors::queryHits(hits)]                                     # take overlapping regions
subject <- mutation.gr[S4Vectors::subjectHits(hits)]

overlap_data <- GenomicRanges::pintersect(query, subject)                         # generate GRanges with data about overlaps
GenomicRanges::mcols(overlap_data) <- NULL

overlap_metadata <- cbind(GenomicRanges::mcols(query),                            # get binIDs and stateIDs
                        GenomicRanges::mcols(subject))

GenomicRanges::mcols(overlap_data) <- overlap_metadata                            # assign metadata

overlap_data <- data.frame(overlap_data)
colnames(overlap_data) <- c("chr", "start", "end", "width", "strand", "binID", "mut_ID")

overlap_data <- overlap_data[, c("binID", "width")]

length_df_agnostic <- overlap_data %>%                                              
                dplyr::group_by(binID) %>%                                        # for each bin, for each state
                dplyr::mutate(length_noncoding_agnostic = sum(width)) %>%                   # compute bp as the overlap between a certain state and a bin
                dplyr::ungroup() %>%
                dplyr::select(-width) %>%                
                dplyr::distinct()

write.table(x = length_df_agnostic,                                               # write output
            file = output_path_mutation, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)