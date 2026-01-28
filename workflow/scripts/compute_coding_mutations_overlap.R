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
                                   mut_ID=mutation$V4,
                                   mut_types=mutation$V5)

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
colnames(overlap_data) <- c("chr", "start", "end", "width", "strand", "binID", "mut_ID", "mut_types")

mut_types <- unique(
                unlist(
                    lapply(X = overlap_data$mut_types, 
                           FUN = function(x){
                            strsplit(x, split = ",")[[1]]
                            }
                        )
                    )
                )

mut_per_category <- matrix(data = 0, nrow = nrow(overlap_data), ncol = length(mut_types))
colnames(mut_per_category) <- mut_types

for (idx in seq(length(overlap_data$mut_types))){
    mut_types_whole <- overlap_data$mut_types[idx]
    mut_types <- strsplit(mut_types_whole, split = ",")[[1]]
    mut_per_category[idx, mut_types] <- mut_per_category[idx, mut_types] + overlap_data[idx, "width"]
}

mut_per_category <- as.data.frame(mut_per_category)
mut_per_category$binID <- overlap_data$binID

overlap_data <- overlap_data[, c("binID", "width")]

length_df_agnostic <- overlap_data %>%                                              
                dplyr::group_by(binID) %>%                                        # for each bin, for each state
                dplyr::mutate(length_coding_agnostic = sum(width)) %>%                   # compute bp as the overlap between a certain state and a bin
                dplyr::ungroup() %>%
                dplyr::select(-width) %>%                
                dplyr::distinct()

length_df_per_category <- mut_per_category %>%                                              
                dplyr::group_by(binID) %>%                                        # for each bin, for each state
                dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(., na.rm = TRUE))) %>%                                           # compute bp as the overlap between a certain state and a bin
                dplyr::ungroup() %>%             
                dplyr::distinct()

colnames(length_df_per_category)[-1] <- paste0("length_coding_", colnames(length_df_per_category)[-1])

output_df <- merge(length_df_agnostic, length_df_per_category, by = "binID")

write.table(x = output_df,                                                        # write output
            file = output_path_mutation, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)