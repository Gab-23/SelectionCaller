options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args


epigenome_path       <- args$epigenome_path
backbone_path        <- args$backbone_path
output_path          <- args$output_path

`%>%` <- magrittr::`%>%`

epigenome <- rtracklayer::import.bed(epigenome_path)
backbone  <- rtracklayer::import.bed(backbone_path)

hits <- GenomicRanges::findOverlaps(query = backbone, subject = epigenome)        # compute overlaps

query <- backbone[S4Vectors::queryHits(hits)]                                     # take overlapping regions
subject <- epigenome[S4Vectors::subjectHits(hits)]

overlap_data <- GenomicRanges::pintersect(query, subject)                         # generate GRanges with data about overlaps
GenomicRanges::mcols(overlap_data) <- NULL

overlap_metadata <- cbind(GenomicRanges::mcols(query),                            # get binIDs and stateIDs
                          GenomicRanges::mcols(subject))

GenomicRanges::mcols(overlap_data) <- overlap_metadata                            # assign metadata

overlap_data <- data.frame(overlap_data)
overlap_data$width <- overlap_data$width - 1                                      # set as proper BED file coordinate (0-based, half-open)

overlap_data <- overlap_data[, c("binID", "stateID", "width")]

length_counts_df <- overlap_data %>%                                              
                     dplyr::group_by(binID, stateID) %>%                          # for each bin, for each state
                     dplyr::mutate(length_counts = sum(width)) %>%                # compute bp as the overlap between a certain state and a bin
                     dplyr::distinct() %>%
                     dplyr::select(-width) %>%
                     tidyr::pivot_wider(names_from = stateID,                     # each row is a bin, each column is a state
                                        values_from = length_counts)

write.table(x = length_counts_df,                                                 # write output
            file = output_path, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = F)

