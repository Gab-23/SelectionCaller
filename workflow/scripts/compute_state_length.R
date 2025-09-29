compute_state_length <- function(epigenome_path, backbone_path, epigenome_output_path){

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
colnames(overlap_data) <- c("chr", "start", "end", "width", "strand", "binID", "stateID")

overlap_data$width <- overlap_data$width - 1                                      # set as proper BED file coordinate (0-based, half-open)

overlap_data <- overlap_data[, c("binID", "stateID", "width")]

length_state_df <- overlap_data %>%                                              
                     dplyr::group_by(binID, stateID) %>%                          # for each bin, for each state
                     dplyr::mutate(length_counts = sum(width)) %>%                # compute bp as the overlap between a certain state and a bin
                     dplyr::ungroup() %>%
                     dplyr::select(-width) %>%                
                     dplyr::distinct() %>%
                     tidyr::pivot_wider(names_from = stateID,                     # each row is a bin, each column is a state
                                        values_from = length_counts,
                                        values_fill = 0)                          # if an overlap between a bin and a state is not found, it's 0

write.table(x = length_state_df,                                                  # write output
            file = epigenome_output_path, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)

}