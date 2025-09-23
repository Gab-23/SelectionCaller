options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args


epigenome_path       <- args$epigenome_path
backbone_path        <- args$backbone_path
output_path          <- args$output_path

`%>%` <- magrittr::`%>%`

epigenome <- read.table(file = epigenome_path,
                        header = F,
                        sep = "\t",
                        col.names = c("chr", "start", "end", "stateID"))


backbone <- read.table(file = backbone_path,
                       header = F,
                       sep = "\t",
                       col.names = c("chr", "start", "end", "binID"))


epigenome <- GenomicRanges::GRanges(epigenome)
backbone  <- GenomicRanges::GRanges(backbone)

hits <- GenomicRanges::findOverlaps(query = backbone, subject = epigenome)

query <- backbone[S4Vectors::queryHits(hits)]
subject <- epigenome[S4Vectors::subjectHits(hits)]

overlap_data <- GenomicRanges::pintersect(query, subject)
GenomicRanges::mcols(overlap_data) <- NULL

overlap_metadata <- cbind(GenomicRanges::mcols(query), 
                          GenomicRanges::mcols(subject))

GenomicRanges::mcols(overlap_data) <- overlap_metadata

overlap_data <- data.frame(overlap_data)

overlap_data$width <- overlap_data$width - 1

overlap_data <- overlap_data[, c("binID", "stateID", "width")]

length_counts_df <- overlap_data %>% 
                     dplyr::group_by(binID, stateID) %>% 
                     dplyr::mutate(length_counts = sum(width)) %>%
                     dplyr::distinct() %>%
                     dplyr::select(-width) %>%
                     tidyr::pivot_wider(names_from = stateID,
                                        values_from = length_counts)

write.table(x = length_counts_df,                                                 # write output
            file = output_path, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = F)

