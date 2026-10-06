do_get_counts_analysis <- function(counts, scaling_factors, correct_for_CNAs){

if (correct_for_CNAs){
    stopifnot(                                                                   # fix: was a bare `all(...)` expression with the result never checked -
        "counts and scaling_factors are not aligned by binID" =                  # merge()/dplyr::summarise() sort by binID by default, so this likely
        all(counts$binID == scaling_factors$binID)                               # already held in practice, but nothing was actually enforcing it -
    )                                                                            # a defensive, zero-cost safety net in case that default ever changes.
    all(colnames(counts[,colnames(scaling_factors)]) == colnames(scaling_factors))
    counts[,colnames(scaling_factors[,-1])] <- round(counts[,colnames(scaling_factors[,-1])] / scaling_factors[,-1])

}

library_sizes <- colSums(counts[,-1]) # in this case it makes sense to compute lib size like this because it's at the bin level. For peaks we need to use the .BAM files
library_sizes_df <- data.frame(SRA_id = names(library_sizes), library_size = library_sizes, row.names = NULL)

group <- metadata$condition
names(group) <- metadata$SRA_id
group <- group[names(library_sizes)]

metadata <- merge(metadata, library_sizes_df, by="SRA_id")

sample_order <- colnames(counts)[-1]
metadata <- metadata[match(sample_order, metadata$SRA_id), ]
all(metadata$SRA_id == colnames(counts)[-1])

count.f <- normATAC(counts = counts,
		    minSamples_acc = minSamples_atacseq,
		    minCount_acc = minCounts_atacseq,
		    minCpm_acc = minCpm_atacseq,
		    group = group,
		    libSize = library_sizes,
		    normalize = F)

normATAC.out <- normATAC(counts = counts,
                         minSamples_acc = minSamples_atacseq,
                         minCount_acc = minCounts_atacseq,
                         minCpm_acc = minCpm_atacseq,
                         group = group,
                         libSize = library_sizes,                       
                         normalize = T)

count.f.norm <- normATAC.out$counts
binIDs <- normATAC.out$binIDs
count.f.norm$binID <- binIDs
count.f.norm <- count.f.norm[, c("binID", setdiff(names(count.f.norm), "binID"))]

########## PRE SVA

pca <- prcomp(t(count.f.norm[,-1]), scale. = TRUE)
pca_variance_vec <- summary(pca)$importance[2,]

pca_df <- data.frame(SRA_id = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2])
pca_df <- merge(pca_df, metadata, by = "SRA_id")

############ SVA

study <- metadata$study_id
condition <- metadata$condition

mm1 <- model.matrix(~condition + study)
mm0 <- model.matrix(~study)

# NOTE (fix): removed a dead `n.sv <- num.sv(count.f.norm[,-1], mm1, method="leek")` call -
# its result was never referenced anywhere below (sva() is always called with the explicit,
# user-provided `SVs` count, not this estimate). Confirmed unused by inspection. Also purely
# incidentally, this specific call is numerically fragile with few features/samples (fails
# with "infinite or missing values" on small toy data) - removing the dead computation removes
# that fragility too, with no change in behavior since its result was never used.

sv <- sva(as.matrix(count.f.norm[,-1]),mm1,mm0, SVs)$sv

mm <- model.matrix(~sv+condition)
rownames(mm) <- metadata$SRA_id

count.f.norm.sva <- as.data.frame(svaBatchCor(count.f.norm[,-1], mm1, mm0, SVs)$corrected)
count.f.norm.sva$binID <- count.f.norm$binID
count.f.norm.sva <- count.f.norm.sva[, c("binID", setdiff(names(count.f.norm.sva), "binID"))]

pca_sva <- prcomp(t(count.f.norm.sva[,-1]), scale. = TRUE)
pca_sva_variance_vec <- summary(pca_sva)$importance[2,]

pca_sva_df <- data.frame(SRA_id = rownames(pca_sva$x), PC1 = pca_sva$x[,1], PC2 = pca_sva$x[,2])
pca_sva_df <- merge(pca_sva_df, metadata, by = "SRA_id")

############ PLOTS

scree_plot                <- pca_scree_plot(pca_variance_vec)
pca_plot_condition        <- pca_plot(pca, pca_df, hue = "condition")
pca_plot_library_size     <- pca_plot(pca, pca_df, hue = "library_size")
pca_plot_sex              <- pca_plot(pca, pca_df, hue = "sex")
pca_plot_subtype          <- pca_plot(pca, pca_df, hue = "subtype")
pca_plot_donor            <- pca_plot(pca, pca_df, hue = "donor_id", shape = "condition")
pca_plot_study            <- pca_plot(pca, pca_df, hue = "study_id", shape = "condition")

scree_plot_sva            <- pca_scree_plot(pca_sva_variance_vec)
pca_sva_plot_condition    <- pca_plot(pca_sva, pca_sva_df, hue = "condition")
pca_sva_plot_library_size <- pca_plot(pca_sva, pca_sva_df, hue = "library_size")
pca_sva_plot_sex          <- pca_plot(pca_sva, pca_sva_df, hue = "sex")
pca_sva_plot_subtype      <- pca_plot(pca_sva, pca_sva_df, hue = "subtype")
pca_sva_plot_donor        <- pca_plot(pca_sva, pca_sva_df, hue = "donor_id", shape = "condition")
pca_sva_plot_study        <- pca_plot(pca_sva, pca_sva_df, hue = "study_id", shape = "condition")

pca_distance_plot         <- plot_distances_PCA(pca, pca_sva, metadata)
pca_distance_panel        <- plot_PCA_SVA_distance_panel(pca, pca_sva, metadata)

out.list = list(
    count.f=count.f,
    library_sizes_df=library_sizes_df,
    mm=mm,
    scree_plot=scree_plot,
    pca_plot_condition=pca_plot_condition,
    pca_plot_library_size=pca_plot_library_size,
    pca_plot_sex=pca_plot_sex,
    pca_plot_subtype=pca_plot_subtype,
    pca_plot_donor=pca_plot_donor,
    pca_plot_study=pca_plot_study,
    scree_plot_sva=scree_plot_sva,
    pca_sva_plot_condition=pca_sva_plot_condition,
    pca_sva_plot_library_size=pca_sva_plot_library_size,
    pca_sva_plot_sex=pca_sva_plot_sex,
    pca_sva_plot_subtype=pca_sva_plot_subtype,
    pca_sva_plot_donor=pca_sva_plot_donor,
    pca_sva_plot_study=pca_sva_plot_study,
    pca_distance_plot=pca_distance_plot,
    pca_distance_panel=pca_distance_panel
)

return(out.list)
    
}