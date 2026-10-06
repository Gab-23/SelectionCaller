options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

`%>%` <- magrittr::`%>%`

input_path_models_DepMap                           <- args$input_path_models_DepMap
input_path_manifest_GDC                            <- args$input_path_manifest_GDC
input_path_samples_COSMIC                          <- args$input_path_samples_COSMIC
input_path_counts_atacseq                          <- args$input_path_counts_atacseq
input_path_metadata_atacseq                        <- args$input_path_metadata_atacseq
input_path_library_sizes                           <- args$input_path_library_sizes
output_path_samples_barplot                        <- args$output_path_samples_barplot
output_path_plot_samples_per_study_by_condition    <- args$output_path_plot_samples_per_study_by_condition
output_path_plot_subtype_distribution              <- args$output_path_plot_subtype_distribution
output_path_plot_grade_distribution                <- args$output_path_plot_grade_distribution
output_path_plot_age_distribution                  <- args$output_path_plot_age_distribution
output_path_plot_library_size_by_study             <- args$output_path_plot_library_size_by_study

# ---- inputs ----
DepMap <- readLines(input_path_models_DepMap)
COSMIC <- readLines(input_path_samples_COSMIC)
GDC <- read.table(input_path_manifest_GDC, header = T, sep = "\t")

# ATAC-seq samples = columns of the final counts matrix (i.e. after samples_to_exclude)
samples_ATAC  <- setdiff(colnames(read.table(input_path_counts_atacseq, header = T, sep = "\t", nrows = 1, check.names = F)), "binID")
metadata_ATAC <- read.table(input_path_metadata_atacseq, header = T, sep = "\t")
library_sizes <- read.table(input_path_library_sizes, header = T, sep = "\t")

missing_cols <- setdiff(c("SRA_id", "donor_id", "condition", "study_id"), colnames(metadata_ATAC))
if (length(missing_cols) > 0) stop("ATAC-seq metadata is missing column(s): ", paste(missing_cols, collapse = ", "))

missing_ATAC <- setdiff(samples_ATAC, metadata_ATAC$SRA_id)
if (length(missing_ATAC) > 0) stop("ATAC-seq samples not in metadata: ", paste(missing_ATAC, collapse = ", "))

# one table, one sample set (the counts-matrix samples) for every ATAC-seq plot below
df <- merge(metadata_ATAC, library_sizes, by = "SRA_id")
df <- df[df$SRA_id %in% samples_ATAC, ]
if (nrow(df) != length(samples_ATAC))
  warning(length(samples_ATAC) - nrow(df), " ATAC-seq samples have no library size and are missing from the descriptive plots")

studies_ATAC    <- as.character(metadata_ATAC$study_id[match(samples_ATAC, metadata_ATAC$SRA_id)])
donors_ATAC     <- metadata_ATAC$donor_id[match(samples_ATAC, metadata_ATAC$SRA_id)]
conditions_ATAC <- metadata_ATAC$condition[match(samples_ATAC, metadata_ATAC$SRA_id)]

if (anyNA(studies_ATAC)) stop(sum(is.na(studies_ATAC)), " ATAC-seq samples have no study_id")
if (anyNA(donors_ATAC)) warning(sum(is.na(donors_ATAC)), " ATAC-seq samples have no donor_id and are not counted in 'Subjects'")
if (!all(conditions_ATAC %in% c("tumor", "NAT"))) warning("ATAC-seq samples with condition other than tumor/NAT: Tumour + NAT will not add up to All")

# unique subjects = unique (study_id, donor_id): robust to the same donor_id being reused across studies
n_subjects_ATAC <- nrow(unique(data.frame(study = studies_ATAC, donor = donors_ATAC)[!is.na(donors_ATAC), ]))

# ---- samples barplot: cohorts panel + ATAC-seq panel ----
cohort_df <- data.frame(Dataset = c("DepMap", "COSMIC", "GDC"),
                        Samples = c(length(DepMap), length(COSMIC), nrow(GDC)), panel = "cohorts")

# ATAC-seq panel, fixed order: All | Subjects | study 1 | study 2 | ... | Tumour | NAT
atac_bar <- function(keep, label){
  data.frame(Dataset = label, Samples = sum(keep), panel = "ATAC")
}
study_levels <- names(sort(table(studies_ATAC), decreasing = TRUE))                # largest study first; set manually for a custom order
atac_df <- rbind(
  atac_bar(rep(TRUE, length(samples_ATAC)),  "All"),
  data.frame(Dataset = "Subjects", Samples = n_subjects_ATAC, panel = "ATAC"),
  do.call(rbind, lapply(study_levels, function(s) atac_bar(studies_ATAC == s, s))),
  atac_bar(conditions_ATAC == "tumor", "Tumour"),
  atac_bar(conditions_ATAC == "NAT",   "NAT")
)
samples_df <- rbind(cohort_df, atac_df)

if (anyDuplicated(samples_df$Dataset) > 0)
  stop("Duplicated bar labels (a study_id clashes with another bar): ", paste(samples_df$Dataset[duplicated(samples_df$Dataset)], collapse = ", "))

cohort_levels <- cohort_df$Dataset[order(-cohort_df$Samples)]                       # cohorts sorted by size, ATAC-seq bars in the fixed order above
samples_df$Dataset <- factor(samples_df$Dataset, levels = c(cohort_levels, "All", "Subjects", study_levels, "Tumour", "NAT"))
samples_df$panel   <- factor(samples_df$panel, levels = c("cohorts", "ATAC"))

study_palette <- c("#0072B2", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7", "#999999")
fill_colors <- c(setNames(c("#F8766D", "#00BA38", "#619CFF")[seq_along(cohort_df$Dataset)], cohort_df$Dataset),
                 "All" = "#E69F00", "Subjects" = "#4D4D4D",
                 setNames(rep_len(study_palette, length(study_levels)), study_levels),
                 "Tumour" = "#b8077d", "NAT" = "#9623b3")

p.barplot_samples <- ggplot2::ggplot(samples_df, ggplot2::aes(x = Dataset, y = Samples, fill = Dataset)) +
  ggplot2::geom_col(color = "black") +
  ggplot2::geom_text(ggplot2::aes(label = Samples), vjust = -0.4, size = 3.5) +
  ggplot2::scale_fill_manual(values = fill_colors) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
  ggplot2::facet_wrap(~panel, scales = "free", strip.position = "bottom",     # separate y axes: cohorts vs ATAC-seq
                      labeller = ggplot2::as_labeller(c(cohorts = "", ATAC = "ATAC-seq"))) +
  ggplot2::labs(title = "", x = "", y = "# Samples / subjects") +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none", strip.placement = "outside", strip.text = ggplot2::element_text(size = 11),
                 panel.grid.minor = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

message("Samples barplot data:\n", paste(utils::capture.output(print(samples_df)), collapse = "\n"))
message("ATAC-seq: ", length(samples_ATAC), " samples from ", n_subjects_ATAC, " subjects, ", length(study_levels), " studies")

# ---- ATAC-seq cohort descriptive plots ----
subtype_distribution <- df %>%
                            dplyr::filter(condition == "tumor") %>%
                            dplyr::count(subtype) %>%
                            ggplot2::ggplot(ggplot2::aes(x = subtype, y = n, fill = subtype)) +
                            ggplot2::geom_bar(stat = "identity") +
                            ggplot2::labs(
                                title = "Tumor subtype distribution",
                                x = "Subtype",
                                y = "Number of samples"
                            ) +
                            ggplot2::theme_minimal() +
                            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

grade_distribution <- df %>%
                        dplyr::filter(condition == "tumor") %>%
                        dplyr::count(grade) %>%
                        ggplot2::ggplot(ggplot2::aes(x = grade, y = n, fill = grade)) +
                        ggplot2::geom_bar(stat = "identity") +
                        ggplot2::labs(
                            title = "Tumor grade distribution",
                            x = "Grade",
                            y = "Number of samples"
                        ) +
                        ggplot2::theme_minimal()

samples_per_study_by_condition <- df %>%
                                    dplyr::distinct(study_id, donor_id, condition) %>%
                                    dplyr::count(study_id, condition) %>%
                                    ggplot2::ggplot(ggplot2::aes(x = study_id, y = n, fill = condition)) +
                                    ggplot2::geom_bar(stat = "identity", position = "dodge") +
                                    ggplot2::labs(
                                        title = "Samples per study by condition",
                                        x = "Study ID",
                                        y = "Number of samples"
                                    ) +
                                    ggplot2::theme_minimal()

age_distribution <- df %>%
                        dplyr::mutate(age_bin = cut(age, breaks = c(40, 50, 60, 70, 80))) %>%
                        dplyr::count(age_bin) %>%
                        ggplot2::ggplot(ggplot2::aes(x = age_bin, y = n, fill = age_bin)) +
                        ggplot2::geom_bar(stat = "identity") +
                        ggplot2::labs(
                            title = "Age distribution",
                            x = "Age range",
                            y = "Number of patients"
                        ) +
                        ggplot2::theme_minimal()

library_size_by_study <- df %>%
    ggplot2::ggplot(ggplot2::aes(x = study_id, y = library_size, fill = study_id)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
    ggplot2::labs(
        title = "Library size by study",
        x = "Study ID",
        y = "Library size"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )

ggplot2::ggsave(output_path_samples_barplot, plot = p.barplot_samples, width = 8, height = 4)
ggplot2::ggsave(output_path_plot_samples_per_study_by_condition, plot = samples_per_study_by_condition, width = 6, height = 5)
ggplot2::ggsave(output_path_plot_subtype_distribution, plot = subtype_distribution, width = 6, height = 5)
ggplot2::ggsave(output_path_plot_grade_distribution, plot = grade_distribution, width = 6, height = 5)
ggplot2::ggsave(output_path_plot_age_distribution, plot = age_distribution, width = 6, height = 5)
ggplot2::ggsave(output_path_plot_library_size_by_study, plot = library_size_by_study, width = 7, height = 5)