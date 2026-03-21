options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

`%>%` <- magrittr::`%>%`

input_path_metadata                                <-  args$input_path_metadata
input_path_library_sizes                           <-  args$input_path_library_sizes
output_path_plot_samples_per_study                 <-  args$output_path_plot_samples_per_study
output_path_plot_samples_per_condition             <-  args$output_path_plot_samples_per_condition
output_path_plot_subtype_distribution              <-  args$output_path_plot_subtype_distribution
output_path_plot_grade_distribution                <-  args$output_path_plot_grade_distribution
output_path_plot_samples_per_study_by_condition    <-  args$output_path_plot_samples_per_study_by_condition
output_path_plot_age_distribution                  <-  args$output_path_plot_age_distribution
output_path_plot_library_size_by_study             <-  args$output_path_plot_library_size_by_study

metadata <- read.table(input_path_metadata, header = T, sep = "\t")
library_sizes <- read.table(input_path_library_sizes, header = T, sep = "\t")

df <- merge(metadata, library_sizes, by = "SRA_id")

samples_per_study <- df %>%
                        dplyr::distinct(study_id, donor_id) %>%
                        dplyr::count(study_id) %>%
                        ggplot2::ggplot(ggplot2::aes(x = study_id, y = n, fill = study_id)) +
                        ggplot2::geom_bar(stat = "identity") +
                        ggplot2::labs(
                            title = "Number of samples per study",
                            x = "Study ID",
                            y = "Number of samples"
                        ) +
                        ggplot2::theme_minimal()

samples_per_condition <- df %>%
                            dplyr::count(condition) %>%
                            ggplot2::ggplot(ggplot2::aes(x = condition, y = n, fill = condition)) +
                            ggplot2::geom_bar(stat = "identity") +
                            ggplot2::labs(
                                title = "Samples by condition",
                                x = "Condition",
                                y = "Number of samples"
                            ) +
                            ggplot2::theme_minimal()

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

ggplot2::ggsave(output_path_plot_samples_per_study, plot = samples_per_study, width = 6, height = 5)
ggplot2::ggsave(output_path_plot_samples_per_condition, plot = samples_per_condition, width = 6, height = 5)
ggplot2::ggsave(output_path_plot_subtype_distribution, plot = subtype_distribution, width = 6, height = 5)
ggplot2::ggsave(output_path_plot_grade_distribution, plot = grade_distribution, width = 6, height = 5)
ggplot2::ggsave(output_path_plot_samples_per_study_by_condition, plot = samples_per_study_by_condition, width = 6, height = 5)
ggplot2::ggsave(output_path_plot_age_distribution, plot = age_distribution, width = 6, height = 5)
ggplot2::ggsave(output_path_plot_library_size_by_study, plot = library_size_by_study, width = 7, height = 5)