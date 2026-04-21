feature_boxplots <- function(df){
    p <- ggplot2::ggplot(df, ggplot2::aes(x = variable, y = value)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Variable", y = "Value") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::coord_cartesian(ylim = c(-3, 6))

    return(p)
}