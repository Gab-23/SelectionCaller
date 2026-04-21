fraction_plot <- function(df, y_axis){
    p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = stats::reorder(variable, -frac), y = frac)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Variable", y = paste0("Fraction of ", y_axis)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

    return(p)
}