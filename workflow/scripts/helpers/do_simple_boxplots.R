do_simple_boxplots <- function(list_of_vecs, x_axis = "", title = ""){

    df <- data.frame(
        value = unlist(list_of_vecs),
        column = rep(names(list_of_vecs), lengths(list_of_vecs))
    )

    #medians <- aggregate(value ~ column, df, median)
    #medians$group <- ifelse(medians$value < median_var, "Below", "Above")
    #df <- merge(df, medians[, c("column", "group")], by = "column")

    p <- ggplot2::ggplot(df, ggplot2::aes(x = column, y = value)) +
        ggplot2::geom_boxplot(color = "black", outlier.alpha = 0.3) +
        #ggplot2::geom_hline(yintercept = median_var, linetype = "dashed", color = "red", linewidth = 1) +
        ggplot2::labs(title = title, x = x_axis, y = "Zero stretch length") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    return(p)
}