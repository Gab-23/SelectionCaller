do_simple_histogram <- function(vec, x_axis, title){
    
    df <- data.frame(x = vec)
    qs <- quantile(df$x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = x)) +
          ggplot2::geom_histogram(binwidth = 0.5, fill = "blue", color = "black") +
          ggplot2::geom_vline(xintercept = qs,
                              linetype = "dashed",
                              color = "red",
                              linewidth = 1) +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = title, x = x_axis, y = "Count")

    return(p)
}
