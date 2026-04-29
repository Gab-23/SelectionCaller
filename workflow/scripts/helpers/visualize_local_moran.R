visualize_local_moran <- function(df, feature){

df <- df %>%
    tidyr::separate(binID, into = c("chr", "pos"), sep = "_") %>%
    dplyr::mutate(chr = as.numeric(chr), pos = as.numeric(pos)) %>%
    dplyr::arrange(chr, pos)

p <- ggplot2::ggplot(df, ggplot2::aes(x = pos, y = I, color = hotspot)) +
        ggplot2::geom_point(size = 1.5) +
        ggplot2::facet_wrap(~ chr, scales = "free_x", nrow = 1) +
        ggplot2::theme_bw() +
        ggplot2::labs(
            title = paste("Local Moran's I (Manhattan-style):", feature),
            x = "Position",
            y = "Local Moran's I"
        )

return(p)

}