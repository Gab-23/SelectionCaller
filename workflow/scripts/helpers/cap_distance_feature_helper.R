cap_distance_feature <- function(df, dist_cap){                                    # df: binID + one distance column
    dist_col <- setdiff(colnames(df), "binID")
    if (length(dist_col) != 1) stop("expected binID + exactly one distance column, got: ", paste(colnames(df), collapse = ", "))
    df[[dist_col]] <- pmin(df[[dist_col]], dist_cap)                                # NA stays NA (na.rm = FALSE)
    return(df)
}

split_binID <- function(binID){                                                    # "17_123" -> chr "17", bin 123
    parts <- do.call(rbind, strsplit(binID, "_"))
    chr_order <- unique(parts[order(suppressWarnings(as.numeric(parts[,1])), parts[,1]), 1])   # 1..22, then X, Y
    return(data.frame(chr = factor(parts[,1], levels = chr_order), bin = as.numeric(parts[,2])))
}

contiguous_run <- function(chr, bin){                                              # run ID per stretch of consecutive bins
    o <- order(chr, bin)                                                           # (bins filtered out of the backbone leave gaps in the index)
    new_run <- c(TRUE, chr[o][-1] != chr[o][-length(o)] | diff(bin[o]) != 1)
    run <- integer(length(bin)); run[o] <- cumsum(new_run)
    return(run)
}

plot_capping_diagnostics <- function(before, after, dist_cap){
    dist_col <- setdiff(colnames(before), "binID")
    long <- rbind(data.frame(split_binID(before$binID), dist = before[[dist_col]], version = "before"),
                  data.frame(split_binID(after$binID),  dist = after[[dist_col]],  version = "after"))

    plot_profile <- function(v){                                                    # one page per version: "after" is unreadable on the "before" scale
        d <- long[long$version == v, ]
        d$run <- contiguous_run(d$chr, d$bin)                                       # no line drawn across a gap in the backbone
        ggplot2::ggplot(d, ggplot2::aes(x = bin, y = dist / 1e6, group = run)) +
            ggplot2::geom_line(linewidth = 0.3, na.rm = TRUE) +
            ggplot2::geom_hline(yintercept = dist_cap / 1e6, linetype = "dashed", colour = "firebrick", linewidth = 0.3) +
            ggplot2::facet_wrap(~chr, scales = "free_x", drop = FALSE) +
            ggplot2::labs(title = paste0(dist_col, ": ", v, " capping (cap = ", dist_cap / 1e6, " Mbp, dashed)"),
                          x = "bin index", y = "distance (Mbp)") +
            ggplot2::theme_bw()
    }

    a <- long[long$version == "after", ]
    per_chr <- do.call(rbind, lapply(split(a, a$chr), function(x){
        data.frame(chr = x$chr[1],
                   metric = c("fraction of bins at cap", "fraction of bins NA (no locus on chr)"),
                   value  = c(mean(x$dist >= dist_cap, na.rm = TRUE), mean(is.na(x$dist))))
    }))
    per_chr$value[is.nan(per_chr$value)] <- NA                                      # all-NA chromosome: no bin can be at cap

    summary <- ggplot2::ggplot(per_chr, ggplot2::aes(x = chr, y = value)) +
                    ggplot2::geom_col(na.rm = TRUE) +
                    ggplot2::facet_wrap(~metric, ncol = 1) +
                    ggplot2::coord_cartesian(ylim = c(0, 1)) +
                    ggplot2::labs(title = paste0(dist_col, ": per-chromosome summary after capping"), x = "chromosome", y = NULL) +
                    ggplot2::theme_bw()

    return(list(before = plot_profile("before"), after = plot_profile("after"), summary = summary))
}