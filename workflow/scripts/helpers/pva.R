pva <- function(pv) cbind("none" = pv, 
    "FDR" = p.adjust(pv, "fdr"), "BY" = p.adjust(pv, "BY"),
    "Bonferroni" = p.adjust(pv, "bonferroni"))