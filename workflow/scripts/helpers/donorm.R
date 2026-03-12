donorm <- function (dataset, method = "TMM", returnCPM = FALSE)
{
    method <- match.arg(method, c("TMM", "RLE", "upperquartile", "linear", "none"))
    if (method == "none")
        return(dataset)
    if (method == "linear") {
        nf <- edgeR::getLinearNormalizers(dataset)
        return(t(t(dataset)/nf))
    }
 
    d <- edgeR::DGEList(counts = dataset, group = colnames(dataset))
    d <- edgeR::calcNormFactors(d, method = method)
    e <- edgeR::cpm(d, normalized.lib.sizes = T)
    if (!returnCPM)
        e <- e * mean(d$samples$lib.size)/1e+06
    return(as.data.frame(e))
}