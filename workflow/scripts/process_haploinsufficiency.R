options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_haplo                     <- args$input_path_haplo
input_path_backbone                  <- args$input_path_backbone
output_path                          <- args$output_path
haplo_thr                            <- as.numeric(args$haplo_thr)
C								     <- as.integer(args$C)


`%>%` <- magrittr::`%>%`

haplo <- read.table(input_path_haplo, header = T, sep = ";")

# questa contiene solo bins con almeno un gene
backbone <- read.table(input_path_backbone, header = T, sep = "\t")
colnames(backbone) <- c("chr", "start", "end", "binID", "annot")

haplo.s <- haplo[, c("transcript","gene","pLI")]
haplo.s.f <- haplo.s[haplo.s$pLI > haplo_thr,]

haplo.s.f.avg <- haplo.s.f %>%
                select(-transcript) %>%
                group_by(gene) %>%
                summarise(gene_avg_pLI = mean(pLI))

backbone_w_haplo <- do.call(rbind,
                            parallel::mclapply(mc.cores = C, 
                                               X = seq(nrow(backbone)),
                                               FUN = function(idx){
                                                row_backbone <- backbone[idx,]
                                                genes <- strsplit(x = row_backbone$annot, split = ",")[[1]]
                                                avg.avg.pLI <- mean(haplo.s.f.avg[haplo.s.f.avg$gene %in% genes, "gene_avg_pLI"]$gene_avg_pLI) # NaN solo se il gene in backbone non è contenuto in haplo
                                                to.return <- row_backbone
                                                to.return$haploinsufficiency_score <- avg.avg.pLI
                                                return(to.return)
                                               }              
                                            )
                                        )

write.table(
    file = output_path,
    x = backbone_w_haplo[,c("binID", "haploinsufficiency_score")],
    sep = "\t",
    col.names = T,
    row.names = F,
    quote = F
)