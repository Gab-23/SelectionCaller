# taken from https://www.biostars.org/p/478444/

require(BSgenome.Hsapiens.UCSC.hg19)
require(BSgenome)
library(Repitools)

backbone <- rtracklayer::import.bed(input_path)
backbone.gr <- makeGRangesFromDataFrame(data.frame(backbone))
gc <- gcContentCalc(backbone.gr , organism=Hsapiens, verbose=F)