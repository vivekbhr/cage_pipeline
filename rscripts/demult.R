#!/usr/bin/Rscript
## demult input data (R script)
# Arg 1 : experiment name
# Arg 2,3 : R1 and R2 for de-multiplexing
# Arg 4 : tsv file with barcodes (col1) and sample names(col2)
# Arg 5 : output dir
# Arg 6 : No. of cores
library(icetea)
args<-commandArgs(TRUE)
df <- read.delim(args[4], stringsAsFactors = FALSE)

colnames(df) <- c("barcode", "sample", "group")
cs <- newCapSet(expMethod = args[1],
                fastq_R1 = args[2],
                fastq_R2 = args[3],
                idxList = df$barcode,
                sampleNames = df$sample)

cs <- demultiplexFASTQ(cs, outdir = args[5], max_mismatch = 1, ncores = args[6])
save(cs, file = file.path(args[5], "CSobject.Rdata"))
