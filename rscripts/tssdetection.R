#!/usr/bin/Rscript
## detect TSS from filtered files

library(icetea)
Args = commandArgs(TRUE)
csobject <- Args[1] # previously saved CS object
sample_info <- Args[2] # sample info tsv file (with "group" information as 3rd column)
outdir <- Args[3] # output dir
cores <- as.numeric(Args[4]) # number of threads
fold_change <- as.numeric(Args[5]) # fold change threshold for TSS detection
prefix <- Args[6] # 'tss'

df <- read.delim(sample_info, stringsAsFactors = FALSE)
colnames(df) <- c("barcode", "sample", "group")

load(csobject) # CSobject -> cs
# detect TSS
message("Detecting TSS")
cs <- detectTSS(cs, 
                groups = df$group,
                outfile_prefix = file.path(outdir, prefix), 
                ncores = cores, 
                foldChange = fold_change )
# export TSS as bed files
exportTSS(cs, pergroup = TRUE, outfile_prefix =  file.path(outdir, prefix))

message("writing CSobject")
save(cs, file = file.path(outdir, "CSobject.Rdata"))
