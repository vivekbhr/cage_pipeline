#!/usr/bin/Rscript
## demult input data (R script)
# Arg 1 : previously saved CS object
# Arg 2 : path to mapped BAMs
# Arg 3 : output dir
# Args 4 : Number of threads
setwd("/data/processing4/jana/cage_pipeline")



library(icetea)
Args = commandArgs(TRUE)
csobject <- Args[1]
mappeddir <- Args[2]
outdir <- Args[3]
cores <- Args[4]

load(csobject) # CSobject -> cs
si <- sampleInfo(cs)
si$mapped_file <- file.path(mappeddir, paste0(si$samples, ".bam"))

message(paste0("Files used: ", si$mapped_file))

# calc num mapped reads
message("Counting mapped reads")
si$num_mapped <- sapply(si$mapped_file, function(x) {
    print(x)
    Rsamtools::countBam(x, param = Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE,
                                  isFirstMateRead = TRUE,
                                  isNotPrimaryRead = FALSE
                                  )) )[,6]
    })

# assign sampleinfo
si -> sampleInfo(cs)
message("removing duplicates")
system.time({
    cs <- filterDuplicates(cs, outdir,
                           ncores = cores,
                           keepPairs = FALSE)
})

message("writing CSobject")
save(cs, file = file.path(outdir, "CSobject.Rdata"))
