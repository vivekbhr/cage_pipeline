
## plot read statistics from the CapSet object
library(icetea)
## INPUT OUTPUT
Args = commandArgs(TRUE)
# Input/Output filename
capset_obeject <- Args[1]
plotFile_prefix <- Args[2]

load(capset_obeject)
# plot read Numbers
plotReadStats(cs, plotType = "dodge", plotValue = "numbers", outFile = paste0(plotFile_prefix, "_numbers.pdf"))
# plot proportions
plotReadStats(cs, plotType = "stack", plotValue = "proportions", outFile = paste0(plotFile_prefix, "_proportions.pdf"))

