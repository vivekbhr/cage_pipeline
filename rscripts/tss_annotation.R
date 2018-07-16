# load the package
library(GenomicFeatures)
library(VariantAnnotation)
library(rtracklayer)
library(magrittr)
library(ggplot2)

#GOAL : annotate the TSS positions by functional categories:

# - known TSS : direct overlap with annotated FlyBase gene TSSs
# - Sense
#    - Intron, CDS, 3UTR, 5UTR, promoter_upstream
# - Antisense
#    - Intron, CDS, 3UTR, 5UTR, promoter_upstream (PROMPT)
# - eRNA : if it overlaps with enhancer annotation (Kvon_2014 data)
# - Intergenic : None of above
#    - repeat, transposon, DHS site, novel

# for promoter, the criteria is <500bp upstream of TSS

## Input files :
#    - Bed file for TSS
#    - Latest flybase GTF
#    - enhancer bed file
#    - repeats bed file
#    - DHS bed file

## Output : .tsv file with bed entries + extra entries for annotation

# From [andersson group paper](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gky244/4962481). :
# PROMPTs : TCs <500 bp upstream of and antisense to annotated FlyBase gene TSSs
# ncRNAs : TCs associated with annotated FlyBase ncRNA TSSs
# eRNAs : TCs associated with gene TSS-distal dCP STARR-seq enhancers


## INPUT OUTPUT
Args = commandArgs(TRUE)
# Input/Output filename
input_bed <- Args[1]
outfile <- Args[2]
plotFile <- Args[3]
# get Genic annotations
annotation_folder <- Args[4]

# Get User supplied external annotations
enhancer_file <- Args[5]
repeat_file <- Args[6]
dhs_file <- NA

# get Genic annotations
flist <- c("TSS", "fiveUTRs", "CDS", "threeUTRs", "introns", 
           "antisense_threeUTRs", "antisense_fiveUTRs", "antisense_transcript", "intergenic")
annoFeatures <- lapply(flist, function(f) {
    print(f)
    path = file.path(annotation_folder, paste0(f, ".bed")) 
    return(import.bed(path))
})
names(annoFeatures) <- flist

## get External annotations
external_annotations <- list(enhancers = import.bed(enhancer_bed),
                             repeats = import.bed(repeats_bed))

if (!(is.na(dhs_bed))) external_annotations$dhs = import.bed(dhs_bed)

## bed file to annotate
bed <- import.bed(input_bed)

## function to annotate features and pass the rest
annotateAndSubset <- function(bed, feature) {
    olap <- findOverlaps(bed, feature)
    olap <- olap[!(duplicated(queryHits(olap)))]
    overlapping <- bed[queryHits(olap)]
    overlapping$name <- feature[subjectHits(olap)]$name
    non_overlapping <- bed[-queryHits(olap)]
    return(list(overlapping = overlapping,
                non_overlapping = non_overlapping))
}

## STEP 1 : GENIC ANNOTATION
# annotate in a loop 
all_annotations <- list()
b <- bed
for (f in names(annoFeatures) ) {
    ann <- annotateAndSubset(bed = b, feature = annoFeatures[[f]])
    ann$overlapping$LOCATION <- f
    all_annotations[[f]] <- ann$overlapping
    b <- ann$non_overlapping
}

# leftOver entries 
b$LOCATION <- "None"
all_annotations$None <- b
## convert to GR
all_annotations <- GRangesList(all_annotations) %>% unlist()

## STEP 2 : EXTERNAL ANNOTATIONS
# annotate in a loop 
all_annotations_step2 <- list()
b <- all_annotations
for (f in names(external_annotations) ) {
    ann <- annotateAndSubset(bed = b, feature = external_annotations[[f]])
    
    if(length(ann$overlapping) != 0) {
        ann$overlapping$LOCATION <- paste(ann$overlapping$LOCATION, f, sep = "-")
        all_annotations_step2[[f]] <- ann$overlapping
        b <- ann$non_overlapping
    } 
}

# leftOver entries 
all_annotations_step2$None <- b
## convert to GR
all_annotations_step2 <- GRangesList(all_annotations_step2) %>% unlist()

### make final
### 
final_annot <- all_annotations_step2
print(table(final_annot$LOCATION))
if (!is.na(plotFile)) {
    df <- as.data.frame(table(final_annot$LOCATION))
    colnames(df) <- c("feature", "number")
    df$percent <- (df$number/sum(df$number))*100
    df$feature <- factor(df$feature, levels = df[order(df$number), "feature"])
    max_number <- round(max(df$number), -1)
    
    ggplot(df, aes(feature, number, fill = "feature")) +
        geom_bar(stat = "identity") +
        scale_fill_brewer(type = "qual") +
        scale_y_continuous(breaks = seq(0, max_number, max_number/10)) +
        coord_flip() +
        theme_bw(base_size = 16) +
        theme(legend.position = "none")
    ggsave(filename = plotFile)
}

write.table(as.matrix(final_annot), file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


