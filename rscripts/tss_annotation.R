

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
# get GTF file
GTFfile <- Args[4]

# Get User supplied annotations
enhancer_file <- Args[5]
repeat_file <- Args[6]
dhs_file <- NA


## PREPARE ANNOTATIONS
# create Txdb of flybase GTF
gtf <- makeTxDbFromGFF(GTFfile)

# Get Annotations from TxDb
columns <- c("tx_name", "gene_id")
known_tss <- unique(promoters(gtf, upstream = 1, downstream = 1, columns = columns))
known_promoters <- unique(promoters(gtf, upstream = 500, downstream = 0, columns = columns))
txbygene <- unique(unlist(transcriptsBy(gtf, "gene")))

# Import other annotations
if (!is.na(enhancer_file)) enhancers <- import.bed(enhancer_file)
if (!is.na(repeat_file)) repeats <- import.bed(repeat_file)
if (!is.na(dhs_file)) dhs <- import.bed(dhs_file)


## DEFINE FUNCTIONS
annotateTSS <- function(tssBED,
                        txdb,
                        ignoreStrand = FALSE,
                        featureRank = c("fiveUTR",
                                        "promoter",
                                        "coding",
                                        "intron",
                                        "spliceSite",
                                        "threeUTR",
                                        "intergenic")
                        ){
    ## resolve 1:many mapping issue by prioritising some features over others
    stopifnot(length(featureRank) == 7)
    rankvec <- seq_len(7)
    names(rankvec) <- featureRank
    tss <- tssBED
    # Annotate
    suppressWarnings({
        suppressMessages({
        db <- locateVariants(
        query = tss,
        subject = txdb,
        ignore.strand = ignoreStrand,
        AllVariants(promoter = PromoterVariants(
            upstream = 500, downstream = 0
            ))
        )
        })
    })

    ## resolve 1:many mapping isues using ranks from rankdf
    db$rank <- vapply(as.character(db$LOCATION),
                      getranks,
                      rank_vec = rankvec,
                      FUN.VALUE = numeric(length = 1))
    return(db)
}

getranks <- function(x, rank_vec) {
                    return(rank_vec[names(rank_vec) == x])
}

# split the entries by rank and only keep unique highest rank per entry
splitranks <- function(x) {
    x <- as.data.frame(x)
    x <- x[,c(1:6,9:10,12:15)]
    l <- lapply(split(x, x$QUERYID), unique)
    l2 <- lapply(l, function(y) {
        return(y[which(y$rank == min(y$rank)), ])
    })
    l3 <- do.call(rbind, l2)
    return(l3)
}


###########################  ANNOTATE  ########################################
### First Annotate TSS

# Import bed and annotate TSS first
bed <- import.bed(input_bed)
bed$name <- paste0("region", 1:length(bed))
ol <- findOverlaps(bed, known_tss)
leftover_bed <- bed[-queryHits(ol)]

## STEP 1 : anotate TSS entries
message("Annotating TSS")
tss_bed <- as.data.frame(bed[queryHits(ol)])
tss_bed$LOCATION <- "TSS"
tss_bed$TXNAME <- as.character(known_tss[subjectHits(ol)]$tx_name)
tss_bed$GENEID <- as.character(known_tss[subjectHits(ol)]$gene_id)
tss_bed$score <- NULL
tss_bed$PRECEDEID <- "None"
tss_bed$FOLLOWID <- "None"

## unique entries left to annotate
left1 <- length(unique(bed$name)) - length(unique(tss_bed$name))
print("entries left after step1")
print(left1)
print("confirm!")
print(length(unique(leftover_bed$name)))


### Then annotate other GTF features
message("Annotating other features from the GTF")
## STEP 2: annotate sense features in genes
bed_annot <- annotateTSS(leftover_bed, txdb = gtf)
bed_annot2 <- splitranks(bed_annot)
# convert TXID to name
bed_annot2$TXID <- select(gtf, keys = as.character(bed_annot2$TXID), columns = "TXNAME", keytype = "TXID")$TXNAME

print("entries left after step2")
done <- as.data.frame(bed_annot2)[1:5] %>% unique() %>% nrow()
left2 <- left1 - done
print(left2)

## unique entries left to annotate
unique(bed_annot$QUERYID) -> qid
x <- 1:length(leftover_bed)
leftover_ids <- x[!(x %in% qid)]
print("confirm!")
print(length(unique(leftover_bed$name)) - length(qid))


## STEP 3: annotate antisense features in genes
message("annotate antisense features in genes")
leftover_bed2 <- leftover_bed[leftover_ids]
# check if antisense RNA is there (if leftover entried overlap on the other strand)
leftover_annot <- annotateTSS(leftover_bed2, txdb = gtf, ignoreStrand = TRUE)
bed_annot3 <- splitranks(leftover_annot)
bed_annot3$TXID <- select(gtf, keys = as.character(bed_annot3$TXID), columns = "TXNAME", keytype = "TXID")$TXNAME

print("entries left after step 3")
done <- as.data.frame(bed_annot3)[1:5] %>% unique() %>% nrow()
left3 <- left2 - done
print(left3)

# unique entries left to annotate?
unique(leftover_annot$QUERYID) -> qid
x <- 1:length(leftover_bed2)
leftover_ids <- x[!(x %in% qid)]
print("confirm!")
print(length(unique(leftover_bed2$name)) - length(qid))

## STEP 4: annotate the leftover by manual intersect with transcript
leftover_bed3 <- leftover_bed2[leftover_ids]
hits <- findOverlaps(leftover_bed3, txbygene )
bed_annot4 <- leftover_bed3[queryHits(hits)]
bed_annot4$LOCATION <- "sense-transcript"


# try to add Tx and gene ID as well
hits2 <- findOverlaps(bed_annot4, txbygene)
hits2 %<>% as.data.frame()
bed_annot4 %<>% as.data.frame()

bed_annot4$name <- NULL
bed_annot4$score <- NULL

split(hits2, list(hits2$queryHits)) %>%
    lapply(function(x) paste(unique(x$TXNAME), collapse = ";")) %>%
    unlist() -> bed_annot4$TXNAME
split(hits2, list(hits2$queryHits)) %>%
    lapply(function(x) paste(unique(x$GENEID), collapse = ";")) %>%
    unlist() -> bed_annot4$GENEID
bed_annot4$PRECEDEID <- "None"
bed_annot4$FOLLOWID <- "None"

# unique entries left to annotate?
print("entries left after step 4")
done <- as.data.frame(bed_annot4)[1:5] %>% unique() %>% nrow()
left4 <- left3 - done
print(left4)

# annotate the rest as intergenic
bed_annot5 <- as.data.frame(leftover_bed3[-queryHits(hits)])
if (nrow(bed_annot5) > 0) {
    print("Annotating the rest as `intergenic`")
    bed_annot5 <- as.data.frame(bed_annot5)[1:5]
    bed_annot5$LOCATION <- "intergenic"
    use_leftover <- TRUE
} else {
    use_leftover <- FALSE
}


### Finally, annotate Intergenics/Introns
message("Further re-annotating intergenic and intronic regions!")

## STEP 5-8 : Annotate intergenic and intronic regions into eRNA/repeat/DHS-rna/novel
bed_intron_intergenic <- bed_annot2[bed_annot2$LOCATION %in% c("intergenic", "intron"),]
# also add antisense
bed_intron_intergenic <- rbind(bed_intron_intergenic,
                               bed_annot3[bed_annot3$LOCATION %in% c("intergenic", "intron"),])

# remove entries from corresponding beds
bed_annot2 <- bed_annot2[!(bed_annot2$LOCATION %in% c("intergenic", "intron")),]
bed_annot3 <- bed_annot3[!(bed_annot3$LOCATION %in% c("intergenic", "intron")),]

## Finally, add sense/antisense definition to them
if (nrow(bed_annot2) > 0) bed_annot2$LOCATION <- paste0("sense-",bed_annot2$LOCATION)
if (nrow(bed_annot3) > 0) bed_annot3$LOCATION <- paste0("antisense-",bed_annot3$LOCATION)


## add bed_annot5 if it's non-empty
if (isTRUE(use_leftover)) {
    bed_annot5$QUERYID <- 0
    bed_annot5$GENEID <- "None"
    bed_annot5$PRECEDEID <- "None"
    bed_annot5$FOLLOWID <- "None"
    bed_annot5$rank <- 0
    bed_intron_intergenic <- rbind(bed_intron_intergenic, bed_annot5)
    bed_annot5 <- NULL
} else {
    bed_annot5 <- NULL
}

print(paste0("Entries used for re-annotation: ", nrow(bed_intron_intergenic) ))

## STEP 5 : eRNA
bed_intron_intergenic <- GRanges(seqnames = bed_intron_intergenic$seqnames,
                                 ranges = IRanges(start = bed_intron_intergenic$start,
                                                  end = bed_intron_intergenic$end),
                                 strand = bed_intron_intergenic$strand,
                                 LOCATION = bed_intron_intergenic$LOCATION,
                                 TXNAME = bed_intron_intergenic$QUERYID,
                                 GENEID = bed_intron_intergenic$GENEID,
                                 PRECEDEID = bed_intron_intergenic$PRECEDEID,
                                 FOLLOWID = bed_intron_intergenic$FOLLOWID,
                                 rank = bed_intron_intergenic$rank)

suppressWarnings(hits <- findOverlaps(bed_intron_intergenic, enhancers))
bed_annot6 <- unique(bed_intron_intergenic[queryHits(hits)])
if(length(bed_annot6) > 0) {
   bed_annot6$LOCATION <- paste0(bed_annot6$LOCATION, "-enhancer")
}

print("eRNA annotated : ")
print(length(bed_annot6))

# get leftOvers
leftover_bed4 <- bed_intron_intergenic[!(bed_intron_intergenic$QUERYID %in% bed_annot6$QUERYID)]


## STEP 6 : annotate repeats
message("Annotating Repeats (if provided)")
suppressWarnings(hits <- findOverlaps(leftover_bed4, repeats))
bed_annot7 <- unique(leftover_bed4[queryHits(hits)])
if (length(bed_annot7) > 0) {
    bed_annot7$LOCATION <- paste0(bed_annot7$LOCATION, "-repeat")
    bed_annot7$TXNAME <- "None"
    bed_annot7$GENEID <- "None"
    bed_annot7$PRECEDEID <- "None"
    bed_annot7$FOLLOWID <- "None"
}

print("repeats annotated : ")
print(length(bed_annot7))

# get leftOvers
leftover_bed5 <- leftover_bed4[!(leftover_bed4$QUERYID %in% bed_annot7$QUERYID)]


## STEP 7 : annotate DHS
message("Annotating DHS sites (if provided)")

if (!is.na(dhs_file) & length(leftover_bed5) > 0) {
    suppressWarnings(hits <- findOverlaps(leftover_bed5, dhs))
    bed_annot8 <- unique(leftover_bed5[queryHits(hits)])
    if(length(bed_annot8) > 0) {
        bed_annot8$LOCATION <- paste0(bed_annot8$LOCATION, "-DHS")
        bed_annot8$TXNAME <- "None"
        bed_annot8$GENEID <- "None"
        bed_annot8$PRECEDEID <- "None"
        bed_annot8$FOLLOWID <- "None"
    }

    print("DHS annotated : ")
    print(length(bed_annot8))
    # get leftOvers (final bed)
    bed_annot9 <- leftover_bed5[!(leftover_bed5$QUERYID %in% bed_annot5$QUERYID)]
} else {
        bed_annot8 <- NULL
        bed_annot9 <- NULL
    }

## FINAL STEP : MERGE ALL ANNOTATIONS
# tss_bed : TSS
# bed_annot2 : sense features
# bed_annot3 : antisense features
# bed_annot4 : unmapped features within transcripts
# bed_annot5 : NULL (leftover marked as intergenics)
# bed_annot6 : eRNAs
# bed_annot7 : repeats
# bed_annot8 : DHS
# bed_annot9 : final leftOver bed (intergenic/intron)
message("Finally, merging the annotated entries!")
all_annot <- list(bed_annot1 = tss_bed,
                  bed_annot2 = bed_annot2,
                  bed_annot3 = bed_annot3,
                  bed_annot4 = bed_annot4,
                  bed_annot5 = bed_annot5,
                  bed_annot6 = bed_annot6,
                  bed_annot7 = bed_annot7,
                  bed_annot8 = bed_annot8,
                  bed_annot9 = bed_annot9)

all_annot %<>% lapply(as.data.frame)
all_annot <- all_annot[sapply(all_annot, function(x) dim(x)[1]) > 0]

## remove queryID column from all entries
all_annot %<>% lapply(function(x){
  x$QUERYID <- NULL
  x$name <- NULL
  x$rank <- NULL
  return(x)
} )

newname <- colnames(all_annot$bed_annot1)
all_annot %<>% lapply(function(x) {
    colnames(x) <- newname
    return(x)
    })

## paste all entries together
final_annot <- Reduce(function(x,y) rbind(x,y), all_annot)
tot <- final_annot[1:5] %>% unique() %>% nrow()

message(paste0("DONE! All unique entries annotated : ", tot))
print("Table of all entries (might be non-unique)")
annot_features <- unique(final_annot[1:6])
print(table(annot_features$LOCATION))
if (!is.na(plotFile)) {
    df <- as.data.frame(table(annot_features$LOCATION))
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

