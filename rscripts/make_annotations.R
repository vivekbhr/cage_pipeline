# load the packages
library(GenomicFeatures)
library(rtracklayer)
library(magrittr)

## load Args
Args <- commandArgs(TRUE)
GTFfile <- Args[1] ## GTF file to extract annotations from
output_dir <- Args[2] ## path to write extracted annotations to

dir.create(output_dir)
## function to remove overlaping region
removeOverlapRegion <- function(query, ref, ignoreStrand = FALSE) {
    hits <- findOverlaps(query, ref, ignore.strand = ignoreStrand)
    grl <- extractList(ref, as(hits, "List"))
    if(isTRUE(ignoreStrand)) strand(grl) <- "*"
    removed <- psetdiff(query, grl) %>% unlist() %>% unique()
    return(removed)
}


## PREPARE ANNOTATIONS
# create Txdb of flybase GTF
gtf <- makeTxDbFromGFF(GTFfile)

## take out genes
dm6genes <- genes(gtf)
export.bed(dm6genes, con = file.path(output_dir, "genes.bed"))

## get intergenic regions
dm6genes.nostrand <- dm6genes
strand(dm6genes) <- "*"
gaps(dm6genes) %>% export.bed(con = file.path(output_dir, "intergenic.bed"))

## get TSS (-10 to +100 bp)
transcriptsBy(gtf, "gene") %>%
    promoters(upstream = 10, downstream = 100) %>%
    unlist() %>% unique() -> prombyTx
export.bed(prombyTx, con = file.path(output_dir, "TSS.bed"))

## get 5 UTR
fiveUTRsByTranscript(gtf) %>% unlist() %>% unique()  -> fiveUTR
# Remove the part of fiveUTRs that overlaps with TSS
removeOverlapRegion(fiveUTR, prombyTx) %>% export.bed(file.path(output_dir, "fiveUTRs.bed"))

### get CDS
cdsBy(gtf, "gene") %>% unlist() -> dm6cds
# Remove the part of CDS that overlaps with 5UTR or TSS
removeOverlapRegion(dm6cds, fiveUTR) %>%
    removeOverlapRegion(prombyTx) %>%
        export.bed(file.path(output_dir, "CDS.bed"))

## get 3 UTR
threeUTRsByTranscript(gtf) %>% unlist() %>% unique() -> threeUTR
# Remove the part of 3UTRs that overlaps with a 5UTR/CDS or TSS
removeOverlapRegion(threeUTR, fiveUTR) %>%
    removeOverlapRegion(dm6cds) %>%
    removeOverlapRegion(prombyTx) %>%
    export.bed(file.path(output_dir, "threeUTRs.bed"))

## get Introns
intronsByTranscript(gtf) %>% unlist() %>% unique() -> dm6introns
# Remove the part of Introns that overlaps with a 5UTR/CDS/3UTR or TSS
removeOverlapRegion(dm6introns, fiveUTR) %>%
    removeOverlapRegion(dm6cds) %>%
    removeOverlapRegion(threeUTR) %>%
    removeOverlapRegion(prombyTx) %>%
    export.bed(file.path(output_dir, "introns.bed"))

################# ANTISENSE FEATURES   ###################
## get antisense 3 UTR
anti_threeUTR <- invertStrand(threeUTR)
# Remove the part of anti_threeUTR that overlaps with a 5UTR/CDS/3UTR or TSS
removeOverlapRegion(anti_threeUTR, fiveUTR) %>%
    removeOverlapRegion(dm6cds) %>%
    removeOverlapRegion(threeUTR) %>%
    removeOverlapRegion(prombyTx) %>%
    export.bed(file.path(output_dir, "antisense_threeUTRs.bed"))

## get antisense 5 UTR
anti_fiveUTR <- invertStrand(fiveUTR)
# Remove the part of anti_fiveUTR that overlaps with a 5UTR/CDS/3UTR/TSS or anti-3primeUTR
removeOverlapRegion(anti_fiveUTR, fiveUTR) %>%
    removeOverlapRegion(dm6cds) %>%
    removeOverlapRegion(threeUTR) %>%
    removeOverlapRegion(prombyTx) %>%
    removeOverlapRegion(anti_threeUTR) %>%
    export.bed(file.path(output_dir, "antisense_fiveUTRs.bed"))

## get antisense to transcripts
transcriptsBy(gtf, "gene") %>%
    unlist() %>% unique() -> dm6transcripts
anti_transcript <- invertStrand(dm6transcripts)
# Remove the part of anti_transcript that overlaps with a transcript/anti-UTR on same strand
removeOverlapRegion(anti_transcript, dm6transcripts) %>%
    removeOverlapRegion(anti_threeUTR) %>%
    removeOverlapRegion(anti_fiveUTR) %>%
    export.bed(file.path(output_dir, "antisense_transcript.bed"))

## get antisense promoter
promoters <- transcriptsBy(gtf, "gene") %>%
    promoters(upstream = 200, downstream = 0) %>%
    unlist() %>% unique()
anti_promoter <- invertStrand(promoters)
# Remove the part of anti_promoter that overlaps with a transcript on either strand
removeOverlapRegion(anti_promoter, dm6transcripts, ignoreStrand = TRUE) %>%
    export.bed(file.path(output_dir, "antisense_promoters.bed"))
