library(getDNAClient)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_supportedGenomes()
   test_getSequenceByLocString()
   test_getSequenceByLocParameters()
   test_fromBedFile()
   test_allGenomes()
   test_getSequencesByLocStrings()
   
} # runTests    
#------------------------------------------------------------------------------------------------------------------------
test_supportedGenomes <- function()
{
   printf("--- test_supportedGenomes")
   checkTrue(all(c("hg18", "hg19", "hg38", "mm8", "mm9", "mm10") %in% DNAClientSupportedGenomes()))   


} # test_supportedGenomes
#------------------------------------------------------------------------------------------------------------------------
test_getSequenceByLocString <- function()
{
    printf("--- test_getSequenceByLocString")

    gdc <- getDNAClient("hg19")
    getSequenceByLocString(gdc, "chr5:1295260-1295262")

    start <- 1295262
    end.loc <- unlist(lapply(1:3,  function(i) start + (3*i)))
    locs <- sprintf("chr5:%d-%d", 1295262, end.loc)

    seqs <- unlist(lapply(locs, function(loc) getSequenceByLocString(gdc, loc)))
    checkEquals(seqs, c("CGGG", "CGGGGCG", "CGGGGCGGGG"))

    seqs <- unlist(lapply(locs, function(loc) getSequenceByLocString(gdc, loc, reverseComplement=TRUE)))
    checkEquals(seqs, c("CCCG", "CGCCCCG", "CCCCGCCCCG"))


} # test_getSequenceByLocString
#------------------------------------------------------------------------------------------------------------------------
test_getSequenceByLocParameters <- function()
{
    printf("--- test_getSequenceByLocParameters")

    gdc <- getDNAClient("hg19")

    chrom <- "chr5"
    start <- 1295262
    end.loc <- unlist(lapply(1:3,  function(i) start + (3*i)))

    seq.1 <- getSequenceByLoc(gdc, chrom, start, end.loc[1])
    seq.1.rc <- getSequenceByLoc(gdc, chrom, start, end.loc[1], reverseComplement=TRUE)

    checkEquals(seq.1, "CGGG")
    checkEquals(seq.1.rc, "CCCG")

    seq.2 <- getSequenceByLoc(gdc, chrom, start, end.loc[2])
    seq.2.rc <- getSequenceByLoc(gdc, chrom, start, end.loc[2], reverseComplement=TRUE)
    checkEquals(seq.2, "CGGGGCG")
    checkEquals(seq.2.rc, "CGCCCCG")

    seq.3 <- getSequenceByLoc(gdc, chrom, start, end.loc[3])
    seq.3.rc <- getSequenceByLoc(gdc, chrom, start, end.loc[3], reverseComplement=TRUE)
    checkEquals(seq.3, "CGGGGCGGGG")
    checkEquals(seq.3.rc, "CCCCGCCCCG")

} # test_getSequenceByLocString
#------------------------------------------------------------------------------------------------------------------------
test_getSequencesByLocStrings <- function()
{
   printf("--- test_getSequencesByLocStrings")

   file <- system.file(package="getDNAClient", "extdata", "sample.bed")
   checkTrue(file.exists(file))
   tbl <- read.table(file, sep="\t", as.is=TRUE, header=TRUE)

   locs <- unlist(lapply(1:6, function(r) sprintf("%s:%d-%d", tbl[r,"chrom"], tbl[r, "chromStart"], tbl[r, "chromEnd"])))

   gdc <- getDNAClient("hg38")
   seqs <- getSequencesByLocStrings(gdc, locs)
   checkEquals(seqs, c("CTTACAGT", "TATGGAATGT", "CATTGCTC", "CTCTAGGA", "TGAGTTAAAT", "AATGAGAA"))
   seqs.rc <- getSequencesByLocStrings(gdc, locs, TRUE)
   checkEquals(seqs.rc, c("ACTGTAAG", "ACATTCCATA", "GAGCAATG", "TCCTAGAG", "ATTTAACTCA", "TTCTCATT"))

   seqs.mixed <- getSequencesByLocStrings(gdc, locs, c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE))
   checkEquals(seqs.mixed, c("ACTGTAAG", "TATGGAATGT", "GAGCAATG", "CTCTAGGA", "ATTTAACTCA", "AATGAGAA"))

   names(locs) <- head(tbl$name)
   seqs.names <- getSequencesByLocStrings(gdc, locs)
   checkEquals(names(seqs.names), names(locs))   
   
} # test_getSequencesByLocString
#------------------------------------------------------------------------------------------------------------------------
test_getSequenceByLoc <- function()
{
   gdc <- getDNAClient("hg19")
   tert.tss <- 1295262
   seq <- getSequenceByLoc(gdc, chrom="chr5", start=tert.tss-2, end=tert.tss+2)
   checkEquals(seq, "GACGG")
   seq.revcomp <- getSequenceByLoc(gdc, chrom="chr5", start=tert.tss-2, end=tert.tss+2, reverseComplement=TRUE)
   checkEquals(seq.revcomp, "CCGTC")

} # getSequence
#------------------------------------------------------------------------------------------------------------------------
test_fromBedFile <- function()
{
    printf("--- test_fromBedFile")
    file <- system.file(package="getDNAClient", "extdata", "sample.bed")
    checkTrue(file.exists(file))
    tbl <- read.table(file, sep="\t", as.is=TRUE, header=TRUE)

      # these segments are all 1167 bsaes long - which is more than we need, and makes
      # our results harder to check.  artificially shorten them here

    set.seed(31)
    tbl$chromEnd <- tbl$chromStart + as.integer(runif(9, 5, 10))

    locStrings <- sprintf("%s:%d-%d", tbl$chrom, tbl$chromStart, tbl$chromEnd)
    names(locStrings) <- tbl$name
    strands <- tbl$strand
    revComp <- strands == "-"
    
    gdc <- getDNAClient("hg38")
    seqs <- lapply(1:length(locStrings), function(i) getSequenceByLocString(gdc, locStrings[i], revComp[i]))
    names(seqs) <- names(locStrings)
    checkTrue(all(nchar(seqs) == 1 + tbl$chromEnd - tbl$chromStart))

} # test_fromBedFile
#------------------------------------------------------------------------------------------------------------------------
test_allGenomes <- function()
{
   printf("--- test_allGenomes")
   locString <- "chr12:57,795,963-57,796,000"
   for(genome in DNAClientSupportedGenomes()){
      client <- getDNAClient(genome); 
      seq <- getSequenceByLocString(client, locString)
      checkEquals(nchar(seq), 38)
      }

}  # test_allGenomes
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
