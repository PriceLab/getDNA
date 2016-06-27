library(getDNAClient)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_getSequenceByLocString()
   test_getSequenceByLocParameters()
   test_fromBedFile()
   
} # runTests    
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")
   gdc <- getDNAClient("hg38")

} # test_constructor
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
    tbl <- read.table(file, sep="\t", as.is=TRUE)


    colnames(tbl) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand",
                       "thickStart", "thickEnd", "itemRgb")

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
test_getSequences.2 <- function()
{
   printf("--- test_getSequences.2")
   gdc <- getDNAClient("hg19")
   
      # transcription start site of the TERT gene
   tss <- 1295262
   locs <- list(list(name="foo", chrom="chr5", start=tss-5, end=tss+5, revcomp=FALSE),
                list(name="bar", chrom="chr5", start=tss-1, end=tss+1, revcomp=FALSE))

   locs <- list(foo=list(chrom="chr5", start=tss-5, end=tss+5, revcomp=FALSE),
                bar=list(name="bar", chrom="chr5", start=tss-1, end=tss+1, revcomp=FALSE))

   seqs <- getDNA(gdc, locs)
   checkEquals(seqs$foo, "CGGGACGGGGC")
   checkEquals(seqs$bar, "ACG")

   seq.1.rv <- getDNA(gdc, list(name="bar", chrom="chr5", start=tss-1, end=tss+1, revcomp=FALSE))


} # test_getSequences.2
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
