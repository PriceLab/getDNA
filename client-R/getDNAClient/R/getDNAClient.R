printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
.GetDNAClient <- setClass("GetDNAClientClass",
                          slots = c(url="character",
                                    genome="character",
                                    quiet="logical"
                                    )
                            )
#------------------------------------------------------------------------------------------------------------------------
setGeneric("getSequenceByLocString", signature="obj", function(obj, locString, reverseComplement=FALSE)
            standardGeneric("getSequenceByLocString"))
setGeneric("getSequenceByLoc", signature="obj", function(obj, chrom, start, end, reverseComplement=FALSE)
            standardGeneric("getSequenceByLoc"))
#------------------------------------------------------------------------------------------------------------------------
DNAClientSupportedGenomes <- function()
{
   c("hg16", "hg17", "hg18", "hg19", "hg38", "mm7", "mm8", "mm9", "mm10")

} # DNAClientSupportedGenomes
#------------------------------------------------------------------------------------------------------------------------
getDNAClient <- function(genome, quiet=FALSE)
{
   stopifnot(genome %in% DNAClientSupportedGenomes())

   url <- "http://genome.ucsc.edu/cgi-bin/das"
   .GetDNAClient(url=url, genome=genome, quiet=quiet)

} # getDNAClient, the constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("getSequenceByLocString", "GetDNAClientClass",

     function(obj, locString, reverseComplement=FALSE){
         s <- gsub(",", "", locString);
         stopifnot(grepl(":", s))
         stopifnot(grepl("-", s))
         tokens <- strsplit(s, ":")[[1]]
         chrom <- tokens[1]
         startEnd.tokens <- strsplit(tokens[2], "-")[[1]]
         start <- as.integer(startEnd.tokens[1])
         end <- as.integer(startEnd.tokens[2])
         getSequenceByLoc(obj, chrom, start, end, reverseComplement)
         })

#------------------------------------------------------------------------------------------------------------------------
setMethod("getSequenceByLoc", "GetDNAClientClass",

    function(obj, chrom, start, end, reverseComplement=FALSE){
       search.url <- sprintf("%s/%s/dna?segment=%s:%d,%d", obj@url, obj@genome,
                                                           chrom, start, end)
       text <- getURL(search.url)
       doc <- xmlParse(text)
       s0 <- xpathSApply(doc, "//DASDNA/SEQUENCE/DNA/text()")
       s1 <- xmlValue(s0[[1]])
       s2 <- gsub("\n", "", s1)
       if(reverseComplement)
          s2 <- toString(reverseComplement(DNAString(s2)))
       toupper(s2)
       }) # getSequenceByLoc

#------------------------------------------------------------------------------------------------------------------------
