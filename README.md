# getDNA
This service offers  easy, programmatic retrieval of (typically, relatively short) DNA sequences, by genone and chromosome location,
from Python and R.

At present (June 2016) this service uses the UCSC DAS (distributed annotation server).  A url of this form

```
http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr1:100000,100020
```

returns an XML document:

```
<DASDNA>
<SEQUENCE id="chr1" start="100000" stop="100020" version="1.00">
<DNA length="21">cactaagcacacagagaataa</DNA>
</SEQUENCE>
</DASDNA>
```

The Python and R clients provided here are simple convenience functions to facilitate the use of
this service from within scripts or programs.  The XML markup is stripped off, the
sequence is converted to upper case, reverse complementation performed if requested.

## Using the R client.  First, install the package and its prerequisites (XML, RCurl, Biostrings):

```
source("http://bioconductor.org/biocLite.R")
biocLite(c("XML", "RCurl", "Biostrings", "devtools"))
library(devtools)
install_github("PriceLab/getDNA", subdir="client-R/getDNAClient")
```

When the package and rerequisites are all installed, try it out:

```
library(getDNAClient)
dnaClient <- getDNAClient("hg19")
seq <- getSequenceByLocString(dnaClient, "chr5:1295260-1295270")
seq.rc <- getSequenceByLocString(dnaClient, "chr5:1295260-1295270", reverseComplement=TRUE)


```

Obtain seqence for regions specified in a bed file:

```
file <- system.file(package="getDNAClient", "extdata", "sample.bed")  # included with package
tbl <- read.table(file, sep="\t", as.is=TRUE, header=TRUE)
locs <- unlist(lapply(1:6, function(r) sprintf("%s:%d-%d", tbl[r,"chrom"], tbl[r, "chromStart"], tbl[r, "chromEnd"])))
   # give each segment a name if you wish
names(locs) <- head(tbl$name)
dnaClient <- getDNAClient("hg38")
seqs <- getSequencesByLocStrings(dnaClient, locs, reverseComplement=FALSE)
seqs.rc <- getSequencesByLocStrings(dnaClient, locs, reverseComplement=TRUE)

strand.status.rc <- tbl$strand == "-"    # creates a vector, one entry for each sequence
seqs.mixed <- getSequencesByLocStrings(dnaClient, locs, reverseComplement=strand.status.rc)
```

## Using the Python 3.x client.  First, install the module (more discussion needed here).


```
from GetDNAClient import *
print(DNAClientSupportedGenomes())

client = GetDNAClient("hg19")

   # implicit reverseComplement=False
client.getSequence("chr5", 1295262, 1295266) # CGGGG

   # explicit reverseComplement=False
client.getSequence("chr5", 1295262, 1295266, False) # CGGGG

   # explicit reverseComplement=True
client.getSequence("chr5", 1295262, 1295266, True) # CCCCG

```
