# getDNA
This service offers  easy, programmatic retrieval of (typically, relatively short) DNA sequences, by genone and chromosome location,
from Python and R.

At present (June 2016) this service uses the UCSC DAS (distributed annotation server).  An url of this form

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

The Python and R clients provided here are simply convenience functions, facilitating the use of
this service from within programs,  stripping off the XML markup, converting to upper case,
and returning the sequence's reverse complement if requested.

# Using the R client.  First, install the package and its prerequisites (XML, RCurl, Biostrings):

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
