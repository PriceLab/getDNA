from lxml import etree
import sys
#------------------------------------------------------------------------------------------------------------------------
class GetDNAClient:

   'Access to reference DNA sequence on a (remote) server'

   baseURL = None
   genome = None

   def __init__(self, genome):
      if(not genome in ["hg18", "hg19", "hg38"]):
         sys.exit()
      self.genome = genome
      self.baseURL = 'http://genome.ucsc.edu/cgi-bin/das/%s/dna?segment=' % genome

   def getSequence(self, chrom, start, end, reverseComplement=False):
     url = self.baseURL + chrom + ':' + str(start) + ',' + str(end)
     doc = etree.parse(url)
     if doc == '':
        'THE SEQUENCE DOES NOT EXIST FOR GIVEN COORDINATES: %s:%d-%d' % (chrom, start, end)
     sequence = doc.xpath('SEQUENCE/DNA/text()')[0].replace('\n','').upper()
     if(reverseComplement):
        sequence = self.revComp(sequence)
     return(sequence)

   def revComp(self, sequence):
     complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
     result = "".join(complement.get(base, base) for base in reversed(sequence))
     return(result)

