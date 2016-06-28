from GetDNAClient import *
from datetime import datetime
#------------------------------------------------------------------------------------------------------------------------
def testSequence():

   print("--- testSequence")
   client = GetDNAClient("hg19")
   assert(client.getSequence("chr5", 1295262, 1295266) == 'CGGGG')
   assert(client.getSequence("chr5", 1295262, 1295266, False) == 'CGGGG')
   
#------------------------------------------------------------------------------------------------------------------------
def testReverseComplementSequence():

   print("--- testReverseComplementSequence")
   client = GetDNAClient("hg19")
   assert(client.getSequence("chr5", 1295262, 1295266, True) == 'CCCCG')
   
#------------------------------------------------------------------------------------------------------------------------
def testAllGenomes():

   print("--- testSequence")
   for genome in DNAClientSupportedGenomes(): 
      client = GetDNAClient(genome)
      seq = client.getSequence("chr12", 57795963, 57796000)
      assert(len(seq) == 38)
   
#------------------------------------------------------------------------------------------------------------------------
testSequence()
testReverseComplementSequence()
testAllGenomes()