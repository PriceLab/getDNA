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
testSequence()
testReverseComplementSequence()