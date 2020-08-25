__author__="Jerry Jenkins, jjenkins@hudsonalpha.org"
__date__ ="Created:  5/24/17"

from hagsc_lib import findTelomereCentromere, iterFASTA

from sys import argv

#==============================================================
def real_main():
    
    assemblyFASTA = argv[1]
    stage         = argv[2]
    
    if ( stage == '1' ):
        # Looking for telomere/centrmere sequence
        patterns  = [ 'TTAGGG', 'TTTAGGG', 'TTTCGGG' ]
        findTelomereCentromere( assemblyFASTA, patterns, repLow=10, repHigh=22 )
    #####
    
    if ( stage == '2' ):
        telomereSet = set( [line[:-1] for line in open('possibleTelomereContainingScaffs.dat')] )
        for r in iterFASTA(open(assemblyFASTA)):
            if ( r.id in telomereSet ):
                tmpSeq = str(r.seq)
                print r.id
                print tmpSeq[:2000]
                print "----------"
                print tmpSeq[-2000:]
                print "======================================"
            #####
        #####
    #####

    if ( False ):
        for r in iterFASTA(open('remaining.fasta')):
            if ( r.id == 'tig00035141' ):
                tmpSeq = str(r.seq)
                print r.id
                print tmpSeq[31282-100:31360+100]
                print "======================================"
            #####196750	196784
        #####
    #####

#==============================================================
if ( __name__ == '__main__' ):
    real_main()
