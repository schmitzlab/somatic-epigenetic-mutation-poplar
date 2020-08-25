__author__="Jerry Jenkins, jjenkins@hudsonalpha.org"
__date__ ="Created:  10/23/15"

from hagsc_lib import iterFASTA, FASTAFile_dict, SeqRecord, writeFASTA

#==============================================================
def real_main():

    breakDict = {}
    for line in open('nisqually.synteny_new_break_file.refined.dat'):
        scaffID, d1, d2, d3, breakPoint = line.split(None)
        try:
            breakDict[scaffID].append( ( int(breakPoint), int(breakPoint) ) )
        except KeyError:
            breakDict[scaffID] = [ ( int(breakPoint), int(breakPoint) ) ]
        #####
    #####
    
    oh = open( 'poplar14_5_snpFixed_phaseFixed.merged.broken.fasta', 'w' )
    for r in iterFASTA( open( 'poplar14_5_snpFixed_phaseFixed.merged.fasta' ) ):
        try: 
            # Sorting the breaks
            tmpList = breakDict[r.id]
            tmpList.sort(reverse=True)
            # Pulling the sequence
            tmpSeq = str(r.seq)
            # Setting the number of breaks
            n = len(tmpList) + 1
            # Making the breaks
            for start, end in tmpList:
                leftSeq  = tmpSeq[:start]
                rightSeq = tmpSeq[end:]
                print "%s_%d"%(r.id,n)
                tmp_r    = SeqRecord( id="%s_%d"%(r.id,n), seq=rightSeq, description='' )
                writeFASTA( [tmp_r], oh )
                n -= 1
                tmpSeq   = leftSeq
            #####
            # Writing the final record
                print "%s_%d"%(r.id,n)
            tmp_r    = SeqRecord( id="%s_%d"%(r.id,n), seq=tmpSeq, description='' )
            writeFASTA( [tmp_r], oh )
        except KeyError:
            # If no break, then write the record
            print r.id
            writeFASTA( [r], oh )
        #####
    #####
    oh.close()
    
#==============================================================
if ( __name__ == '__main__' ):
    real_main()
