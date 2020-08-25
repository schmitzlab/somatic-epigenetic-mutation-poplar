__author__="Jerry Jenkins, jjenkins@hudsonalpha.org"
__date__ ="Created:  2/4/15"

from hagsc_lib import iterFASTA, writeFASTA, revComp, SeqRecord, FASTAFile_dict

from sys import argv

#==============================================================
def real_main():
    
    oh = open( 'poplar14_5_snpFixed_phaseFixed.merged.broken.joined.fasta', 'w' )
    
    # Indexing the FASTA
    indexedFASTA = FASTAFile_dict( 'poplar14_5_snpFixed_phaseFixed.merged.broken.fasta' )
    
    nameList = [ 'Chr%s' % ( str(n).zfill(2) ) for n in xrange(1,20) ]
    
    scaffSet  = set()
    sepString = 10000 * 'N'
    n = 0
    for line in open( 'nisqually.synteny.broken_joins.refined.dat' ):
        joinList = [item.strip() for item in line.split(' 10000 ')]
        seqList = []
        for scaffNum in joinList:
            if ( scaffNum.count('rc') == 1 ):
                scaffID = scaffNum.replace('rc','')
                scaffSet.add( scaffID )
                seqList.append( revComp(str(indexedFASTA[scaffID].seq).upper()) )
            else:
                scaffID = scaffNum
                scaffSet.add( scaffID )
                seqList.append( str(indexedFASTA[scaffID].seq).upper() )
            #####
        #####
        r = SeqRecord( id=nameList[n], seq=sepString.join(seqList), description='' )
        writeFASTA( [r], oh )
        n += 1
    #####
    
    for scaffID, r in indexedFASTA.iteritems():
        if ( r.id in scaffSet ): continue
        writeFASTA( [r], oh )
    #####
    
    oh.close()

#==============================================================
if ( __name__ == '__main__' ):
    real_main()
