__author__="Jerry Jenkins, jjenkins@hudsonalpha.org"
__date__ ="6/15/18"

from phaseFix_lib import createGroupDictionary

from hagsc_lib import iterFASTA

from sys import argv, stderr

from os.path import join

import re

seqParser = re.compile( r'([ACTG]+)\(([ACTG]+)/([ACTG]+)\)([ACTG]+)' ).finditer

#==============================================================
def real_main():

    # Genome FASTA
    genomeFASTA       = argv[1]
    minCov            = float(argv[2])
    maxCov            = float(argv[3])
    primaryScaffsFile = argv[4]
    basePath          = argv[5]
    RUNID             = argv[6]

    # Header Text
    headerText = 'Scaffold\tPosition\tSequence Context\tZygosity\n'
    
    # Reading in the primary scaffs file
    group_to_scaffList, scaff_to_group, finalGroup_ID = createGroupDictionary( primaryScaffsFile )
    
    #==================================================================
    # Reading in all of the single clone alignment vcf files
    stderr.write(  "\t-Reading outlier files\n" )
    finalFixingDict = {}
    for GRP in group_to_scaffList.iterkeys():
        
        stderr.write( "Analyzing Group: %s\n" % GRP )
        
        # Making the group variant file
        outlier_file = join( basePath, '%s.%s.outlier.dat' % ( RUNID, GRP ) )
        for line in open( outlier_file ):

            # Reading in the line
            chrID, snpInfo, d_outlier, d_nonOutlier = line.split(None)
            
            # Screening total coverage
            cov = int(d_outlier) + int(d_nonOutlier)
            if ( (cov < minCov) or (cov > maxCov) ): continue
            
            # Screening the coverage ratio
            try:
                r = float(d_outlier) / float(d_nonOutlier)
            except ZeroDivisionError:
                r = 2.0
            #####
            if ( r < 3.0 ): continue
            
            # Splitting the snp information
            snpPos, ref, alt = snpInfo.split(';')
            snpPos = int(snpPos)
            
            # Screening the depth ratio on the het
            try:
                finalFixingDict[chrID].append( (snpPos, ref, alt, "RAW") )
            except KeyError:
                finalFixingDict[chrID] =     [ (snpPos, ref, alt, "RAW") ]
            #####
            
        #####
        
    #####
    
    #==================================================================
    # Writing the output file
    stderr.write(  "\t-Writing fixing file\n" )
    flankSize = 100 # How much sequence is on either side of the flanking region
    het_outFile = join( basePath, '%s.HetPollution.Heterozygous.context.dat'%RUNID )
    oh_het      = open( het_outFile, 'w' )
    oh_het.write( headerText )
    for r in iterFASTA( open( genomeFASTA ) ):
        print r.id
        scaffID = r.id
        tmpSeq  = str( r.seq )
        nSize   = len(tmpSeq)
        #==========================================
        # Writing the true hets
        try:
            finalFixingDict[scaffID].sort()
            for snpPos, refBase, altBase, tmpType in finalFixingDict[scaffID]:
                # Computing the start position
                refLen   = len(refBase)
                start    = max( 0, (snpPos - flankSize) )
                end      = min( (snpPos + refLen + flankSize), nSize )
                leftSeq  = tmpSeq[           start : snpPos ]
                rightSeq = tmpSeq[ (snpPos+refLen) : end    ]
                if ( refBase != tmpSeq[snpPos:snpPos+len(refBase)].upper() ):
                    print 'something went wrong'
                    print r.id
                    print refBase, tmpSeq[snpPos:snpPos+len(refBase)], snpPos
                    print tmpSeq[snpPos-10:snpPos+10]
                    assert False
                #####
                oh_het.write( '%s\t%d\t%s(%s/%s)%s\t%s\n'%( scaffID, (snpPos+1), leftSeq.upper(), refBase, altBase, rightSeq.upper(), tmpType ) )
            #####
        except KeyError:
            continue
        #####
    #####
    oh_het.close()
    
#==============================================================
if ( __name__ == '__main__' ):
    real_main()
    