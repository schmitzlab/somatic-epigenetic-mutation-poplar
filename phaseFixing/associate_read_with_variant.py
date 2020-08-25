__author__="Jerry Jenkins, jjenkins@hudsonalpha.org"
__date__ ="Created:  3/23/16"

from hagsc_lib import iterCounter, IntervalTree, IntervalClass
from hagsc_lib import FASTAFile_dict, revComp

from phaseFix_lib import compute_COV, compute_ApproximateMatch, vcfSplitter
from phaseFix_lib import cigarClass

from os.path import join, splitext

from sys import stderr, stdout, argv

import subprocess

import gzip

import re

#==============================================================
def real_main():
    
    # Reading in the stage
    bamFile     = argv[1]
    FASTA_file  = argv[2]
    minReadSize = int(argv[3])   # Minimum read size to be considered for haplotyping
    minCOV      = float(argv[4]) # Minimum coverage for accepting a read alignment
    vcf_FOFN    = argv[5]

    # Setting the debug state
    debug   = False

    #=====================================================
    # Initializing the inputs
    vcfFiles   = [ line[:-1] for line in open(vcf_FOFN) ]
    bamFiles   = [ bamFile ]
    
    #=====================================================
    # Initializing the outputs
    READ_VAR     = '%s.read_to_variant.dat' % splitext(bamFile)[0]
    
    #=====================================================
    # Pulling the scaffold from the FASTA file
    stderr.write( 'INDEXING THE ASSEMBLY\n' )
    indexFASTA   = FASTAFile_dict( FASTA_file )
    
    #=====================================================
    # PULL HETEROZYGOUS SNP AND INDEL CALLS FROM VCF FILES AND ASSOCIATE READS WITH VARIANTS

    #----------------------------------------------------------
    # Finding the best read alignment for each read
    stderr.write( "----------------------------\n" )
    stderr.write( "FINDING BEST READ ALIGNMENTS\n" )
    readSeqDict = {}  # key:readBase   value: Read sequence and associated cigar string
    bestHitDict = {}  # key: scaffID  value: List of alignment score, readBase, readStart, cigarString
    x           = iterCounter(1000000)
    testDict    = {True:None}
    for BAM_file in bamFiles:
        
        # Finding the bestHit for the read
        p              = subprocess.Popen( 'samtools view %s | cut -f1,2,3,4,6,10' % BAM_file, shell=True, stdout=subprocess.PIPE )
        tmpBestHitDict = {}
        for line in p.stdout:
            # Reading in a new line
            readBase, bitwiseFlag, scaffID, readStart, cigarString, readSeq = line.split(None)
            # Screening out the secondary alignments
            try:
                if ( bin(int(bitwiseFlag))[::-1][8] == "1" ): continue
            except IndexError:
                pass
            #####
            readStart = int(readStart) - 1 # This is shifted by 1 position because it is zero based
            nBases    = len(readSeq)
            # Screen read size
            if ( nBases < minReadSize ): continue
            # Screening for read coverage
            COV = compute_COV( cigarString, nBases )
            if ( COV < minCOV ): continue
            x()
            # Storing the read sequence in a dictionary
            readSeqDict[readBase] = readSeq
            # Computing total matches
            totalMatch = compute_ApproximateMatch( cigarString )
            # Storing the read sequence in a dictionary
            try:
                pm, ps, pcs, pid = tmpBestHitDict[readBase]
                try:
                    testDict[ totalMatch > pm ]
                    tmpBestHitDict[readBase] = ( totalMatch, readStart, cigarString, scaffID )
                except KeyError:
                    pass
                #####
            except KeyError:
                tmpBestHitDict[readBase] = ( totalMatch, readStart, cigarString, scaffID )
            #####
        #####
        p.poll()
        
        # Translating the read information to scaffolds
        for readBase, tmpTuple in tmpBestHitDict.iteritems():
            totalMatch, readStart, cigarString, scaffID = tmpTuple
            try:
                bestHitDict[scaffID].append( (totalMatch, readBase, readStart, cigarString) )
            except KeyError:
                bestHitDict[scaffID] = [ (totalMatch, readBase, readStart, cigarString) ]
            #####
        #####
        
        if ( debug ): break
        
    #####

    del tmpBestHitDict
    
    #------------------------
    # Reading in the VCF file
    stderr.write( "-------------------------------------------\n" )
    stderr.write( "READING IN THE HET CALLS FROM THE VCF FILES\n" )
    x       = iterCounter(100000)
    hetDict = {}  # Dictionary containing key:scaffID and value:(snpPos, ref, and alt bases)
    for vcfFile in vcfFiles:
        for line in open( vcfFile ):
            # Screening the comment lines
            if line[0] == '#': continue
            # Parsing the VCF line
            try:
                currentSNP = vcfSplitter(line)
                x()
            except ValueError:
                print line
                assert False
            #####
            # Shifting the current snp position
            snpPos = currentSNP.scaffPos - 1 # This is shifted by 1 position because it is zero based
            # This is the case where it is not 0/1, 1/1, or 1/2, but something entirely different
            if ( not currentSNP.goodSNP ): continue
            # Screening the depth ratio on the het
            if ( currentSNP.isHet ):
                if ( currentSNP.is_1_2 ): # The case of a consensus error
                    try:
                        hetDict[currentSNP.scaffID].append( (snpPos, currentSNP.refBase, currentSNP.altBase_1) )
                    except KeyError:
                        hetDict[currentSNP.scaffID] = [ (snpPos, currentSNP.refBase, currentSNP.altBase_1) ]
                    #####
                    hetDict[currentSNP.scaffID].append( (snpPos, currentSNP.refBase, currentSNP.altBase_2) )
                else: # The case of a true het
                    try:
                        hetDict[currentSNP.scaffID].append( (snpPos, currentSNP.refBase, currentSNP.altBase) )
                    except KeyError:
                        hetDict[currentSNP.scaffID] = [ (snpPos, currentSNP.refBase, currentSNP.altBase) ]
                    #####
                #####
            #####
        #####
    #####
    
    #----------------------------------------------------------
    # Associating variant calls with a readID
    stderr.write( "----------------------------------\n" )
    stderr.write( "ASSOCIATE VARIANT CALLS WITH READS\n" )
    read_to_variant_Dict = {}  # key:readID   value: List of indels
    readCIGAR_dict       = {}  # key:readID   value: List of (cigarString, readStart)
    read_to_scaff_map    = {}  # key:readID   value: scaffID
    x                    = iterCounter(1000000)
    for scaffID, tmpBestHitList in bestHitDict.iteritems():
        
        #-----------------
        # Build the interval tree
        try:
            hetDict[scaffID].sort()
        except KeyError:
            continue
        #####
        
        print scaffID
        
        indelList = []
        for pos_n, ref_n, alt_n in hetDict[scaffID]:
            pos_m = pos_n + 1
            indelList.append( IntervalClass( pos_n, pos_m, '%d;%s;%s'%(pos_n,ref_n,alt_n) ) )
        #####
        tmpTree = IntervalTree( indelList )
        
        #-----------------
        # Pulling the scaffold sequence
        scaffSeq = str( indexFASTA[scaffID].seq )
        
        #-----------------
        # Looping over the best hits for the scaffold
        for totalMatch, readBase, readStart, cigarString in tmpBestHitList:
            
            # Associating reads with scaffolds
            read_to_scaff_map[readBase] = scaffID
            
            # Pulling the read sequence
            readSeq = readSeqDict[readBase]
            
            # Parsing the cigar string
            tmpClass = cigarClass( readBase, cigarString, readStart, readSeq, scaffSeq )
            
            # Finding the actual end of the read on the reference using the cigar string
            valList = [ item for item in tmpClass.read_to_ref.itervalues() if item != '-' ]
            try:
                readEnd = sorted( valList )[-1]
            except IndexError:
                print scaffID
                print readBase
                print readStart
                print valList
                print tmpClass.displayAlignment()
                print "You should not be here!!"
                assert False
            #####

            # Storing the read cigar string
            readCIGAR_dict[readBase] = ( cigarString, readStart, readEnd, scaffID )
            
            # Find all variants that are associated with this read
            tmpVars = tmpTree.find( readStart, readEnd )
            for item in tmpVars:
                indelID = item.label
                try:
                    read_to_variant_Dict[readBase].append( indelID )
                except KeyError:
                    read_to_variant_Dict[readBase] = [ indelID ]
                #####
            #####
            
            x()
            
        #####
        
    #####
    
    #----------------------------------
    # Writing the information to a file
    print "-------------------------------------------------"
    print "WRITING THE READ TO VARIANT INFORMATION TO A FILE"
    oh = open( READ_VAR, 'w' )
    x           = iterCounter(1000000)
    for readBase, tmpList in read_to_variant_Dict.iteritems():
        cigarString, readStart, readEnd, scaffID = readCIGAR_dict[readBase]
        readStart = int(readStart)
        readEnd   = int(readEnd)
        readSeq   = readSeqDict[readBase]
        oh.write( ">%s\t%s\t%d\t%d\t%s\n"%( readBase, scaffID, readStart, readEnd, readSeq ) )
        oh.write( '<%s\t%s\n' % ( readBase, cigarString ) )
        for indelID in tmpList:
            oh.write( '%s\t%s\n'%(readBase, indelID) )
            x()
        #####
    #####
    oh.close()
    
#==============================================================
if ( __name__ == '__main__' ):
    real_main()
    