__author__="Jerry Jenkins, jjenkins@hudsonalpha.org"
__date__ ="Created:  3/23/16"

from hagsc_lib import iterCounter, IntervalTree, IntervalClass
from hagsc_lib import FASTAFile_dict, revComp

from os.path import join, splitext

from sys import stderr, stdout, argv

import subprocess

import gzip

import re

#================================
# Cigar string parsing regular expression
# cigarParse  = re.compile( r'([0-9]+)([SHIP=XMDN])' ).findall
cigarParse    = re.compile( r'(\d+)(\w)' ).findall
parseSoftClip = re.compile( r'(\d+)S' ).findall
parseMatches  = re.compile( r'(\d+)M' ).findall

#==============================================================
class vcfSplitter:
    def __init__(self, line):
        # Pulling the initial information
        self.s        = line.split(None)
        self.scaffID  = self.s[0]
        self.scaffPos = int( self.s[1] )
        self.regSet   = set(["0/1", "1/1"])
        
        # Building the depth dictionary
        GT_value = self.s[9].split(':')[0]
	    
        # preComputing call depth
        self.goodSNP = True
        self.is_1_2  = False
        if ( GT_value in self.regSet ):
            # Computing zygosity
            # GT:  0/1 is a heterozygous and 1/1 homozygous
            self.isHet = ( GT_value == "0/1" )
            # Generating the ref and alt bases
            self.refBase  = self.s[3]
            self.altBase  = self.s[4]
            # Screening ref/alt bases on size
            if ( (len(self.refBase) > 3) or (len(self.altBase) > 3) ): self.goodSNP = False
        elif ( GT_value == "1/2" ):
            #================================
            # Computing zygosity
            # Strictly speaking this is a het call (GT:1/2), but it is actually just
            # an error in the consensus that needs to be fixed
            self.isHet  = True
            self.is_1_2 = True
            # Setting the type
            self.refBase    = self.s[3]
            splitBase       = self.s[4].split(',')
            self.altBase_1  = splitBase[0]
            self.altBase_2  = splitBase[1]
            # Screening ref/alt bases on size
            if ( (len(self.refBase) > 3) or \
                 (len(self.altBase_1) > 3) or \
                 (len(self.altBase_2) > 3) ): self.goodSNP = False
        else:
            self.goodSNP = False
        #####
	    
#=====================================================================
def compute_COV( cigarString, nBases ):
    nSoftClipped = sum( map( int, parseSoftClip(cigarString) ) )
    return 100.0 * float(nBases - nSoftClipped) / float(nBases)

#=====================================================================
def compute_ApproximateMatch( cigarString ):
    return sum( map( int, parseMatches(cigarString) ) )

#=====================================================================
class cigarClass( object ):

    def __init__( self, cigarString, readStart, readSeq, refSeq ):
        # Inputs
        self.cigarString  = cigarString
        self.readStart    = readStart
        self.readSeq      = readSeq
        self.nBases       = len(readSeq)
        
        # Calculated values
        self.align_refSeq  = []
        self.alignment     = []
        self.align_readSeq = []
        self.ref_to_read   = {} # key:refPos  value:readPos
        self.read_to_ref   = {} # key:readPos value: refPos
        self.totalMismatch = 0
        self.totalMatch    = 0
        self.nSoftClipped  = 0
    
        # Initial parse of string
        self.parsedString = cigarParse( cigarString )
        
        # Computing everything else
        self.parseCigarString( refSeq )
        
    def parseCigarString(self, refSeq):
        # Initializing the variables
        refPos   = self.readStart
        alignPos = 0
        readPos  = 0
        
        # Parsing the string
        for N, alignType in self.parsedString:

            N = int(N)
            
            if ( alignType == "M" ): # CONSUMES BOTH QUERY AND TARGET BASES
                
                # Determining the Match/Mismatch
                try:
                    t1       =       refSeq[  refPos : refPos  + N ]
                    t2       = self.readSeq[ readPos : readPos + N ]
                    tmpAlign = ''.join( [ "|" if t1[m] == t2[m] else "X" for m in xrange(N) ] )
                except IndexError:
                    break
                #####

                #----------------------------
                # Performing the mapping
                for m in xrange( N ):
                    self.ref_to_read[refPos  + m] = readPos + m
                    self.read_to_ref[readPos + m] = refPos  + m
                #####
                
                # Incrementing the reference
                self.align_refSeq.append( refSeq[ refPos : refPos + N ] )
                refPos += N
                
                # Incrementing the alignment
                self.alignment.append( str(tmpAlign) )
                alignPos += N
    
                # Incrementing the read
                self.align_readSeq.append( self.readSeq[ readPos : readPos + N ] )
                readPos += N
    
            elif ( alignType in "X=SHIP" ): # CONSUMES QUERY BASES
                
                # Performing the mapping
                for m in xrange(N): self.read_to_ref[ readPos + m ] = '-'
                
                # Incrementing the relative positions
                self.align_refSeq.append( N * "-" )
                
                self.alignment.append( N * " " )
                alignPos += N
    
                self.align_readSeq.append( self.readSeq[readPos:readPos+N] )
                readPos += N
                
                # Computing the soft clipping
                self.nSoftClipped += N if (alignType == "S") else 0
    
            elif ( alignType in "DN" ): # CONSUMES TARGET BASES
    
                # Performing the mapping
                for m in xrange(N): self.ref_to_read[ refPos + m ] = '-'
    
                # Incrementing the relative positions
                self.align_refSeq.append( refSeq[ refPos : refPos + N ] )
                refPos += N
                
                self.alignment.append( N * " " )
                alignPos += N
    
                self.align_readSeq.append( N * "-" )
    
            #####
            
        #####
        
        self.totalMismatch = ''.join(self.alignment).count('X')
        self.totalMatch    = ''.join(self.alignment).count('|')
        
    def displayAlignment( self, wrap=80 ):
        newRef   = ''.join( self.align_refSeq )
        newAlign = ''.join( self.alignment )
        newRead  = ''.join( self.align_readSeq )
        for n in xrange( len(newRef) / wrap + 1 ):
            start = n       * wrap
            end   = min( len(newRef), (n + 1) * wrap )
            print "READ  %d   %s   %d" % ( start, newRead [ start : end ], end - 1 )
            print "ALIGN %d   %s   %d" % ( start, newAlign[ start : end ], end - 1 )
            print "REF   %d   %s   %d" % ( start, newRef  [ start : end ], end - 1 )
            print " "
        #####
        
#==============================================================
class readClass(object):
    
    def __init__( self, readBase, readSeq, readStart, readEnd ):
        self.readBase    = readBase
        self.seq         = readSeq
        self.readStart   = readStart
        self.readEnd     = readEnd
        self.cigarString = None
        self.indelSet    = set()
        self.indelDict   = {}
    
    def addCigar( self, cigarString ):
        self.cigarString = cigarString
    
    def addIndel( self, indelID ):
        self.indelSet.add( indelID )
    
    def addPhaseCall( self, indelID, pos, phaseCall ):
        # Storing the final set of information
        self.indelDict[indelID] = ( pos, phaseCall )
    
    def writeAndEcho( self, oh, outputString ):
#         stdout.write( '%s\n' % outputString )
        oh.write( '%s\n' % outputString )
    
    def showPhasing( self, oh, scaffID ):
        self.writeAndEcho( oh, "------------------------------" )
        self.writeAndEcho( oh, self.readBase )
        DSU = [ ( tmpTuple[0], indelID, tmpTuple[1] ) for indelID, tmpTuple in self.indelDict.iteritems() ]
        DSU.sort()
        for pos, indelID, outcome in DSU:
            self.writeAndEcho( oh, '%s\t%s\t%s' % ( scaffID, indelID, outcome ) )
        #####
    
    def findOutliers( self, P_outlier ):
        # Performing the voting
        excludeSet = set( [ 'AMBIGUOUS_AND_DIFFERENT', \
                            'AMBIGUOUS_AND_EQUAL', \
                            'NO_REF_READ_POS' ] )
        voteDict   = {'REF':0, 'ALT':0}
        nTotal     = 0
        for indelID, tmpTuple in self.indelDict.iteritems():
            pos, outcome = tmpTuple
            if ( outcome in excludeSet ): continue
            try:
                voteDict[outcome] += 1
            except KeyError:
                continue
            #####
            nTotal += 1
        #####
        DSU    = [ (v,k) for k,v in voteDict.iteritems() ]
        DSU.sort()
        
        # Declaring a winner
        winningVote, winner = DSU[-1]
        
        # Computing the win ratio
        try:
            r = float( winningVote ) / float( nTotal )
        except ZeroDivisionError:
            r = 1.0
        #####
        
        # Find all variants that are outliers
        outlierVariantSet    = set()
        nonOutlierVariantSet = set()
        if ( r >= P_outlier ):
            for indelID, tmpTuple in self.indelDict.iteritems():
                pos, outcome = tmpTuple
                if ( outcome in excludeSet ): continue
                if ( outcome != winner ):
                    outlierVariantSet.add( indelID )
                else:
                    nonOutlierVariantSet.add( indelID )
                #####
            #####
        #####
        return ( outlierVariantSet, nonOutlierVariantSet )
        
#==============================================================
class outlierVariantClass( object ):

    def __init__( self, varID ):
        self.varID       = varID
        self.outliers    = 0
        self.nonOutliers = 0
    
    def addOutlier(self):
        self.outliers += 1

    def addNonOutliers(self):
        self.nonOutliers += 1
        
#==============================================================
def real_main():
    
    # Reading in the stage
    stage    = int( argv[1] )
    
    # Reading in the bam file for stage 1
    if ( stage == 1 ):
        bamFile  = argv[2]
        bamFiles = [ bamFile ]
    #####
    
    # Setting the debug state
    debug   = False
    
    P_outlier   = 0.65 # Reads have to have this fraction of one call to be counted for outliers
    minReadSize = 1000 # Minimum read size to be considered for haplotyping
    minCOV      = 95.0 # Minimum coverage for accepting a read alignment
    
    #=====================================================
    # Main Variables
    basePath = '/projectb/scratch/jjenkins/Poplar_14_5_phaseFixing'
    
    #=====================================================
    # Initializing the inputs
    FASTA_file = join( basePath, 'combined.snpFixed.fasta' )
    
    vcf_FOFN   = join( basePath, 'vcf.fofn' )
    vcfFiles   = [ line[:-1] for line in open(vcf_FOFN) ]
    
    R2V_FOFN   = join( basePath, 'read_to_variant.fofn' )
    R2V_files  = [ line[:-1] for line in open(R2V_FOFN) ]

    #=====================================================
    # Initializing the outputs
    if ( stage == 1 ):
        READ_VAR     = '%s.read_to_variant.snpFixed.dat' % splitext(bamFile)[0]
    #####
    FULL_FIXING  = 'full_fixing.snpFixed.dat'
    OUTLIER_FILE = 'outlierResults.snpFixed.dat'
    
    #=====================================================
    # Pulling the scaffold from the FASTA file
    stderr.write( 'INDEXING THE ASSEMBLY\n' )
    indexFASTA   = FASTAFile_dict( FASTA_file )
    
    #=====================================================
    # STEP 1:  PULL HETEROZYGOUS SNP AND INDEL CALLS FROM VCF FILES AND ASSOCIATE READS WITH VARIANTS
    if ( stage == 1 ):

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
            p              = subprocess.Popen( 'samtools view %s | cut -f1,3,4,6,10' % BAM_file, shell=True, stdout=subprocess.PIPE )
            tmpBestHitDict = {}
            for line in p.stdout:
                # Reading in a new line
                readBase, scaffID, readStart, cigarString, readSeq = line.split(None)
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
        # Associating indel calls with a readID
        print "--------------------------------"
        print "ASSOCIATE INDEL CALLS WITH READS"
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
                tmpClass = cigarClass( cigarString, readStart, readSeq, scaffSeq )
                
                # Finding the actual end of the read on the reference using the cigar string
                valList = [ item for item in tmpClass.read_to_ref.itervalues() if item != '-' ]
                readEnd = sorted( valList )[-1]
                
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
        
    #####

    #========================================================
    # STAGE 2:  Separate out the
    if ( stage == 3 ):

    
    #========================================================
    # STEP 2:  Screen the read pairs for ones that overlap a 
    # heterozygous snp where the read matches the reference.
    if ( stage == 3 ):
        
        #----------------------
        # Read in the read sets
        stderr.write( 'READING IN THE READ/VARIANT ASSOCIATIONS\n' )
        read_to_indel   = {}
        readVariantDict = {}
        scaff_to_reads  = {}
        read_to_scaff   = {}
        y               = iterCounter(1000000)
        for READ_VAR in R2V_files:
            print READ_VAR
            for line in open( READ_VAR ):
                y()
                # Pulling the read sequence
                if ( line[0] == ">" ):
                    readBase, scaffID, readStart, readEnd, readSeq = line[1:].split(None)
                    readStart = int(readStart)
                    readEnd   = int(readEnd)
                    readVariantDict[readBase] = readClass( readBase, readSeq, readStart, readEnd )
                    # Mapping a scaffold to a set of reads
                    try:
                        scaff_to_reads[scaffID].add( readBase )
                    except KeyError:
                        scaff_to_reads[scaffID] = set( [ readBase ] )
                    #####
                    # Mapping a read to a scaffold
                    read_to_scaff[ readBase ] = scaffID
                elif ( line[0] == "<" ):
                    readBase, cigarString = line[1:].split(None)
                    readVariantDict[readBase].addCigar( cigarString )
                else:
                    # Filling the indel information
                    readBase, indelID = line.split(None)
                    readVariantDict[readBase].addIndel( indelID )
                #####
            #####
        #####
        
        #-------------------------
        # Finding variants in reads
        stderr.write( 'Calling REF or ALT on the read\n')
        stderr.write( 'Number of scaffolds %d\n' % len(scaff_to_reads.keys()) )
        z = iterCounter( 1000000 )
        # Looping over the indel ID's in the chromosome
        nn = 0
        for scaffID, readList in scaff_to_reads.iteritems():
            
            nn += 1
            print scaffID, nn
            
            chrSeq = str( indexFASTA[scaffID].seq )
            
            for readBase in readList:
                
                # Pulling the read sequence
                readSeq    = str(readVariantDict[readBase].seq)
                nReadBases = len(readSeq)
                readCigar  = readVariantDict[readBase].cigarString
                readStart  = readVariantDict[readBase].readStart
                readEnd    = readVariantDict[readBase].readEnd
                
                # Parse the cigarString
                tmpClass = cigarClass( readCigar, readStart, readSeq, chrSeq )
                
                # Looping over the relevant variants
                for indelID in readVariantDict[readBase].indelSet:
                    
                    # Splitting the variant ID
                    pos, ref, alt = indelID.split(';')
                    pos           = int(pos)
                    
                    # Setting the print variable
                    printAlign = True
                    
                    # Evaluation of the reference base
                    try:
                        rstart = tmpClass.ref_to_read[ pos ]
                        rend   = rstart + len(ref)
                    except TypeError:
                        phaseCall = "NO_REF_READ_POS"
                        printAlign = False
                    except KeyError:
                        phaseCall = "NO_REF_READ_POS"
                        printAlign = False
                    #####
                    
                    # Evaluation of the alternate base
                    try:
                        astart = tmpClass.ref_to_read[pos]
                        aend   = astart + len(alt)
                    except TypeError:
                        phaseCall = "NO_REF_READ_POS"
                        printAlign = False
                    except KeyError:
                        phaseCall = "NO_REF_READ_POS"
                        printAlign = False
                    #####
                    printAlign = True if printAlign else printAlign
                    
                    # Printing the calling
                    if ( printAlign ):
                        if ( ref == readSeq[ rstart:rend ] ):
                            if ( alt == readSeq[ astart:aend ] ):
                                phaseCall = "AMBIGUOUS_AND_EQUAL"
                            else:
                                phaseCall = "REF"
                            #####
                        else:
                            if ( alt == readSeq[ astart:aend ] ):
                                phaseCall = "ALT"
                            else:
                                phaseCall = "AMBIGUOUS_AND_DIFFERENT"
                            #####
                        #####
                    else:
                        pass
                    #####
                    
                    # Adding the phase call
                    readVariantDict[readBase].addPhaseCall( indelID, pos, phaseCall )
                    
                #####
                
            #####
            
        #####
        
        stderr.write( 'WRITING PHASING INFORMATION\n' )
        oh = open( FULL_FIXING, 'w' )
        for readBase, tmpClass in readVariantDict.iteritems():
            scaffID = read_to_scaff[readBase]
            tmpClass.showPhasing( oh, scaffID )
        #####
        oh.close()

        #-------------------------
        # Find outlier variants
        stderr.write( 'WRITING OUTLIER FILE\n' )
        outlierVariantDict = {}
        for readBase, tmpClass in readVariantDict.iteritems():
            
            scaffID = read_to_scaff[readBase]
            
            # Pulling outliers
            outlierVariantSet, nonOutlierVariantSet = tmpClass.findOutliers( P_outlier )
            
            # Adding outliers
            for varID in outlierVariantSet:
                try:
                    outlierVariantDict[varID].addOutlier()
                except KeyError:
                    outlierVariantDict[varID] = outlierVariantClass( varID )
                #####
            #####

            # Counting non-outliers
            for varID in nonOutlierVariantSet:
                try:
                    outlierVariantDict[varID].addNonOutliers()
                except KeyError:
                    outlierVariantDict[varID] = outlierVariantClass( varID )
                #####
            #####

        #####
        
        oh = open( OUTLIER_FILE, 'w' )
        for outlierVariant, tmpClass in outlierVariantDict.iteritems():
            oh.write( '%s\t%s\t%d\t%d\n' % (scaffID, outlierVariant, tmpClass.outliers, tmpClass.nonOutliers ) )
        #####
        oh.close()
    
    #####
    
#==============================================================
def profile_main():
    from cProfile import Profile
    from pstats import Stats
    prof  = Profile().runctx("real_main()", globals(), locals())
    stats = Stats( prof ).sort_stats("time").print_stats(60)
    return

#==============================================================
if ( __name__ == '__main__' ):
    real_main()
#     profile_main()
