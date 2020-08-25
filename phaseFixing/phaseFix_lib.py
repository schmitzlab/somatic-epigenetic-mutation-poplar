__author__="Jerry Jenkins, jjenkins@hudsonalpha.org"
__date__ ="Created:  7/28/15"

# from hagsc_lib import isGzipFile, isBzipFile, commify, IntervalClass, IntervalTree
# from hagsc_lib import iterFASTA, histogramClass, fileHeader_Class
# 
# from gp_lib import jobFinished, submitJob, create_sh_file
# 
# from math import ceil, floor, log
# 
# from numpy import array
# 
# from os import system
# 
# import subprocess

from os.path import join, isfile

from sys import stdout, stderr

import re

#================================
# Cigar string parsing regular expression
# cigarParse  = re.compile( r'([0-9]+)([SHIP=XMDN])' ).findall
cigarParse    = re.compile( r'(\d+)(\w)' ).findall
parseSoftClip = re.compile( r'(\d+)S' ).findall
parseHardClip = re.compile( r'(\d+)H' ).findall
parseMatches  = re.compile( r'(\d+)M' ).findall

#==============================================================
# Default values
runPath    = '/global/dna/projectdirs/plant/geneAtlas/HAGSC_TOOLS/jwj_pythonScripts'

#==============================================================
# Setting default values
def_MIN_READ_SIZE = 1000
def_ALIGN_COV     = 95.0
def_P_OUTLIER     = 0.65
def_ALIGN_MEM     = '10G'
def_MAX_READ_COV  = None
def_ALIGN_ID      = 0.80

#==============================================================
# Possible job types
jobTypeSet = set( [ 'align_reads', \
                    'read_to_variant', \
                    'partition_groups', \
                    'call_bases', \
                    'make_fixing_file' ] )

#==============================================================
def writeAndEcho( tmpStr, cmdLogFile ):
    stdout.write( '%s\n'%tmpStr )
    oh = open( cmdLogFile, 'a' )
    oh.write( '%s\n'%tmpStr )
    oh.close()

#==============================================================
def echoAndExecute( tmpCmd, cmdLogFile, execute=True ):
    stdout.write( '%s\n'%tmpCmd )
    oh = open( cmdLogFile, 'a' )
    oh.write( '%s\n'%tmpCmd )
    oh.close()
    if ( execute ): system( tmpCmd )

#==============================================================
def throwError( tmpString ):
    stderr.write( '-------------------------\n' )
    stderr.write( 'ERROR:  %s\n'%tmpString )
    stderr.write( 'Halting execution\n' )
    stderr.write( '-------------------------\n' )
    assert False
    
#=====================================================================
def compute_COV( cigarString, nBases ):
    nSoftClipped = sum( map( int, parseSoftClip(cigarString) ) )
    nHardClipped = sum( map( int, parseHardClip(cigarString) ) )
    # Allowing for the hard clipped bases:  Remember that the nBases already has hard clipping removed
    # but to get a proper COV you have to add them on in the denominator
    return 100.0 * float(nBases - nSoftClipped) / ( float(nBases) + nHardClipped )
    # Previous equation, changed 5/9/2019
#     return 100.0 * float(nBases - nSoftClipped) / float(nBases)

#=====================================================================
def compute_ApproximateMatch( cigarString ):
    return sum( map( int, parseMatches(cigarString) ) )

#========================
def parseConfigFile( configFile, basePath ):
    
    # Checking for the config file
    if ( not isfile(configFile) ):
        throwError( 'parseConfigFile could not find config file %s'%configFile )
    #####

    # Initializing the key set
    keySet = set(['RUNID', \
                  'REFERENCE', \
                  'READS', \
                  'EMAIL', \
                  'NUM_READS', \
                  'NUM_JOBS', \
                  'EXCLUDE', \
                  'MIN_READ_SIZE', \
                  'ALIGN_COV', \
                  'VCF_FOFN', \
                  'P_OUTLIER', \
                  'ALIGN_MEM', \
                  'MAX_READ_COV', \
                  'ALIGN_ID' ])

    # Parsing the config file and pulling the key:value relationships
    configDict = {}
    for line in open( configFile ):
        
        # Breaking out the key:value pair
        key, value   = line.split(None)
        
        # Adding the key:value pair
        if ( key in set( [ 'NUM_READS', 'NUM_JOBS', 'MIN_READ_SIZE' ] ) ): 
            configDict[key] = int( value )
        elif ( key in set( [ 'ALIGN_COV', 'P_OUTLIER', 'MAX_READ_COV', 'ALIGN_ID' ] ) ):
            configDict[key] = float(value)
        else:
            configDict[key] = value
        #####
        
        print key, value
        
        # Chedking for any keys that are mislabeled
        try:
            keySet.remove(key)
        except KeyError:
            throwError( 'parseConfigFile encountered an unrecognized key: %s'%key )
        #####
        
    #####
    
    #-----------------------------------    
    # Creating the command log filename
    try:
        cmd_log = join( basePath, '%s_cmdLog.dat'%configDict['RUNID'] )
        if ( not isfile(cmd_log) ):
            # Clearing the command log file
            oh = open( cmd_log, 'w' )
            oh.close()
        #####
        # Writing to the command log
        writeAndEcho( '- EXECUTING: parseConfigFile', cmd_log )
    except KeyError:
        throwError( 'Unable to read RUNID' )
    #####

    #-----------------------------------    
    # Checking the ALIGN_ID
    try:
        x = configDict['ALIGN_ID']
    except KeyError:
        # Default P_OUTLIER
        configDict['ALIGN_ID'] = def_P_OUTLIER
        writeAndEcho( '- Using the default ALIGN_ID %.3f' % def_ALIGN_ID, cmd_log )
        keySet.remove('ALIGN_ID')
    #####

    #-----------------------------------    
    # Checking the P_OUTLIER
    try:
        x = configDict['P_OUTLIER']
    except KeyError:
        # Default P_OUTLIER
        configDict['P_OUTLIER'] = def_P_OUTLIER
        writeAndEcho( '- Using the default P_OUTLIER %.3f' % def_P_OUTLIER, cmd_log )
        keySet.remove('P_OUTLIER')
    #####

    #-----------------------------------    
    # Checking the exclusion file
    try:
        x = configDict['EXCLUDE']
    except KeyError:
        # Default exclude file is None
        configDict['EXCLUDE'] = None
        writeAndEcho( '- Using the default empty excluded reads file ', cmd_log )
        keySet.remove('EXCLUDE')
    #####

    #-----------------------------------    
    # Checking the MIN_READ_SIZE
    try:
        x = configDict['MIN_READ_SIZE']
    except KeyError:
        # Default value
        configDict['MIN_READ_SIZE'] = def_MIN_READ_SIZE
        writeAndEcho( '- Using the default MIN_READ_SIZE = %d '%def_MIN_READ_SIZE, cmd_log )
        keySet.remove('MIN_READ_SIZE')
    #####

    #-----------------------------------    
    # Checking the ALIGN_COV
    try:
        x = configDict['ALIGN_COV']
    except KeyError:
        # Default value
        configDict['ALIGN_COV'] = def_ALIGN_COV
        writeAndEcho( '- Using the default ALIGN_COV = %.1f '%def_ALIGN_COV, cmd_log )
        keySet.remove('ALIGN_COV')
    #####

    #-----------------------------------    
    # Checking the ALIGN_MEM
    try:
        x = configDict['ALIGN_MEM']
    except KeyError:
        # Default value
        configDict['ALIGN_MEM'] = def_ALIGN_MEM
        writeAndEcho( '- Using the default ALIGN_MEM = %s '%def_ALIGN_MEM, cmd_log )
        keySet.remove('ALIGN_MEM')
    #####

    #-----------------------------------    
    # Checking the MAX_READ_COV
    try:
        x = configDict['MAX_READ_COV']
    except KeyError:
        # Default value
        configDict['MAX_READ_COV'] = def_MAX_READ_COV
        writeAndEcho( '- Using the default MAX_READ_COV = %s '%def_MAX_READ_COV, cmd_log )
        keySet.remove('MAX_READ_COV')
    #####
    
    #-----------------------------------    
    # Checking for missing keys
    if ( len(keySet) > 0 ):
        strList = ['parseConfigFile detected that the following keys are missing in your phaseFix.config file\n']
        for key in keySet:
            strList.append( '\t%s\n'%key )
        #####
        strList.append( 'Please add these keys to your config file and resubmit the snp calling job.')
        throwError( ''.join(strList) )
    #####
    
    #-----------------------------------    
    # Checking to see if the reference exists
    if ( not isfile(configDict['REFERENCE']) ):
        throwError( 'REFERENCE genome can not be located %s'%(configDict['REFERENCE']) )
    #####

    #-----------------------------------    
    # Checking to see if the reads exist
    if ( not isfile(configDict['READS']) ):
        throwError( 'READS file can not be located %s'%(configDict['READS']) )
    #####

    return ( configDict, cmd_log )

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
class cigarClass( object ):

    def __init__( self, readBase, cigarString, readStart, readSeq, refSeq ):
        # Inputs
        self.cigarString  = cigarString
        self.readStart    = readStart
        self.readSeq      = readSeq
        self.nBases       = len(readSeq)
        self.readBase     = readBase
        
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
        
        # Removing the 5' and 3' hard clipped bases:  Added 5/9/2019
        if ( self.parsedString[0][1]  == "H" ): self.parsedString.pop(0)
        if ( self.parsedString[-1][1] == "H" ): self.parsedString.pop(-1)
        
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
                    print self.readBase
                    print self.cigarString
                    print N, alignType
                    print len(refSeq)
                    print len(self.readSeq)
                    print refPos, refPos+N
                    print self.align_readSeq
                    print readPos, readPos+N
                    print self.align_refSeq
                    print self.alignment
                    print "Should not be here"
                    assert False
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
        excludeSet = set( [ 'AMBIGUOUS_AND_DIFFERENT_SNP', \
                            'MISCALLED_DEL_INDEL', \
                            'MISCALLED_INS_INDEL', \
                            'NO_REF_READ_POS_SNP' ] )
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
def createGroupDictionary( primaryScaffsFile ):

    # Loading the group dictionary
    group_to_scaffList = {}
    scaff_to_group     = {}
    for line in open( primaryScaffsFile ):
        scaffID, grpNum = line.split( None )
        scaff_to_group[scaffID] = grpNum
        try:
            group_to_scaffList['GRP_%s'%grpNum].append( scaffID )
        except KeyError:
            group_to_scaffList['GRP_%s'%grpNum] = [ scaffID ]
        #####
    #####
    # Finding the final group
    DSU = [ (int(item.split('_')[1]),item) for item in group_to_scaffList.iterkeys() ]
    DSU.sort()

    # Return to the user
    return group_to_scaffList, scaff_to_group, DSU[-1][1]

