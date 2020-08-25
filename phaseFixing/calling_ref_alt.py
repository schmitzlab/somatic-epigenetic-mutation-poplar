__author__="Jerry Jenkins, jjenkins@hudsonalpha.org"
__date__ ="7/12/18"

from hagsc_lib import FASTAFile_dict, iterCounter

from phaseFix_lib import readClass, cigarClass, outlierVariantClass

from sys import argv, stderr

#==============================================================
def evaluateSNP( tmpClass, readSeq, pos, ref, alt, readCigar, readBase ):
    
    # Evaluation of the reference base
    try:
        rstart = tmpClass.ref_to_read[ pos ]
    except KeyError:
        # This represents a SNP where the alignment is out of bounds for the read
        return "NO_REF_READ_POS_SNP"
    #####
    
    if ( rstart == "-" ):
        # This indicates that the read has an
        # error in it, and we can not make a call
        # of ref or alt using this read
        return "NO_REF_READ_POS_SNP"
    else:
        rend = rstart + len(ref)
    #####
    
    # Evaluation of the alternate base
    astart = tmpClass.ref_to_read[pos]
    if ( astart == "-" ):
        # This indicates that the read has an
        # error in it, and we can not make a call
        # of ref or alt using this read
        return "NO_REF_READ_POS_SNP"
    else:
        aend = astart + len(alt)
    #####
    
    if ( ref == readSeq[ rstart:rend ] ):
        if ( alt == readSeq[ astart:aend ] ):
            # Realisticaly, we should not ever make it here with a SNP call
            return "AMBIGUOUS_AND_EQUAL_SNP"
        else:
            return "REF"
        #####
    else:
        if ( alt == readSeq[ astart:aend ] ):
            return "ALT"
        else:
            # This is another example of where a read has an
            # error in it.
            return "AMBIGUOUS_AND_DIFFERENT_SNP"
        #####
    #####

#==============================================================
def evaluateINDEL( tmpClass, readSeq, pos, ref, alt, readCigar, readBase ):

    #-----------------
    # TYPING THE INDEL
    if ( len(alt) > len(ref) ):  # Example A/AG
        # INS = the read indicates more bases than the reference
        indelType = 'INS'
    elif ( len(ref) > len(alt) ):  # Example AG/A
        # DEL = the read indicates fewer bases than the reference
        indelType = 'DEL'
    else:
        print "SHOULD NEVER MAKE IT HERE"
        print ref, alt
        assert False
    #####
    
    # Pulling the read position based on the reference position
    try:
        readPos = tmpClass.ref_to_read[pos]
    except KeyError:
        # This is a siutation where the read is in error.
        return "MISCALLED_INS_INDEL"
    #####
    
    #-------------------------------
    # Evaluation of indel insertions
    if ( indelType == "INS" ):
    
        try:
            refRange = [tmpClass.read_to_ref[item] for item in xrange(readPos,readPos+len(alt))]
        except TypeError:
            # This is a siutation where the read is in error.
            return "MISCALLED_INS_INDEL"
        except KeyError:
            # This is a siutation where the read is in error.
            return "MISCALLED_INS_INDEL"
        #####
        
        # If the indel is an insertion, then we should
        # expect dashes in the reference alignment
        # That means it is an ALT call
        if ( '-' in set(refRange) ):
            # Does the ALT call match the read sequence?
            readEnd = readPos + len(alt)
            if ( alt == readSeq[readPos:readEnd] ):
                # If the read sequence matches the alternate, then it is a ALT call
                return 'ALT'
            else:
                # This is a siutation where the read is in error.
                return "MISCALLED_INS_INDEL"
            #####
        elif ( len(refRange) > 0 ): # Making sure we don't need to skip this step
            # Now check to see if the read sequence matches the reference
            readEnd = readPos + len(ref)
            if ( ref == readSeq[readPos:readEnd] ):
                # If the read sequence matches the reference, then it is a REF call
                return 'REF'
            else:
                # This is a siutation where an indel was called in error
                # Generally, this is a variant that should likely have been
                # called a SNP.
                return "MISCALLED_INS_INDEL"
            #####
        #####

    #-------------------------------
    # Evaluation of indel deletions
    elif ( indelType == "DEL" ):

        try:
            readRange = [tmpClass.ref_to_read[item] for item in xrange(pos,pos+len(ref))]
        except TypeError:
            # This is a siutation where the read is in error.
            print "Should not make it to here"
            assert False
            return "MISCALLED_DEL_INDEL"
        except KeyError:
            # This would be a case where the read alignment
            # did not extend through the ref bases.
            return "MISCALLED_DEL_INDEL"
        #####
        
        # If the indel is a deletion, then we should
        # expect dashes in the read alignment
        # That means it is an ALT call
        if ( '-' in set(readRange) ):
            # Is the indel miscalled already?
            try:
                readEnd = readPos + len(alt)
            except TypeError:
                # This is the case where the reference position maps to
                # a dash in the read alignment
                return "MISCALLED_DEL_INDEL"
            #####
            # Does the ALT call match the read sequence?
            if ( alt == readSeq[readPos:readEnd] ):
                # If the read alignment has dashes in it, and
                # the read sequence matches the alternate, 
                # then it is a ALT call
                return 'ALT'
            else:
                # This is a situation where the read base doesn't match
                # the read base
                return 'MISCALLED_DEL_INDEL'
            #####
        elif ( len(readRange) > 0 ):
            # Now check to see if the read sequence matches the reference
            readEnd = readPos + len(ref)
            if ( ref == readSeq[readPos:readEnd] ):
                # If the read sequence matches the reference, then it is a REF call
                return 'REF'
            else:
                # This is a siutation where an indel was called in error
                # Generally, this is a variant that should likely have been
                # called a SNP.
                return "MISCALLED_DEL_INDEL"
            #####
        #####
        
    ####

#==============================================================
def real_main():

    # Reading in the stage
    ARV_file    = argv[1]
    FASTA_file  = argv[2]
    P_outlier   = float( argv[3] )
    
    # Setting the output files
    PHASE_INFO   = ARV_file.replace( 'read_to_variant', 'phase_info' )
    OUTLIER_FILE = ARV_file.replace( 'read_to_variant', 'outlier' )

    #=====================================================
    # Pulling the scaffold from the FASTA file
    stderr.write( 'INDEXING THE ASSEMBLY\n' )
    indexFASTA   = FASTAFile_dict( FASTA_file )

    #----------------------
    # Read in the read sets
    stderr.write( 'READING IN THE READ/VARIANT ASSOCIATIONS\n' )
    read_to_indel   = {}
    readVariantDict = {}
    scaff_to_reads  = {}
    read_to_scaff   = {}
    var_to_scaff    = {}
    y               = iterCounter(100000)
    for line in open( ARV_file ):
    
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
            y()

        elif ( line[0] == "<" ):

            readBase, cigarString = line[1:].split(None)
            readVariantDict[readBase].addCigar( cigarString )

        else:

            # Filling the indel information
            readBase, indelID = line.split(None)
            readVariantDict[readBase].addIndel( indelID )
            
            var_to_scaff[indelID] = scaffID
            
        #####

    #####
    
    #-------------------------
    # Finding variants in reads
    stderr.write( 'Calling REF or ALT on the read\n' )
    z = iterCounter( 1000 )
    for scaffID, readList in scaff_to_reads.iteritems():
        
        chrSeq = str( indexFASTA[scaffID].seq )
        
        for readBase in readList:
            
            z()
            
            # Pulling the read sequence
            readSeq    = str(readVariantDict[readBase].seq)
            nReadBases = len(readSeq)
            readCigar  = readVariantDict[readBase].cigarString
            readStart  = readVariantDict[readBase].readStart
            readEnd    = readVariantDict[readBase].readEnd
            
            # Parse the cigarString
            tmpClass = cigarClass( readBase, readCigar, readStart, readSeq, chrSeq )
            
            # Looping over the relevant variants
            for indelID in readVariantDict[readBase].indelSet:
                
                # Splitting the variant ID
                pos, ref, alt = indelID.split(';')
                pos           = int(pos)
                
                # Are we working with a SNP or an INDEL
                if ( len(ref) * len(alt) == 1 ):
                    varType = "SNP"
                else:
                    varType = "INDEL"
                #####
                
                #-----------------------------------------
                # Making a call based on the variant type
                if ( varType == "SNP" ):
                    phaseCall = evaluateSNP( tmpClass, readSeq, pos, ref, alt, readCigar, readBase )
                elif ( varType == "INDEL" ):
                    phaseCall = evaluateINDEL( tmpClass, readSeq, pos, ref, alt, readCigar, readBase )
                #####
                
                # Adding the phase call
                readVariantDict[readBase].addPhaseCall( indelID, pos, phaseCall )
                
            #####
            
        #####
        
    #####
    
    #-------------------------
    # Write the phasing file
    stderr.write( 'WRITING PHASING INFORMATION\n' )
    oh = open( PHASE_INFO, 'w' )
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
    
    #-------------------------
    # Writing the outlier variants to a file
    oh = open( OUTLIER_FILE, 'w' )
    for outlierVariant, tmpClass in outlierVariantDict.iteritems():
        scaffID = var_to_scaff[outlierVariant]
        if ( tmpClass.outliers == 0 ): continue
        oh.write( '%s\t%s\t%d\t%d\n' % (scaffID, outlierVariant, tmpClass.outliers, tmpClass.nonOutliers ) )
    #####
    oh.close()

#==============================================================
if ( __name__ == '__main__' ):
    real_main()
