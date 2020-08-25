__author__="Jerry Jenkins, jjenkins@hudsonalpha.org"
__date__ ="Created:  10/17/13"

from hagsc_lib import  QUALFile_dict, commify
from hagsc_lib import FASTAFile_dict

from hagsc_lib import SeqRecord, QUALRecord, writeFASTA, writeQUAL

from sys import stderr

import subprocess

import re

seqParser = re.compile( r'([ACTG]+)\(([ACTG]+)/([ACTG]+)\)([ACTG]+)' ).finditer

#==============================================================
def real_main():
    
    # Inputs to the process
    fastaBase       = '/projectb/scratch/jjenkins/Paspalum/phaseFixing/Assembly/combined.paspalum'
    phaseFixingFile = '/projectb/scratch/jjenkins/Paspalum/phaseFixing/phaseFixing/PASP.HetPollution.Heterozygous.context.dat'
    
    # Indexing assem
    indexAssemblyFASTA = FASTAFile_dict( '%s.fasta' % fastaBase )
    
    # Opening the output FASTA file
    oh_f = open( '%s.hetFixed.fasta' % fastaBase, 'w' )
    
    # Log file
    oh_log = open( 'hetPollution_fixing_log.dat', 'w' )
    
    # Reading in the snps from the snp calling file
    stderr.write( 'Reading in the snps\n' )
    oh_log.write( "-Reading in snps\n" )
    snp_indel_Dict = {}
    
    # Reading in the indels from the snp calling file
    stderr.write( 'Reading in the indels\n' )
    oh_log.write( "-Reading in indels\n" )
    n         = 1
    for line in open(phaseFixingFile):
        n += 1
        try:
            scaffID, pos, SNP, tmpType   = line.split(None)
        except ValueError:
            continue
        #####
        try:
            leadSeq, ref, alt, tailSeq = seqParser(SNP).next().groups()
        except StopIteration:
            continue
        #####
        refSeq = ''.join([leadSeq,ref,tailSeq])
        altSeq = ''.join([leadSeq,alt,tailSeq])
        pos    = int(pos)
        try:
            snp_indel_Dict[scaffID].append( (pos, refSeq, altSeq, ref, alt, tmpType) )
        except KeyError:
            snp_indel_Dict[scaffID] =     [ (pos, refSeq, altSeq, ref, alt, tmpType) ]
        #####
    #####
    
    # Reading in the hetAnalysis
    stderr.write( 'Reading in the full het fixing file\n' )
    
    # Removing SNPs and indels that are too close to one another
    stderr.write( 'Removing snps and indels that are too close to one another\n' )
    modifiedSet = set()
    for scaffID, tmpList in snp_indel_Dict.iteritems():
        
        # Storing the modified supers
        modifiedSet.add( scaffID )
        
        # Testing for overlap of snps and indels
        stderr.write( 'TESTING SNP AND INDEL OVERLAP %s\n'%scaffID )
        tmpList.sort()
        removalSet  = set()
        removeItems = False
        tmpSeq      = str( indexAssemblyFASTA[scaffID].seq ).upper()
        for n in xrange( len(tmpList) - 1 ):
            # Reading in a pair
            pos_n, refSeq_n, altSeq_n, ref_n, alt_n, type_n = tmpList[n]
            for m in xrange( (n+1), len(tmpList) ):
                # Reading in a pair
                pos_m, refSeq_m, altSeq_m, ref_m, alt_m, type_m = tmpList[m]
                # If there is an overlap, then act
                if ( ((pos_n+len(alt_n)) >= pos_m) or ((pos_n+len(ref_n)) >= pos_m) ):
                    if ( type_n == 'RAW' ):
                        removeItems = True
                        removalSet.add(n)
                    elif ( type_m == 'RAW' ):
                        removeItems = True
                        removalSet.add(m)
                    else:
                        # This is the case where both are clones.
                        # In this case, I randomly remove one of them
                        removeItems = True
                        removalSet.add(n)
                    #####
                #####
                # Breaking after we are done
                if ( ( pos_m > (pos_n+len(alt_n)) ) and (pos_m > (pos_n+len(ref_n)) ) ): break
            #####
        #####
        stderr.write( 'done\n')
        
        # Removing items that are too close together
        if ( removeItems ):
            for tmpIndex in sorted( removalSet, reverse=True ):
                pos, refSeq, altSeq, ref, alt, itemType = tmpList[tmpIndex]
                oh_log.write( '--------------------\n' )
                oh_log.write( 'REMOVED\t%s\t%s\n'%(itemType,scaffID) )
                oh_log.write( '%s\t%s\n'%(scaffID,refSeq) )
                oh_log.write( '%s\t%s\n'%(scaffID,altSeq) )
                oh_log.write( '%s\t%s\t%s\n'%(itemType,ref,alt) )
                oh_log.write( 'Position of the %s:  %s\n'%(itemType,commify(pos)) )
                oh_log.write( 'Reference Bases:  %s\n'%tmpSeq[pos-len(ref):pos] )
                tmpList.pop(tmpIndex)
            #####
        #####
        
        # Making the fixes on the remaining items
        tmpList.sort( reverse=True )
        for pos, refSeq, altSeq, ref, alt, itemType in tmpList:
            # Writing to the log file
            oh_log.write( '--------------------\n' )
            oh_log.write( 'MODIFIED\t%s\t%s\n'%(itemType,scaffID) )
            oh_log.write( '%s\t%s\n'%(scaffID,refSeq) )
            oh_log.write( '%s\t%s\n'%(scaffID,altSeq) )
            oh_log.write( '%s\t%s\t%s\n'%(itemType,ref,alt) )
            oh_log.write( 'Position of the %s:  %s\n'%(itemType,commify(pos)) )
            oh_log.write( 'Prior to Change Bases:  %s\n'%tmpSeq[(pos-50):(pos+50)] )
            # Making a check on the reference sequence
            if ( tmpSeq[(pos - 1):(pos + len(ref) - 1)] != ref ):
                print "WE HAVE A PROBLEM"
                print pos
                print refSeq
                print altSeq
                print ref, alt, itemType
                print tmpSeq[(pos - 50):(pos + 50)]
                print tmpSeq[pos:(pos + len(ref))]
                print ref
                assert False
            #####
            # Modifying the sequence
            frontSeq = tmpSeq[:pos-1]
            rearSeq  = tmpSeq[pos + (len(ref) - 1):]
            tmpSeq   = ''.join( [frontSeq, alt, rearSeq] )
            oh_log.write( 'Post Change Bases:      %s\n'%tmpSeq[(pos-50):(pos+50)] )
        #####
        
        # Writing the fixed SEQ/QUAL to a file
        fr = SeqRecord(  id=scaffID, seq=tmpSeq,   description='' )
        writeFASTA( [fr], oh_f )
        
    #####

    print 'Writing the remaining FASTA'
    for scaffID, record in indexAssemblyFASTA.iteritems():
        if ( scaffID in modifiedSet ): continue
        oh_log.write( 'Writing remaining FASTA:  %s\n'%record.id )
        writeFASTA( [record], oh_f )
    #####
    
    oh_f.close()
    oh_log.write( '--------------------\n' )
    
#==============================================================
if ( __name__ == '__main__' ):
    real_main()
