__author__="jjenkins"
__date__ ="$Jul 20, 2010 11:02:39 PM$"

# python ../../src/buildChromosomes_after_JOIN_2.py ; 
# python ../../src/chromo_synteny_plot.py ; 
# python ../../src/leveragingSyntenicSequence.py

from hagsc_lib import isOverlapped
from hagsc_lib import iterBestHit
from hagsc_lib import FASTAFile_dict
from hagsc_lib import iterCounter

from os.path import join

import re

from sys import argv

#=====================================================================
def real_main():

    # Assembly File
    assemblyFile  = 'poplar14_5_snpFixed_phaseFixed.merged.broken.fasta'
    markerOutFile = 'merged.broken.bestHit'
    pre           = 'nisqually.synteny.broken'
    markerSet     = pre

    excludeSet = set([line[:-1] for line in open('excluded.Nisqually.dat')])
    singlePlacScaffs = set([line[:-1] for line in open('singlePlacScaffs_%s.out'%pre)])

    # Pulling the number of bases
    numBases = {}
    addOne = iterCounter(10000)
    fastaIndex = FASTAFile_dict( assemblyFile )
    for scaffID, record in fastaIndex.iteritems():
        numBases[record.id] = record.basePairs
        addOne()
    #####
    
    # Setting up the counters
    unMarkedScaffs = set(fastaIndex.keys())
    markedScaffs   = set()
    singleScaffs   = set()
    
    # Pulling marker placements
    recordDict = {}
    for record in iterBestHit( open(markerOutFile,'r') ):

        # Screening for quality
        if ( (record.per_ID < 98.0) or (record.per_coverage < 98.0) ): continue
        if ( record.BAC_recordName in excludeSet ): continue

        # Parsing the marker ID
        splitName = record.BAC_recordName.split('|')
        mrkrNum   = record.BAC_recordName
        if ( len(splitName) == 3 ):
            grpNum  = splitName[1]
            mapLoc    = float(splitName[0].split('-')[1]) / 1.0e6
        else:
            grpNum  = splitName[2]
            mapLoc    = float(splitName[3]) / 1.0e6
        #####

        superName  = record.scaffold
        scaffStart = record.scaffStart
        scaffEnd   = record.scaffEnd
        cov        = record.per_coverage
        ID         = record.per_ID
        
        # Gathering statistics
        if ( (superName in singlePlacScaffs) and (superName in unMarkedScaffs) ):
            unMarkedScaffs.remove(superName)
            singleScaffs.add(superName)
        elif ( superName in unMarkedScaffs ):
            unMarkedScaffs.remove(superName)
            markedScaffs.add(superName)
        #####
        
        # Skipping if single placement
        if ( superName in singlePlacScaffs ): continue

        # Building the output string
        outputString = []
        outputString.append( superName )
        outputString.append( '%d'%numBases[superName] )
        outputString.append( '%d'%record.scaffSize )
        outputString.append( '%d'%scaffStart )
        outputString.append( '%d'%scaffEnd )
        outputString.append( grpNum )
        outputString.append( record.placDir )
        outputString.append( ('%5.2f'%ID).strip() )
        outputString.append( ('%5.2f'%cov).strip() )
        outputString.append( mrkrNum )
        outputString.append( '%s'%mapLoc )
        outputString.append( '%s\n'%markerSet )
        finalTuple = ( int(scaffStart), '\t'.join( outputString ) )
        # Partitioning the information out
        if ( recordDict.has_key( grpNum ) ):
            if ( recordDict[grpNum].has_key(superName) ):
                recordDict[grpNum][superName].append( finalTuple )
            else:
                recordDict[grpNum][superName] = [finalTuple]
            #####
        else:
            recordDict[grpNum] = {}
            recordDict[grpNum][superName] = [finalTuple]
        #####
    #####
    
    # Deduping the gene and marker placements
    for grpNum, tmpDict in recordDict.iteritems():
        for superName in tmpDict.iterkeys():

            tmpList        = recordDict[grpNum][superName]
            tmpList.sort()
            nItems         = len(tmpList)
            removeIndicees = set()
            
            for n in xrange( nItems - 1 ):
                start_n, line_n = tmpList[n]
                s_n     = line_n.split(None)
                start_n = int( s_n[3] )
                end_n   = int( s_n[4] )
                for m in xrange( n+1, nItems ):
                    start_m, line_m = tmpList[m]
                    s_m     = line_m.split(None)
                    start_m = int( s_m[3] )
                    end_m   = int( s_m[4] )
                    if ( isOverlapped(start_n, end_n, start_m, end_m) ):
                        removeIndicees.add( m )
                    #####
                    if ( start_m > end_n ): break
                #####
            #####
            
            # Removing the offending indicees
            map( recordDict[grpNum][superName].pop, sorted(removeIndicees, reverse=True) )
            
        #####
    #####
    
    # Partitioning the reads
    for grpNum, superDict in recordDict.items():
        tmpHandle = open( '%s_Group_%s.dat'%(pre,grpNum), 'w' )
        DSU = [ (numBases[item],item) for item in superDict.keys()]
        DSU.sort()
        for superNum,superName in DSU:
            tupleList = superDict[superName]
            tupleList.sort()
            for tmpStart,tmpString in tupleList:
                # Git-er-dun
                tmpHandle.write( tmpString )
            #####
            tmpHandle.write( '--------------------\n' )
        #####
        tmpHandle.close()
    #####

    # Commafication
    commify = re.compile( r'(?<=\d)(?=(?:\d\d\d)+(?!\d))')
    
    # Counting the bases and scaffolds
    print '%d scaffolds with markers'%len(markedScaffs)
    print '%d scaffolds with no markers'%len(unMarkedScaffs)
    print '%d scaffolds are single placement'%len(singleScaffs)
    
    totalBases = 0
    for scaffID, record in fastaIndex.iteritems(): totalBases += len( record.seq )
    print 'Total bases: %s'%( commify.sub( ',', '%d'%(totalBases) ) )
    
    n_marked_Bases = 0
    for scaffID in markedScaffs: n_marked_Bases += len( fastaIndex[scaffID].seq )
    print '%s marked bases'%( commify.sub( ',', '%d'%(n_marked_Bases) ) )
    print '%.3f%% of total'%( 100.0 * float(n_marked_Bases) / float(totalBases) )
    
    n_unMarked_Bases = 0
    for scaffID in unMarkedScaffs: n_unMarked_Bases += len( fastaIndex[scaffID].seq )
    print '%s unmarked bases'%( commify.sub( ',', '%d'%(n_unMarked_Bases) ) )
    print '%.3f%% of total'%( 100.0 * float(n_unMarked_Bases) / float(totalBases) )

    n_singleScaff_Bases = 0
    for scaffID in singleScaffs: n_singleScaff_Bases += len( fastaIndex[scaffID].seq )
    print '%s single marker bases'%( commify.sub( ',', '%d'%(n_singleScaff_Bases) ) )
    print '%.3f%% of total'%( 100.0 * float(n_singleScaff_Bases) / float(totalBases) )

if ( __name__ == '__main__' ):
    real_main()
