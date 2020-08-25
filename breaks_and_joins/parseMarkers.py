__author__="jjenkins"
__date__ ="3-Feb-2011"

from hagsc_lib import histogramClass, isOverlapped
from hagsc_lib import iterBestHit
from hagsc_lib import iterCounter
from hagsc_lib import iterFASTA
from hagsc_lib import commify

from math import log10

from sys import argv

#==============================================================
def real_main():
    
    # Parameters
    bestHitFile = 'merged.broken.bestHit'
    fastaFile   = 'poplar14_5_snpFixed_phaseFixed.merged.broken.fasta'
    preFix      = 'nisqually.synteny.broken'
    
    excludeSet = set([line[:-1] for line in open('excluded.Nisqually.dat')])
    
    # Screening the marker placements
    groupDict    = {}
    outputString = 8*['']
    for record in iterBestHit( open(bestHitFile, 'r') ):
        # Screening for quality
        if ( (record.per_ID < 98.0) or (record.per_coverage < 98.0) ): continue
        if ( record.BAC_recordName in excludeSet ): continue

        # Parsing the marker ID
        splitName = record.BAC_recordName.split('|')
        markerID  = record.BAC_recordName
        if ( len(splitName) == 3 ):
            mapGroup  = splitName[1]
            mapPos    = float(splitName[0].split('-')[1]) / 1.0e6
        else:
            mapGroup  = splitName[2]
            mapPos    = float(splitName[3]) / 1.0e6
        #####

        # Generating the output string
        outputString[0] = mapGroup
        outputString[1] = '%.4f'%mapPos
        outputString[2] = record.scaffold
        outputString[3] = '%d'%record.scaffSize
        outputString[4] = '%d'%record.scaffStart
        outputString[5] = '%d'%record.scaffEnd
        outputString[6] = record.placDir
        outputString[7] = '%s\n'%markerID
        finalString = '\t'.join(outputString)
        # Storing markers
        try:
            groupDict[mapGroup].append( (mapPos,record.scaffold,record.scaffSize,record.scaffStart,mapGroup,finalString) )
        except KeyError:
            groupDict[mapGroup] = [(mapPos,record.scaffold,record.scaffSize,record.scaffStart,mapGroup,finalString)]
        #####
    #####
    
    # Deduping the gene and marker placements
    for mapGroup in groupDict.iterkeys():
        tmpList = groupDict[mapGroup]
        tmpList.sort()
        nItems  = len(tmpList)
        removeIndicees = set()
        for n in xrange( nItems - 1 ):
            mapPos_n, scaff_n, size_n, start_n, group_n, line_n = tmpList[n]
            s_n     = line_n.split(None)
            start_n = int( s_n[4] )
            end_n   = int( s_n[5] )
            for m in xrange( n+1, nItems ):
                mapPos_m, scaff_m, size_m, start_m, group_m, line_m = tmpList[m]
                s_m     = line_m.split(None)
                start_m = int( s_m[4] )
                end_m   = int( s_m[5] )
                if ( isOverlapped(start_n, end_n, start_m, end_m) ):
                    removeIndicees.add( m )
#                     print line_n,
#                     print line_m,
#                     print "========================="
                #####
                if ( start_m > end_n ): break
            #####
        #####
        
        # Removing the offending indicees
        map( groupDict[mapGroup].pop, sorted(removeIndicees, reverse=True) )
        
    #####
    
    # Writing the results to an output file sorted by map and scaffold
    outputHandle_Map   = open( 'sortByMap_%s.out'%preFix, 'w' )
    outputHandle_Scaff = open( 'sortByScaff_%s.out'%preFix, 'w' )
    
    # Sorting by map
    DSU = [item for item in groupDict.keys()]
    DSU.sort()
    scaffSizes  = {}
    scaffPlacs  = {}
    for mapGroup in DSU:
        mapList = groupDict[mapGroup]
        mapList.sort()
        outputHandle_Map.write( '===========================================\n')
        for mapPos, scaffID, scaffSize, scaffStart, mapGroup, finalString in mapList:
            # Storing the scaffold size
            scaffSizes[scaffID] = scaffSize
            # Storing the scaffold
            try:
                scaffPlacs[scaffID].append( (scaffStart,mapGroup,finalString) )
            except KeyError:
                scaffPlacs[scaffID] = [(scaffStart,mapGroup,finalString)]
            #####
            outputHandle_Map.write(finalString)
        #####
    #####
    
    # Sorting by scaffold            
    DSU_scaff = [ (value,key) for key,value in scaffSizes.iteritems() ]
    DSU_scaff.sort(reverse=True)
    for scaffSize, scaffID in DSU_scaff:
        # Pulling the placement list
        placList = scaffPlacs[scaffID]
        placList.sort()
        # Looking for the most abundant map group
        tmpCmp = {}
        for scaffStart,mapGroup,finalString in placList:
            try:
                tmpCmp[mapGroup] += 1
            except KeyError:
                tmpCmp[mapGroup] = 1
            #####
        #####
        DSU_cmp = [(value,key) for key, value in tmpCmp.iteritems()]
        DSU_cmp.sort(reverse=True)
        targetGrp = DSU_cmp[0][1]
        # Write the list
        outputHandle_Scaff.write( '===========================================\n')
        for scaffStart, mapGroup, finalString in placList:
            if ( mapGroup != targetGrp ):
                outputHandle_Scaff.write( ''.join( ['**',finalString] ) )
            else:
                outputHandle_Scaff.write( finalString )
            #####
        #####
    #####
    
    # Determining the percentage of the assembly covered by the markers
    scaffSet = set( [scaffID for scaffSize, scaffID in DSU_scaff] )
    totalBases  = 0
    mappedBases = 0
    x = iterCounter(1000)
    scaffBases_dict = {}
    unmappedScaffSizes = []
    mappedScaffSizes   = []
    for record in iterFASTA( open( fastaFile, 'r' ) ):
        tmpSeq = str(record.seq).upper()
        nA = tmpSeq.count('A')
        nC = tmpSeq.count('C')
        nG = tmpSeq.count('G')
        nT = tmpSeq.count('T')
        scaffBases = nA + nC + nG + nT
        totalBases += scaffBases
        scaffBases_dict[record.id] = scaffBases
        if ( record.id in scaffSet ): 
            mappedBases += scaffBases
            mappedScaffSizes.append(log10(scaffBases))
        else:
            unmappedScaffSizes.append(log10(scaffBases))
        #####
        x()
    #####

    # Generating histogram of unmapped scaffolds
    x = histogramClass( min(unmappedScaffSizes), max(unmappedScaffSizes), 100 )
    map( x.addData, unmappedScaffSizes )
    print x.generateOutputString()

    # Generating histogram of mapped scaffolds 
    x = histogramClass( min(mappedScaffSizes), max(mappedScaffSizes), 100 )
    map( x.addData, mappedScaffSizes )
    print x.generateOutputString()
    
    print preFix
    print '%d total bases'%totalBases
    print '%d mapped bases'%mappedBases
    tmpPer = 100.0 * float(mappedBases) / float(totalBases)
    print '%.3f%% of total bases'%tmpPer
    print '======================='
    
    # Generating scaffold placement statistics
    orientable       = [0,0]
    unorientable     = [0,0]
    g10kbScaffs      = [0,0]
    singlePlacScaffs = set()
    for scaffID, placList in scaffPlacs.iteritems():
        if ( len(placList) == 1 ):
            unorientable[0] += 1
            unorientable[1] += scaffBases_dict[scaffID]
            singlePlacScaffs.add( scaffID )
        else:
            orientable[0] += 1
            orientable[1] += scaffBases_dict[scaffID]
        #####
        if ( scaffBases_dict[scaffID] >= 10000 ):
            g10kbScaffs[0] += 1
            g10kbScaffs[1] += scaffBases_dict[scaffID]
        #####
    #####
    
    print '%d scaffolds containing %d bases are orientable   (>=1 marker)'%(orientable[0],orientable[1])
    print '%d scaffolds containing %d bases are unorientable (1 marker)'%(unorientable[0],unorientable[1])
    print '%d scaffolds containing %d bases are >10kb        (>= 1 marker)'%(g10kbScaffs[0],g10kbScaffs[1])
    print "Single Placement Scaffolds:"
    for item in singlePlacScaffs:
        print "\t-%s\t%s" % ( item, commify(scaffBases_dict[item]) )
    #####
    
    # Writing the single placement scaffolds to a file
    oh = open( 'singlePlacScaffs_%s.out'%preFix, 'w' )
    for item in singlePlacScaffs:
        oh.write( '%s\n'%item )
    #####
    oh.close()
        
    # Closing the output files
    outputHandle_Map.close()
    outputHandle_Scaff.close()

#==============================================================
if ( __name__ == '__main__' ):
    real_main()
