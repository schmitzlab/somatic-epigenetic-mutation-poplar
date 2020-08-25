__author__="Jerry Jenkins, jjenkins@hudsonalpha.org"
__date__ ="$${date} ${time}$"

from hagsc_lib import isOverlapped, IntervalTree, IntervalClass

from os.path import join

from sys import argv, stdin

import re

import subprocess

#==============================================================
def writeHeader( breakType, oh, line, prev_line, intervalDict, oh_break ):
    
    # Determining whether there is a break that overlaps
    LG_1, mapPos_1, scaffID_1, size_1, start_1, end_1, strand_1, markerID_1 = line.split(None)
    LG_2, mapPos_2, scaffID_2, size_2, start_2, end_2, strand_2, markerID_2 = prev_line.split(None)
    
    x        = sorted(map(int,[start_1,end_1,start_2,end_2]))
    start    = x[1]
    stop     = x[2]
    scaffINT = scaffID_1
    try:
        OB = intervalDict[scaffID_1].find( start, stop )
    except KeyError:
        OB = []
    #####
    
    # Writing the header
    oh.write( 'BREAK_TYPE: %s\n'%breakType )
    oh.write( '%s'%prev_line )
    oh.write( '%s'%line )
    if ( OB == [] ):
        oh.write( 'NO_SCAFFOLD_GAPS\n' )
        oh_break.write( '%s %d %d  **\n'%(scaffINT,start,stop) )
    else:
        for tmpClass in OB:
            oh.write( 'SCAFFOLD_BREAK:\t%d\t%d\n'%(tmpClass.start,tmpClass.stop) )
        #####
        # Finding the cloesest one
        midpoint = ( start + stop )/ 2
        DSU = [ ( abs((item.stop+item.start)/2-midpoint),item) for item in OB ]
        DSU.sort()
        winningClass = DSU[0][1]
        oh_break.write( '%s %d %d\n'%(scaffINT,winningClass.start,winningClass.stop) )
    #####
    oh.write( '-----------------------------------------------------\n')
    return

#==============================================================
def real_main():
    
    # Pulling inputs
    pre          = argv[1]
    jumpDist     = float( argv[2] )
    minMarkers   = int( argv[3] )
    intervalDict = {}
    
    # Loading in the BAC and FOS files
    parsedMarkerFile = 'sortByScaff_%s.out'%pre
    
    # Opening the output file
    oh = open( '%s_new_breaks.dat'%pre, 'w' )
    oh_break = open( '%s_new_break_file.dat'%pre, 'w' )
    
    # Loading in the parsed marker file
    prev_LG         = None
    prev_pos        = None
    prev_scaffStart = None
    prev_line       = None
    mainList        = [[]]
    for line in open( parsedMarkerFile ):
    
        # Did we hit a new line?
        if ( line[0] == '=' ):
            
            # Looking for breaks
            if ( mainList != [[]] ):
                
                # Pulling out the lists that are too short
                for index in sorted([n for n in xrange(len(mainList)) if (len(mainList[n])<minMarkers)], reverse=True):
                    mainList.pop( index )
                #####

                # Has the linkage group changed?
                for n in xrange( len(mainList) - 1 ):
                
                    LG_n,   pos_n,   scaffID_n, line_n   = mainList[n][-1]
                    LG_np1, pos_np1, scaffID_m, line_np1 = mainList[n+1][0]
                    
                    if ( LG_n != LG_np1 ):
                        writeHeader( 'LG_CHANGE', oh, line_np1, line_n, intervalDict, oh_break )
                    #####

                    if ( (LG_n==LG_np1) and (abs( pos_n - pos_np1 ) > jumpDist) ):
                        writeHeader( 'MAP_POS_JUMP', oh, line_np1, line_n, intervalDict, oh_break )
                    #####
                    
                #####
                
            #####
            
            prev_LG         = None
            prev_scaffStart = None
            prev_line       = None
            prev_pos        = None
            mainList        = [ [] ]
            continue

        #####
        
        # Parse the line
        splitLine  = line.split(None)
        LG         = splitLine[0].replace('**','')
        
        mapPos     = float( splitLine[1] )
        scaffStart = int(   splitLine[4] )
        scaffID    = splitLine[2]
        
        # Is this the first line in a set?
        if ( prev_LG == None ):
            prev_LG         = LG
            prev_scaffStart = scaffStart
            prev_line       = line
            prev_pos        = mapPos
            mainList[-1].append( (LG,mapPos,scaffID,line) )
            continue
        #####
        
        # Looking for breaks
        if ( LG != prev_LG ):  mainList.append( [] )

        # Looking for breaks
        if ( (LG == prev_LG) and (abs(prev_pos-mapPos) > jumpDist) ):  mainList.append( [] )
        
        # Adding the previous line
        try:
            mainList[-1].append( (LG,mapPos,scaffID,line) )
        except IndexError:
            pass
        #####
        
        # Updating the previous variables
        prev_LG         = LG
        prev_scaffStart = scaffStart
        prev_line       = line
        prev_pos        = mapPos
        
    #####
    
    oh.close()
    
    oh_break.close()

#==============================================================
if ( __name__ == '__main__' ):
    real_main()
