__author__="jjenkins"
__date__ ="$Mar 1, 2011 2:13:00 PM$"

from hagsc_lib import FASTAFile_dict
from hagsc_lib import iterBestHit
from hagsc_lib import SeqRecord
from hagsc_lib import writeFASTA
from hagsc_lib import revComp

from os.path import split, splitext, isfile

from math import sqrt

from sys import stdout, argv

import re

#====================================================
# For Xenopus I ran the commands in this order:
#
# run splitBrokenMarkerPlacements.py to generate the initial Group files
#
# python ../src/buildingChromosomes.py ; python ../src/chromo_synteny_plot.py ; python ../src/leveragingSyntenicSequence.py
#====================================================

#=========================================================================
def pearsonsCorrelation( X, Y ):
    
    try:
        # Computing mean values
        mean_X = float( sum(X) ) / float( len(X) )
        mean_Y = float( sum(Y) ) / float( len(Y) )
    except ZeroDivisionError:
        print X, Y
        assert False
    #####

    # Computing STD
    sigma_X  = 0.0
    sigma_Y  = 0.0
    sigma_XY = 0.0
    for n in xrange(len(X)):
        dx = X[n] - mean_X
        dy = Y[n] - mean_Y
        sigma_X  += dx * dx
        sigma_Y  += dy * dy
        sigma_XY += dx * dy
    #####

    # Returning values
    if ( (sigma_X == 0.0) or (sigma_Y == 0.0) ):
        return 0.0
    #####
    return ( sigma_XY / (sqrt(sigma_X) * sqrt(sigma_Y) ) )

#####

#=========================================================================
class markerPlac_class( object ):

    def __init__(self, parsedLine):
        self.scaffID    = parsedLine[0]
        scaffStart      = int( parsedLine[3] )
        self.scaffStart = scaffStart
        scaffLength     = int( parsedLine[2] )
        self.revStart   = scaffLength - scaffStart + 1
        self.markerPos  = float( parsedLine[10] )

    def outputString(self):
        return '%s\t%d\t%d\t%.3f\n'%( self.scaffID, self.scaffStart, \
                                                 self.revStart, self.markerPos )
#####

#=========================================================================
class clusterClass(object):
	def __init__(self,initialPoint,initialIndex):
		self.points	 = [ (initialPoint,initialIndex) ]
		self.sum	 = initialPoint
		self.nPoints = 1.0
		self.mean	 = initialPoint

	def addPoint(self,point,tmpIndex):
		self.points.append( (point,tmpIndex) )
		self.sum += point
		self.nPoints += 1.0
		self.mean = self.sum / self.nPoints

	def __repr__(self):
		return self.points.__repr__()

	def __str__(self):
		return self.__repr__()
        
#=========================================================================
def simpleCluster( data, outlierThreshold ):
    
    # Copying and sorting the data
    tmpData = [(data[n],n) for n in xrange(len(data))]
    tmpData.sort()
    nPoints = len(data)

    # Selecting the first two points in the cluster
    clusters    = [ clusterClass( tmpData[0][0], tmpData[0][1] ) ]
    
    for n in xrange(1,nPoints):

        # Pull the next point
        tmpPoint, globalIndex = tmpData[n]
        
        # Which cluster are we closest to
        nClusters = len(clusters)
        d         = nClusters * [None]
        for m in xrange(nClusters):
            d[m] = ( tmpPoint - clusters[m].mean, m )
        #####
        
        # Which cluster are we closest to?
        d.sort()
        if ( d[0][0] < outlierThreshold ):
            # Add the point to the closest cluster
            clusters[d[0][1]].addPoint( tmpPoint, globalIndex )
        else:
            # Invoke a new cluster to handle the outlier
            clusters.append( clusterClass( tmpPoint, globalIndex ) )
        #####

    #####
    
	# Setting up the assignment
	clusterAssign = nPoints * [None]
	for n in xrange(len(clusters)):
	    for v_n, i_n in clusters[n].points:
	        clusterAssign[i_n] = n
	    #####
	#####
    
    return [clusters, clusterAssign]

#=========================================================================
class localScaffold( object ):

    def __init__(self, outlierThreshold ):
        self.markerPlacs = []
        self.outliers    = []

        self.length   = 0
        self.numBases = 0
        self.scaffID  = ''

        self.markerMidpoint = 0.0

        self.pearCorrel = 0.0
        self.revCorrel  = 0.0

        self.orientation = ''
        
        self.outlierThreshold = outlierThreshold
        
    def addPlac(self, parsedLine):
        self.length   = int( parsedLine[2] )
        self.numBases = int( parsedLine[1] )
        self.scaffID  = parsedLine[0]
        self.markerPlacs.append( markerPlac_class(parsedLine) )
        return

    def findOutliers_new(self):
        """
        If a marker is further away than outlierThreshold away from the rest of the
        markers, then it is considered an outlier.
        
        Not considered in:
          (a) midpoint calculation
          (b) orientation calculation
          (c) pearsons correlation
          
        Will be considered in:
          (a) marker positions output on the map.

        """

        tmpPlacs = [ tmpPlac.markerPos for tmpPlac in self.markerPlacs ]
                                             
        # Must have 2 markers to make a cluster of 2
        numPlacs = len(tmpPlacs)
        if ( numPlacs <= 2 ): return

        # Performing the clustering
        clusters, clusterAssign = simpleCluster( tmpPlacs, self.outlierThreshold )
        nClusters = len(clusters)
        
        # Forming the outlier list
        if ( nClusters > 1 ):
            # Finding the largest cluster
            DSU = [ (clusters[n].nPoints, n) for n in xrange(nClusters) ]
            DSU.sort(reverse=True)
            largeClusterIndex = DSU[0][1]
            for n in xrange(nClusters):
                if ( n == largeClusterIndex ): continue
                for v_n, i_n in clusters[n].points:
                    self.outliers.append( i_n )
                #####
            #####
        #####
        
#         if ( self.scaffID == 'Sh_240L12_contig-1' ):
#             print tmpPlacs
#             print clusters
#             print clusterAssign
#             print self.outliers
#             assert False
#         #####

        return

    def findOutliers_old(self):
        """
        If a marker is further away than outlierThreshold away from the rest of the
        markers, then it is considered an outlier.
        
        Not considered in:
          (a) midpoint calculation
          (b) orientation calculation
          (c) pearsons correlation
          
        Will be considered in:
          (a) marker positions output on the map.

        """

        from scipy.cluster.vq import kmeans2
        from numpy import array
        from numpy.linalg.linalg import LinAlgError

        tmpPlacs = array( [ [tmpPlac.markerPos] \
                                             for tmpPlac in self.markerPlacs ] )
                                             
        # Must have 2 markers to make a cluster of 2
        numPlacs = len(tmpPlacs)
        if ( numPlacs <= 2 ): return

        # Performing the clustering
        try:
            clustered = kmeans2( tmpPlacs, 2 ) # The 2 is for the number of clusters
        except LinAlgError:
            return
        #####
        bounds = [item[0] for item in clustered[0].tolist()]
        minVal = min(bounds)
        maxVal = max(bounds)
        
        if ( self.scaffID == 'Sh_215E18' ):
            print tmpPlacs
            print clustered
            assert False
        #####
        
        # Do we have any outliers?
        diff = (maxVal - minVal)
        if ( diff > self.outlierThreshold ):
            outlierList = clustered[1].tolist()
            if ( sum(outlierList) > (numPlacs/2) ):
                outlierID = 0
            else:
                outlierID = 1
            #####
            for n in xrange(numPlacs):
                if ( outlierList[n] == outlierID ):
                    self.outliers.append( n )
                #####
            #####
        #####

        return

    def computeMarkerMidpoint(self):
        tmpMid = 0.0
        for n in xrange( len(self.markerPlacs) ):
            # Excluding outliers
            if ( n not in self.outliers ):
                tmpMid += self.markerPlacs[n].markerPos
            #####
        #####
        numMarkers = len(self.markerPlacs) - len(self.outliers)
        self.markerMidpoint = tmpMid / float( numMarkers )
        return

    def getScreenedFwdPositions(self):
        numMarkers = len(self.markerPlacs) - len(self.outliers)
        X = numMarkers * [0]
        Y = numMarkers * [0]
        m = 0
        for n in xrange( len(self.markerPlacs) ):
            if ( n not in self.outliers ):
                X[m] = self.markerPlacs[n].markerPos
                Y[m] = self.markerPlacs[n].scaffStart
                m += 1
            #####
        #####
        return X, Y

    def getAllFwdPositions(self):
        numMarkers = len(self.markerPlacs)
        X = numMarkers * [0]
        Y = numMarkers * [0]
        for n in xrange( len(self.markerPlacs) ):
            X[n] = self.markerPlacs[n].markerPos
            Y[n] = self.markerPlacs[n].scaffStart
        #####
        return X, Y

    def getScreenedRevPositions(self):
        numMarkers = len(self.markerPlacs) - len(self.outliers)
        X = numMarkers * [0]
        Y = numMarkers * [0]
        m = 0
        for n in xrange( len(self.markerPlacs) ):
            if ( n not in self.outliers ):
                X[m] = self.markerPlacs[n].markerPos
                Y[m] = self.markerPlacs[n].revStart
                m += 1
            #####
        #####
        return X, Y

    def getAllRevPositions(self):
        numMarkers = len(self.markerPlacs)
        X = numMarkers * [0]
        Y = numMarkers * [0]
        for n in xrange( len(self.markerPlacs) ):
            X[n] = self.markerPlacs[n].markerPos
            Y[n] = self.markerPlacs[n].revStart
        #####
        return X, Y

    def computeOrientation(self, constraints):
        numMarkers = len(self.markerPlacs) - len(self.outliers)
        if ( numMarkers == 1 ):
            self.pearCorrel  = 1.0
            self.revCorrel   = 1.0
            self.orientation = 'fwd'
        else:
            # Forward Correlation
            X, Y = self.getScreenedFwdPositions()
            self.pearCorrel = pearsonsCorrelation(X,Y)

            # Reverse correlation
            X, Y = self.getScreenedRevPositions()
            self.revCorrel = pearsonsCorrelation(X,Y)

            # Selecting the best orientation
            if ( self.pearCorrel > self.revCorrel ):
                self.orientation = 'fwd'
            elif ( self.pearCorrel < self.revCorrel ):
                self.orientation = 'rev'
            else:
                self.orientation = 'fwd'
            #####
            # Evaluating the constraints
            if ( self.scaffID in constraints.keys() ):
                self.orientation = constraints[self.scaffID]
            #####
        #####
        return

#####

########################################################################
class buildChromosome( object ):

    def __init__(self, markerFile, \
                       markerSet, \
                       buildMode, \
                       assemblyFile, \
                       outlierThreshold, \
                       constraints={}, \
                       startScaffold=None, \
                       ignoreSet = set() ):

        self.markerFile       = markerFile
        self.markerSet        = markerSet
        self.constraints      = constraints
        self.buildMode        = buildMode
        self.startScaffold    = startScaffold
        self.ignoreSet        = ignoreSet
        self.scaffolds        = {}
        self.buffer           = 10000
        self.scaffoldLength   = {}
        self.outlierThreshold = outlierThreshold
        
        # Indexed FASTA for scaffold length
        self.indexedAssembly = assemblyFile

        # Building the output file
        tmpName = splitext( split(self.markerFile)[1] )[0]
        self.outputFile = '%s_locallyOrdered.dat'%tmpName

        # Building the chromosome
        self.getItDone()

    def getItDone(self):

        # (1) Read in contents of the file
        self.readMarkerFile()

        # (2) orderAndOrient
        if ( self.buildMode == 'manual' ):
            self.orderAndOrient_2()
        elif ( self.buildMode == 'auto' ):
            self.orderAndOrient()
        #####

        return

    def readMarkerFile(self):
    
        # Associating all markers with a scaffold
        for line in open(self.markerFile):
            # Parsing the line
            parsedLine = line.split(None)
            # Skipping spacer lines
            if ( parsedLine[0][0] == '-' ): continue
            # Screening for marker set
            if (self.markerSet != None):
                if (parsedLine[11] != self.markerSet):
                    continue
                #####
            #####
            scaffID     = parsedLine[0]
            scaffLength = int(parsedLine[1])
            # Ignoring appropriate scaffolds
            if ( scaffID in self.ignoreSet ): continue
            # Adding to the set
            try:
                self.scaffolds[scaffID].addPlac(parsedLine)
            except KeyError:
                self.scaffolds[scaffID] = localScaffold(self.outlierThreshold)
                self.scaffolds[scaffID].addPlac(parsedLine)
            #####
            # Reading in the scaffold length
            self.scaffoldLength[scaffID] = scaffLength
        #####
        
        # Computing the marker midpoints and orientation for all scaffolds
        for scaffID in self.scaffolds.keys():
#             print scaffID
            self.scaffolds[scaffID].findOutliers_new()
            self.scaffolds[scaffID].computeMarkerMidpoint()
            self.scaffolds[scaffID].computeOrientation( self.constraints )
        #####
        return

    def orderAndOrient(self):

        # (A) Construct an ordered list of scaffolds by midpoint, and by marker
        #     position and use this to determine the first scaffold.
        midPtOrderedList = []
        for scaffID, scaffClass in self.scaffolds.items():
            midPtOrderedList.append( (scaffClass.markerMidpoint, scaffID) )
        #####
        midPtOrderedList.sort()
        
#         for item in midPtOrderedList: print item
#         assert False
        
        # Just order it based on midpoint.  This seems to be the best way.
        # The statistical ordering tried to make everything linear, eventhough
        # it was not.
        scaffOrder = []
        for tmpMid, scaffID in midPtOrderedList:
            scaffOrder.append( (scaffID, self.scaffolds[scaffID].orientation) )
        #####

        # (D) Writing output file
        tmpHandle = open( self.outputFile, 'w' )
        for tuple in scaffOrder:
            tmpHandle.write( '%s\t%s\n'%(tuple[0],tuple[1]) )
        #####
        tmpHandle.write( '---------\n')
        X_in, Y_in = self.markerGraph( scaffOrder, True )
        
        ##################################
        r_1  = pearsonsCorrelation( X_in, Y_in )
        
        tmpHandle.write( 'Pearsons Correlation = %.5f\n'%r_1 )
        tmpHandle.write( '---------\n')
        X_all, Y_all = self.markerGraph( scaffOrder, False )
        for n in xrange( len(X_all) ):
            tmpHandle.write( '%.3f\t%d\n'%(X_all[n], Y_all[n]) )
        #####
        tmpHandle.close()


    def orderAndOrient_2(self):
        # Read in all of the information from the output file
        scaffOrder = []
        for line in open(self.outputFile, 'r'):
            parsedLine = line.split(None)
            try:
                if ( parsedLine[0][:8] == 'scaffold' ):
                    scaffID = parsedLine[0]
                    orient  = parsedLine[1]
                    scaffOrder.append( (scaffID,orient) )
                #####
            except IndexError:
                continue
            #####
        #####
        # Writing output file
        tmpHandle = open( self.outputFile, 'w' )
        for tuple in scaffOrder:
            tmpHandle.write( '%s\t%s\n'%(tuple[0],tuple[1]) )
        #####
        tmpHandle.write( '---------\n')
        X_in, Y_in = self.markerGraph( scaffOrder, True )
        r_1  = pearsonsCorrelation( X_in, Y_in )
        tmpHandle.write( 'Pearsons Correlation = %.5f\n\n'%r_1 )
        tmpHandle.write( '---------\n')
        X_all, Y_all = self.markerGraph( scaffOrder, False )
        for n in xrange( len(X_all) ):
            tmpHandle.write( '%.3f\t%d\n'%(X_all[n], Y_all[n]) )
        #####
        tmpHandle.close()

    def markerGraph( self, tmpScaffOrder, useScreenedPositions ):
        X_graph = []
        Y_graph = []
        offset  = 0
        for scaffID, scaffOrient in tmpScaffOrder:
#             print scaffID
            try:
                # Pulling positions
                if ( scaffOrient == 'fwd' ):
                    if ( useScreenedPositions ):
                        X, Y = self.scaffolds[scaffID].getScreenedFwdPositions()
                    else:
                        X, Y = self.scaffolds[scaffID].getAllFwdPositions()
                    #####
                else:
                    if ( useScreenedPositions ):
                        X, Y = self.scaffolds[scaffID].getScreenedRevPositions()
                    else:
                        X, Y = self.scaffolds[scaffID].getAllRevPositions()
                    #####
                #####
                # Building the marker graph
                for n in xrange( len(X) ):
                    X_graph.append( X[n] )
                    Y_graph.append( offset + Y[n] )
                #####
                # Incrementing the offset
                offset += self.buffer + self.scaffolds[scaffID].length
            except KeyError:
                # This is a scaffold that has no marker files
                # Need to find the length from the FASTA file
                offset += self.buffer + len(self.indexedAssembly[scaffID].seq)
            #####
        #####
        return X_graph, Y_graph

    def getScaffoldLength(self, scaffID):
        for record in iterBestHit( open(self.syntenyFile,'r') ):
            if ( record.scaffold == scaffID ): return record.scaffSize
        #####
        # Raise an error
#         print self.outputFile
        raise ValueError('Scaffold ID %s not found in %s'%(scaffID,\
                                                              self.syntenyFile))

#####


########################################################################
def real_main():

# Examples of what I have done with constraints, ignoreSet, and setting the start
#
#         markerFile = 'JOIN_2_Group_1_new.dat'
#         ignore_set = set()
#         ignore_set.add('super_360')
#         ignore_set.add('super_144')
#         ignore_set.add('super_37')
#         ignore_set.add('super_323')
#         ignore_set.add('super_351')
#         ignore_set.add('super_6')
#         ignore_set.add('super_653')
#         ignore_set.add('super_10')
#         ignore_set.add('super_28')
#         ignore_set.add('super_443')
#         ignore_set.add('super_11412')
#         ignore_set.add('super_11416')
#         ignore_set.add('super_744')
#         constraints = {}
#         constraints['super_1'] = 'rev'
#         buildChromosome( markerFile, constraints, ignoreSet=ignore_set )

# Handy chromosome building commands
#        grep 'super_11345\b' JOIN_2_Group_5b.dat
#
#        cat JOIN_2_Group_3_locallyOrdered.dat | grep 'super_11347\b' -B1 -A1
#
#        cat gg_synteny_JOIN_2.out | grep 'super_375\b' | sort -k8 -n

    if ( True ):
    
        # Setting the scaffold build mode
        # Two choices:  
        # auto:  automated ordering of the scaffolds with markers
        # manual:   manual manipulation of the scaffolds through the locallyOrdered files
        buildMode = 'auto'

        # Setup parameters        
        assemblyFile  = 'poplar14_5_snpFixed_phaseFixed.merged.broken.fasta'
        pre           = 'nisqually.synteny.broken'
        markerSet     = pre
        assemblyFile = FASTAFile_dict( assemblyFile )
        
        outlierThreshold = 10
        
        # Setting the number of chromosomes
        nameList = [ 'Chr%s' % ( str(n).zfill(2) ) for n in xrange(1,20) ]
        
        for chrID in nameList:

            markerFile  = '%s_Group_%s.dat'%(pre,chrID)
            excludeFile = '%s_excludeList.dat'%chrID
            ignore_set = set()
            if ( isfile(excludeFile) ):
                ignore_set = set( [line.strip() for line in open(excludeFile)] )
            #####
            constraints = {}
            buildChromosome( markerFile, markerSet, buildMode, assemblyFile, \
                             outlierThreshold, constraints, ignoreSet=ignore_set )
        #####
        
    #####

    # Finding out how many files a scaffold occurs in
    scaffoldDict     = {}
    for tmpName in nameList:
        outputFile = '%s_Group_%s_locallyOrdered.dat'%(pre,tmpName)
        for line in open(outputFile, 'r'):
            if ( line[0] == '-' ): break
            parsedLine = line.split(None)
            scaffID    = parsedLine[0]
            try:
                scaffoldDict[scaffID].append(tmpName)
            except KeyError:
                scaffoldDict[scaffID] = [tmpName]
            #####
        #####
    #####

    # Finding out how many times a scaffold occurs in a file
    scaffCounterDict = {}
    for tmpName in nameList:
        outputFile = '%s_Group_%s.dat'%(pre,tmpName)
        scaffCounterDict[tmpName] = {}
        for line in open(outputFile,'r'):
            parsedLine = line.split(None)
            if ( parsedLine[0][0] == '-' ): continue
            scaffID = parsedLine[0]
            # Adding to the set
            try:
                scaffCounterDict[tmpName][scaffID] += 1
            except KeyError:
                scaffCounterDict[tmpName][scaffID] = 1
            #####
        #####
    #####

    # Writing the full list
    print 'Full list of differences:'
    for scaffID, fileList in scaffoldDict.items():
        nFiles = len(fileList)
        if ( nFiles > 1 ):
            outputString = []
            outputString.append('%s'%scaffID)
            for n in xrange( len(fileList) ):
                outputString.append( '%s[%d]'%(fileList[n], \
                                   scaffCounterDict[fileList[n]][scaffID]) )
            #####
            print '\t'.join(outputString)
        #####
    #####

    # Putting it all together
    ignoreDict   = {}
    ambiguousSet = set()
    for scaffID, fileList in scaffoldDict.items():
        # If a scaffold occurs in multiple files then build the ignore list
        nFiles = len(fileList)
        if ( nFiles >= 2 ):
            # Counting the number of times a scaffold
            DSU = [ (scaffCounterDict[F][scaffID],F) for F in fileList ]
            DSU.sort(reverse=True)
            countSet = set( [ scaffCounterDict[F][scaffID] for F in fileList ] )
            if ( DSU[0][0] == DSU[1][0] ):
                for F in fileList:
                    try:
                        ignoreDict[F].append( scaffID )
                    except KeyError:
                        ignoreDict[F] = [ scaffID ]
                    #####
                #####
            else:
                for nCount,F in DSU[1:]:
                    try:
                        ignoreDict[F].append( scaffID )
                    except KeyError:
                        ignoreDict[F] = [ scaffID ]
                    #####
                #####
            #####
        #####
    #####

    print 'Ignore sets'
    DSU = [ (tmpFileName,scaffList) for tmpFileName, scaffList in ignoreDict.items()] 
    DSU.sort()
    for tmpFileName, scaffList in DSU:
        oh_exclude = open( '%s_excludeList.dat'%tmpFileName, 'w' )
        print tmpFileName
        for scaff in scaffList:
            print '        ignore_set.add(\'%s\')'%scaff
            oh_exclude.write( '%s\n'%scaff )
        #####
        print '==========='
        oh_exclude.close()
    #####

    print 'Ambiguous List'
    for scaffID, fileList in scaffoldDict.items():
        if ( scaffID in ambiguousSet ):
            nFiles = len(fileList)
            if ( nFiles > 1 ):
                outputString = []
                outputString.append('%s'%scaffID)
                for n in xrange( len(fileList) ):
                    outputString.append( '%s[%d]'%(fileList[n], \
                                       scaffCounterDict[fileList[n]][scaffID]) )
                #####
                print '\t'.join(outputString)
            #####
        #####
    #####

    # Making the join file
    findNum = re.compile( r'_([0-9]+)' ).findall
    tmpHandle = open( '%s_joins.dat'%pre, 'w' )
    for tmpName in nameList:
        outputFile = '%s_Group_%s_locallyOrdered.dat'%(pre,tmpName)
        outputString = []
        for line in open(outputFile, 'r'):
            if ( line[0] == '-' ): break
            parsedLine = line.split(None)
            scaffID    = parsedLine[0]
            index      = scaffID
            if ( parsedLine[1] == 'fwd' ):
                outputString.append( index )
            else:
                outputString.append( 'rc%s'%index )
            #####
        #####
        tmpHandle.write( '%s\n'%(' 10000 '.join(outputString) ) )
    #####
    tmpHandle.close()

    # Building the individual chromosomes
    if ( False and (buildMode == 'manual') ): 

       # ./gepardcmd.sh -seq1 /home/t3c2/AIP4/JWJ_ANALYSIS/Poplar/mapIntegration/joining_BROKEN_2/BROKEN_2_Group_1_chromosome_1.fasta -seq2 /home/t3c2/AIP4/JWJ_ANALYSIS/Poplar/mapIntegration/joining_BROKEN_2/BROKEN_2_Group_1_chromosome_1.fasta -outfile /home/t3c2/AIP4/JWJ_ANALYSIS/test.png -matrix matrices/edna.mat 
    
        # Building the chromosomes
        sepString = 10000 * 'N'
        stdout.write( 'BUILDING CHROMOSOMES\n')
        sizeCutoff = 5000000  #Breaking the chromosome into 5MB chuncks
        nCounter = 1
        for tmpName in nameList:
            # Opening the chromosome handle
            chrHandle = open( '%s_Group_%s_chromosome_%d.fasta'%(pre,tmpName,nCounter), 'w' )
            outputString = []
            for line in open('%s_Group_%s_locallyOrdered.dat'%(pre,tmpName), 'r'):
                parsedLine = line.split(None)
                try:
                    if ( parsedLine[0][:8] == 'scaffold' ):
                        scaffID = parsedLine[0]
                        tmpSeq = str(assemblyFile[scaffID].seq).upper()
                        if ( parsedLine[1] == 'fwd' ):
                            outputString.append( tmpSeq )
                            outputString.append( sepString )
                        else:
                            # Reverse complement the sequence
                            outputString.append( revComp(tmpSeq) )
                            outputString.append( sepString )
                        #####
                    #####
                except IndexError:
                    continue
                #####
                # Computing the size
                if ( sum([len(piece) for piece in outputString]) > sizeCutoff ):
                    # Writing the chromosome to file
                    stdout.write( '\tWriting %s_chromosome_%s_%d\n'%(pre,tmpName,nCounter) )
                    record = SeqRecord( seq = ''.join(outputString[:-1]), \
                                        id  = '%s_chromosome_%s_%d'%(pre,tmpName,nCounter), \
                                        description = '' )
                    writeFASTA( [record], chrHandle )
                    chrHandle.close()
                    # Increment the counter
                    nCounter += 1
                    chrHandle = open( '%s_Group_%s_chromosome_%d.fasta'%(pre,tmpName,nCounter), 'w' )
                    outputString = []
                #####
            #####
            # Writing the chromosome to file
            stdout.write( '\tWriting %s_chromosome_%s_%d\n'%(pre,tmpName,nCounter) )
            record = SeqRecord( seq = ''.join(outputString[:-1]), \
                                id  = '%s_chromosome_%s_%d'%(pre,tmpName,nCounter), \
                                description = '' )
            writeFASTA( [record], chrHandle )
            chrHandle.close()
            nCounter = 1
        #####
    #####


if ( __name__ == '__main__' ):
    real_main()
