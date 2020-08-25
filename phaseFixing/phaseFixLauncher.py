__author__ = "Jerry Jenkins, jjenkins@hudsonalpha.org"
__date__   = "Created:  7/10/18"

from phaseFix_lib import parseConfigFile, writeAndEcho, runPath, throwError
from phaseFix_lib import createGroupDictionary

from gp_lib import create_sh_file, submitJob, genePool_submitter

from hagsc_lib import iterFASTA, commify, writeFASTA

from os.path import join, isfile, splitext

from sys import argv, stderr

from math import sqrt

from os import chdir

#==============================================================
def real_main():
    
    # Read the basic information
    basePath   = argv[1]
    configFile = argv[2]
    phaseStep  = argv[3]
    
    # Changing to the base directory
    chdir( basePath )
    
    # Read the config file
    configDict, cmd_log = parseConfigFile( configFile, basePath )
    
    # Important output files
    primaryScaffsFile = join( basePath, 'primaryScaffolds.dat' )
    bamFile_FOFN      = join( basePath, '%s_bamFiles.fofn'   % configDict['RUNID'] )
    ARV_FOFN          = join( basePath, '%s_ARV_files.fofn'  % configDict['RUNID'] )
    CALL_FOFN         = join( basePath, '%s_CALL_files.fofn' % configDict['RUNID'] )
    
    # Run the commands
    if ( phaseStep == 'align_reads' ):
        
        # Setting the next step
        nextStep = 'read_to_variant'
        
        #-----------------------------
        # Order the scaffolds by size
        stderr.write( '- Computing the total number of bases in the genome\n' )
        DSU = [(len(r.seq),r.id) for r in iterFASTA(open(configDict['REFERENCE']))]
        DSU.sort(reverse=True)

        # Computing the target file size to balance the load
        targetNumBases = int( float(sum([item[0] for item in DSU])) / float(configDict['NUM_JOBS']) )
        
        # Group scaffolds starting with the largest
        oh           = open( primaryScaffsFile, 'w' )
        sizeCounter  = 0
        groupCounter = 1
        for scaffSize, scaffID in DSU:
            # Writing the scaffold to the current file
            oh.write( '%s\t%d\n'%(scaffID,groupCounter) )
            # Incrementing the size
            sizeCounter += scaffSize
            if ( sizeCounter > targetNumBases ):
                stderr.write( '\t\tGroup: %d, Size: %s, Target: %s\n'%( \
                              groupCounter, commify(sizeCounter), commify(targetNumBases)) )
                sizeCounter = 0
                groupCounter += 1
            #####
        #####
        stderr.write( '\t\tGroup: %d, Size: %s, Target: %s\n'%( groupCounter, commify(sizeCounter), commify(targetNumBases)) )
        oh.close()

        # Parameterization
        supressOutput   = False
    
        # Setting up the front command
        frontCommand  = 'ngmlr --bam-fix -t 4 -x pacbio -i %.2f --skip-write' % configDict['ALIGN_ID']
        
        # Generating the base output file name
        base_outFileName = join( basePath, configDict['RUNID'] )
        
        stderr.write( '\n\t-Running alignments on reads\n')
        
        # Reading in the chloroplast set
        if ( configDict['EXCLUDE'] == None ):
            excludeReadSet = set( )
        else:
            excludeReadSet = set( [ line[:-1] for line in open(configDict['EXCLUDE']) ] )
        #####
        
        # Setting up the bam file output list
        oh_bam = open( bamFile_FOFN, 'w' )
        
        #--------------------------------------
        # Running the jobs
        recordCounter = 0
        fileCounter   = 1
        tmpFileName   = '%s_%d.fasta' % ( base_outFileName, fileCounter )
        tmpHandle     = open( tmpFileName, 'w' )
        runningJobSet = set()
        runTime       = "5:00:00"
        runMemory     = configDict['ALIGN_MEM']
        targetRecords = configDict['NUM_READS']
        
        for record in iterFASTA( open( configDict['READS'] ) ):
            
            # Screening for chloroplast reads
            if ( record.id in excludeReadSet ): continue
            
            # Writing the scaffold to the current file
            writeFASTA( [record], tmpHandle )
            
            # Incrementing the size
            recordCounter += 1
            if ( recordCounter >= targetRecords ):
                
                # Closing the file
                stderr.write( '\t\tfileCounter: %d, Size: %d, Target: %d\n'%( fileCounter, recordCounter, targetRecords) )
                tmpHandle.close()

                # Creating the output bam file                
                localOutputBam = '%s.bam' % tmpFileName

                # Recording the bam file name
                oh_bam.write( '%s\n' % localOutputBam )
                
                # Create the command list and submit the job
                cmdList = []
                cmdList.append( 'source activate /global/dna/projectdirs/plant/geneAtlas/HAGSC_TOOLS/ANACONDA_ENVS/NGMLR_ENV/' )
                cmdList.append( 'rm -f %s'%( localOutputBam ) )
                cmdList.append( '%s -r %s -q %s | samtools view -F4 -Sbh - > %s'%( frontCommand, configDict['REFERENCE'], tmpFileName, localOutputBam ) )
                cmdList.append( 'rm -f %s' % tmpFileName )
                
                # Adding the job to a file
                shFileName = join( basePath, "%s_gp_NGMLR_%d.sh" % ( configDict['RUNID'], fileCounter ) )
#                 cmdList.append( 'rm -f %s' % shFileName )
                create_sh_file( cmdList, runTime, runMemory, shFileName, supressOutput )
                
                # Log the running job
                runningJobSet.add( submitJob( shFileName ) )
                stderr.write( "Starting %s, runningJobSet %d\n" % ( shFileName, len(runningJobSet) ) )
    
                # Opening the next file
                fileCounter += 1
                tmpFileName = '%s_%d.fasta'%(base_outFileName,fileCounter)
                tmpHandle   = open( tmpFileName, 'w' )
                recordCounter = 0
                
            #####
            
        #####
    
        #------------------------------
        # Writing the job set to a file
        writeAndEcho( '\n-----------------------\n- Writing jobIDs to the jobID file:', cmd_log )
        jobID_file = join( basePath, '%s_align_and_call_jobIDs.dat' % configDict['RUNID'] )
        oh_jobs    = open( jobID_file, 'w' )
        for jobID in runningJobSet:
            oh_jobs.write( '%s\n' % jobID )
        #####
        oh_jobs.close()

        # Handling the edge case of the last file
        if ( recordCounter < targetRecords ):
            
            # Closing the file
            stderr.write( '\t\tfileCounter: %d, Size: %d, Target: %d\n'%( fileCounter, recordCounter, targetRecords) )
            tmpHandle.close()

            # Creating the output bam file                
            localOutputBam = '%s.bam' % tmpFileName

            # Recording the bam file name
            oh_bam.write( '%s\n' % localOutputBam )
            
            # Create the command list and submit the job
            cmdList = []
            cmdList.append( 'source activate /global/dna/projectdirs/plant/geneAtlas/HAGSC_TOOLS/ANACONDA_ENVS/NGMLR_ENV/' )
            cmdList.append( 'rm -f %s'%( localOutputBam ) )
            cmdList.append( '%s -r %s -q %s | samtools view -F4 -Sbh - > %s'%( frontCommand, configDict['REFERENCE'], tmpFileName, localOutputBam ) )
            cmdList.append( 'rm -f %s' % tmpFileName )
            
            # Adding the job to a file
            shFileName = join( basePath, "%s_gp_NGMLR_%d.sh" % ( configDict['RUNID'], fileCounter ) )
#             cmdList.append( 'rm -f %s' % shFileName )
            
            # Adding the commands for the next stage
            cmdList.append( 'python %s/snpCalling/waitingForJobs.py %s'%( runPath, jobID_file ) )
            cmdList.append( 'python %s/phaseFixing/phaseFixLauncher.py %s %s %s'%(runPath, basePath, configFile, nextStep ) )
            runTime   = '24:00:00'
            runMemory = configDict['ALIGN_MEM']
            create_sh_file( cmdList, runTime, runMemory, shFileName, supressOutput )
            stderr.write( "Starting %s, runningJobSet %d\n" % ( shFileName, len(runningJobSet) ) )

            # Log the running job
            runningJobSet.add( submitJob( shFileName ) )
            stderr.write( "Starting %s, runningJobSet %d\n" % ( shFileName, len(runningJobSet) ) )
        
        else:
            
            #---------------------------
            # Launching the waiting job
            runTime   = '20:00:00'
            runMemory = '1G'
            cmdList = [ 'python %s/snpCalling/waitingForJobs.py %s'%( runPath, jobID_file ), \
                        'python %s/phaseFixing/phaseFixLauncher.py %s %s %s'%(runPath, basePath, configFile, nextStep ) ]
            sh_File = join( basePath, '%s.align_reads.waiting.sh'%configDict['RUNID'] )
            genePool_submitter( cmdList, runTime, runMemory, sh_File, supressOutput=False, delete_sh_file=False, waitUntilDone=False, emailAddress=configDict['EMAIL'], ncores=1 )

        #####
        
        # Closing the bamFile FOFN
        oh_bam.close()
        
    elif ( phaseStep == 'read_to_variant' ):
        
        # Setting the next step
        nextStep = 'partition_groups'

        # Assignment base command
        ARV_command = "python %s/phaseFixing/associate_read_with_variant.py" % runPath

        # Setting up the bam file output list
        oh_ARV   = open( ARV_FOFN, 'w' )

        # Looping over the bam files
        runningJobSet = set()
        runTime       = "2:00:00"
        runMemory     = "10G"
        supressOutput = False
        for line in open( bamFile_FOFN ):
            
            bamFile = line[:-1]
            
            #========================
            # Assign reads to variants command
            cmdList = ['source activate /global/dna/projectdirs/plant/geneAtlas/HAGSC_TOOLS/ANACONDA_ENVS/NGMLR_ENV/' ]
            cmdList.append( '%s %s %s %d %.2f %s' % ( ARV_command, 
                                                      bamFile, \
                                                      configDict['REFERENCE'], \
                                                      configDict['MIN_READ_SIZE'], \
                                                      configDict['ALIGN_COV'], \
                                                      configDict['VCF_FOFN'] ) )

            # Adding the job to a file
            shFileName = join( basePath, "%s.ARV.sh" % ( splitext( bamFile )[0] ) )
#             cmdList.append( 'rm -f %s' % shFileName )
            create_sh_file( cmdList, runTime, runMemory, shFileName, supressOutput )
            
            # Log the running job
            runningJobSet.add( submitJob( shFileName ) )
            stderr.write( "Starting %s, runningJobSet %d\n" % ( shFileName, len(runningJobSet) ) )
            
            # Writing the read_to_variant files to an FOFN
            READ_VAR = '%s.read_to_variant.dat' % splitext( bamFile )[0]
            oh_ARV.write( '%s\n' % READ_VAR )
            
        #####
        
        # Closing the ARV FOFN file
        oh_ARV.close()

        #------------------------------
        # Writing the job set to a file
        writeAndEcho( '\n-----------------------\n- Writing jobIDs to the jobID file:', cmd_log )
        jobID_file = join( basePath, '%s_read_to_variant_jobIDs.dat' % configDict['RUNID'] )
        oh_jobs    = open( jobID_file, 'w' )
        for jobID in runningJobSet:
            oh_jobs.write( '%s\n' % jobID )
        #####
        oh_jobs.close()

        #---------------------------
        # Launching the waiting job
        runTime   = '20:00:00'
        runMemory = '1G'
        cmdList = [ 'python %s/snpCalling/waitingForJobs.py %s'%( runPath, jobID_file ), \
                    'python %s/phaseFixing/phaseFixLauncher.py %s %s %s'%(runPath, basePath, configFile, nextStep ) ]
        sh_File = join( basePath, '%s.read_to_variant.waiting.sh'%configDict['RUNID'] )
        genePool_submitter( cmdList, runTime, runMemory, sh_File, supressOutput=False, delete_sh_file=False, waitUntilDone=False, emailAddress=configDict['EMAIL'], ncores=1 )
        
    elif ( phaseStep == 'partition_groups' ):
        
        # Setting the next stage
        nextStep = 'call_bases'
        
        # Read in the primary scaffolds file
        group_to_scaffList, scaff_to_group, finalGroup_ID = createGroupDictionary( primaryScaffsFile )
        
        # Open the output files for the primary scaffs file
        oh_out = {}
        
        # Loop over all of the read_to_variant files and partition the data out to each one
        currGrp = None
        for line in open( ARV_FOFN ):
            
            ARV_file = line[:-1]
            
            # Checking to see if the file exists
            if ( not isfile(ARV_file) ):
                throwError( '-------\nERROR:  %s ARV file does not exist\n-------\n' % ARV_file )
            #####
            
            for line in open( ARV_file ):

                # Computing the group ID
                if ( line[0] == ">" ):
                    readBase, scaffID, readStart, readEnd, readSeq = line[1:].split(None)
                    currGrp = scaff_to_group[scaffID]
                #####

                # Writing to a file
                try:
                    oh_out[currGrp].write( line )
                except KeyError:
                    oh_out[currGrp] = open( join( basePath, '%s.GRP_%s.read_to_variant.dat' % ( configDict['RUNID'], currGrp ) ), 'w' )
                    oh_out[currGrp].write( line )
                    print currGrp
                #####

            #####

        #####
        
        # Closing the files
        for currGrp in oh_out.iterkeys(): oh_out[currGrp].close()
        
        #---------------------------
        # Launching the waiting job
        runTime   = '20:00:00'
        runMemory = '1G'
        cmdList   = [ 'python %s/phaseFixing/phaseFixLauncher.py %s %s %s'%(runPath, basePath, configFile, nextStep ) ]
        sh_File   = join( basePath, '%s.call_bases.sh'%configDict['RUNID'] )
        genePool_submitter( cmdList, runTime, runMemory, sh_File, supressOutput=False, delete_sh_file=False, waitUntilDone=False, emailAddress=configDict['EMAIL'], ncores=1 )

    elif ( phaseStep == 'call_bases' ):
        
        # Setting the next step
        nextStep = 'make_fixing_file'

        # Read in the primary scaffolds file
        group_to_scaffList, scaff_to_group, finalGroup_ID = createGroupDictionary( primaryScaffsFile )

        # Assignment base command
        CALL_command = "python %s/phaseFixing/calling_ref_alt.py" % runPath

        # Looping over the bam files
        runningJobSet = set()
        runTime       = "5:00:00"
        runMemory     = "20G"
        supressOutput = False
        for GRP in group_to_scaffList.iterkeys():
        
            print GRP
            
            # Making the group variant file
            GRP_VAR_file = join( basePath, '%s.%s.read_to_variant.dat' % ( configDict['RUNID'], GRP ) )
            
            #========================
            # Assign reads to variants command
            cmdList = [ '%s %s %s %.3f' % ( CALL_command, 
                                            GRP_VAR_file, \
                                            configDict['REFERENCE'], \
                                            configDict['P_OUTLIER'] ) ]
            
            # Adding the job to a file
            shFileName = join( basePath, "%s.CALL_REF_ALT.%s.sh" % ( configDict['RUNID'], GRP ) )
#             cmdList.append( 'rm -f %s' % shFileName )
            create_sh_file( cmdList, runTime, runMemory, shFileName, supressOutput )
            
            # Log the running job
            runningJobSet.add( submitJob( shFileName ) )
            stderr.write( "Starting %s, runningJobSet %d\n" % ( shFileName, len(runningJobSet) ) )
            
        #####
        
        #------------------------------
        # Writing the job set to a file
        writeAndEcho( '\n-----------------------\n- Writing jobIDs to the jobID file:', cmd_log )
        jobID_file = join( basePath, '%s_call_bases_jobIDs.dat' % configDict['RUNID'] )
        oh_jobs    = open( jobID_file, 'w' )
        for jobID in runningJobSet:
            oh_jobs.write( '%s\n' % jobID )
        #####
        oh_jobs.close()

        #---------------------------
        # Launching the waiting job
        runTime   = '20:00:00'
        runMemory = '1G'
        cmdList = [ 'python %s/snpCalling/waitingForJobs.py %s'%( runPath, jobID_file ), \
                    'python %s/phaseFixing/phaseFixLauncher.py %s %s %s'%(runPath, basePath, configFile, nextStep ) ]
        sh_File = join( basePath, '%s.call_bases.waiting.sh'%configDict['RUNID'] )
        genePool_submitter( cmdList, runTime, runMemory, sh_File, supressOutput=False, delete_sh_file=False, waitUntilDone=False, emailAddress=configDict['EMAIL'], ncores=1 )

    elif ( phaseStep == 'make_fixing_file' ):
        
        # Read in the primary scaffolds file
        group_to_scaffList, scaff_to_group, finalGroup_ID = createGroupDictionary( primaryScaffsFile )
        
        # Assignment base command
        CALL_command = "python %s/phaseFixing/makeFixingFile.py" % runPath
        
        # Loop over the outlier files and compute coverage statistics
        x  = 0.0
        x2 = 0.0
        N  = 0
        for GRP in group_to_scaffList.iterkeys():
            outlier_file = join( basePath, '%s.%s.outlier.dat' % ( configDict['RUNID'], GRP ) )
            # Computing the coverage on a per snp basis
            for line in open(outlier_file):
                scaffID, varID, d1, d2 = line.split(None)
                cov = int(d1) + int(d2)
                x += cov
                x2 += cov * cov
                N += 1
            #####
        #####
        
        # Computing the coverage statistics
        x_bar  = float(x) / float(N)
        x2_bar = float(x2) / float(N)
        stdDev = sqrt( x2_bar - x_bar * x_bar )
        
        # Setting the bounds on the coverage to accept
        alpha  = 1.3
        minCov = max( 4.0, x_bar - alpha * stdDev )
        if ( configDict['MAX_READ_COV'] == None ):
            maxCov = x_bar + alpha * stdDev
        else:
            maxCov = configDict['MAX_READ_COV']
        #####
        
        # Checking to see if the max coverage is larger than min coverage
        if ( minCov >= maxCov ):
            throwError( 'Minimum coverage %.3f is >= maximum coverage of %.3f\nPlease Reset MAX_READ_COV in the phaseFix.config file.' % (minCov, maxCov) )
        #####
        
        writeAndEcho( 'avgCov = %.3f' % x_bar,  cmd_log )
        writeAndEcho( 'stdCov = %.3f' % stdDev, cmd_log )
        writeAndEcho( 'minCov = %.3f' % minCov, cmd_log )
        writeAndEcho( 'maxCov = %.3f' % maxCov, cmd_log )
        
        #---------------------------
        # Launching the fixing file
        runTime   = '2:00:00'
        runMemory = '2G'
        
        #========================
        # Assign reads to variants command
        cmdList = [ '%s %s %.3f %.3f %s %s %s' % ( CALL_command, 
                                                   configDict['REFERENCE'], \
                                                   minCov, maxCov, \
                                                   primaryScaffsFile, \
                                                   basePath, \
                                                   configDict['RUNID'] ) ]
        
        sh_File = join( basePath, '%s.make_fixing_file.sh'%configDict['RUNID'] )
        genePool_submitter( cmdList, runTime, runMemory, sh_File, supressOutput=False, delete_sh_file=False, waitUntilDone=False, emailAddress=configDict['EMAIL'], ncores=1 )
        
    #####
    
#==============================================================
if ( __name__ == '__main__' ):
    real_main()
