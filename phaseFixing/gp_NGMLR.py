#!/usr/common/usg/languages/python/2.7-anaconda/bin/python
__author__="Jerry Jenkins, jjenkins@hudsonalpha.org"
__date__ ="Created:  7/24/15"

from hagsc_lib import generateTmpDirName, baseFileName, deleteDir
from hagsc_lib import testDirectory, iterCounter, iterFASTA, writeFASTA

from gp_lib import create_sh_file, submitJob

from os.path import abspath, curdir, join

from sys import stderr

from os import mkdir

import subprocess

import time

#==============================================================
def real_main():

    # Parameterization
    targetRecords   = 30000
    jobsInQueue     = 200
    queueID         = 'all.q'
    targetFile      = '/projectb/scratch/jjenkins/Poplar_14_5_phaseFixing/combined.snpFixed.fasta'
    queryFile       = '/projectb/scratch/jjenkins/Poplar_14_5_phaseFixing/poplar_14.5.fasta'
    supressOutput   = False

    # Setting up the front command
    frontCommand  = 'ngmlr --bam-fix -t 4 -x pacbio -i 0.80 --skip-write'
    
    # Generating the directory information
    basePath = abspath(curdir)
    tmpDir   = 'phaseFixingDirectory'
    tmpPath  = join( basePath, tmpDir )

    # Writing the job IDs to a file
    job_ID_File = join( tmpPath, 'jobIDs.dat')

    # Generating the base output file name
    base_outFileName = join( tmpPath, 'combined.snpFixed' )

    stderr.write( '\n\t-Creating tmp directory\n')
    
    # Reading in the chloroplast set
    chloroReadSet = set([line[:-1] for line in open('/projectb/scratch/jjenkins/Poplar_14_5_phaseFixing/chloroplastNames.txt')])
    
    #--------------------------------------
    # Running the jobs
    recordCounter = 0
    fileCounter   = 1
    tmpFileName   = '%s_%d.fasta' % ( base_outFileName, fileCounter )
    tmpHandle     = open( tmpFileName, 'w' )
    x             = iterCounter( 1000000 )
    runningJobSet = set()
    jobNum        = 0
    timeInc       = 10
    runTime       = "5:00:00"
    runMemory     = "5G"
#     p_fasta       = subprocess.Popen( 'cat %s | bunzip2 -c' % queryFile, shell=True, stdout=subprocess.PIPE )
    for record in iterFASTA( open( queryFile ) ):
        
        # Screening for chloroplast reads
        if ( record.id in chloroReadSet ): continue
        
        # Writing the scaffold to the current file
        writeFASTA( [record], tmpHandle )
        x()
        
        # Incrementing the size
        recordCounter += 1
        if ( recordCounter >= targetRecords ):
            
            # Closing the file
            stderr.write( '\t\tfileCounter: %d, Size: %d, Target: %d\n'%( \
                          fileCounter, recordCounter, targetRecords) )
            tmpHandle.close()
            
            # Create the command list and submit the job
            localOutputFile = '%s.bam' % tmpFileName
            cmdList = []
            cmdList.append( 'source activate /global/dna/projectdirs/plant/geneAtlas/HAGSC_TOOLS/ANACONDA_ENVS/NGMLR_ENV/' )
            cmdList.append( 'rm -f %s'%( localOutputFile ) )
            cmdList.append( '%s -r %s -q %s | samtools view -F4 -Sbh - > %s'%( frontCommand, targetFile, tmpFileName, localOutputFile ) )
            cmdList.append( 'rm -f %s' % tmpFileName )

            #========================
            cmdList.append( 'python associate_reads_with_variants.py' ) # Need to add this as part of the pipeline
            #========================
            
            # Adding the job to a file
            shFileName = join( tmpPath, "gp_NGMLR_%d.sh" % fileCounter )
            cmdList.append( 'rm -f %s' % shFileName )
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

        # Keeping the number of jobs in the queue
        while ( len(runningJobSet) > jobsInQueue ):
        
            stderr.write( "%d running jobs\n" % len(runningJobSet) ) 
            
            # Wait a specific amount of time
            time.sleep(timeInc)
            
            # Extracting all job IDs
            p_job = subprocess.Popen( 'squeue', shell=True, stdout=subprocess.PIPE )
            checkerSet = set( [line.split(None)[0] for line in p_job.stdout if (line[0] not in 'j-')] )
            p_job.poll()
            
            # Computing the shared IDs
            sharedIDs     = runningJobSet.intersection(checkerSet)
            IDs_to_delete = runningJobSet - sharedIDs
            
            # If the number of shared IDs is too low
            if ( len(IDs_to_delete) > 0 ):
                for item in IDs_to_delete: stderr.write( "\t-JobID %s Completed\n"%item )
                # Pulling the jobIDs from the running IDs
                map( runningJobSet.remove, IDs_to_delete )
                # Exit the loop and launch more jobs
                break
            #####
            
        #####
        
    #####

#     p_fasta.poll()

    tmpHandle.close()
    stderr.write( '\t\tfileCounter: %d, Size: %d, Target: %d\n'%( \
                                fileCounter, recordCounter, targetRecords) )

#==============================================================
if ( __name__ == '__main__' ):
    real_main()
