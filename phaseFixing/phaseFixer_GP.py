#!/usr/common/software/python/2.7-anaconda-2019.07/bin/python
__author__="Jerry Jenkins, jjenkins@hudsonalpha.org"
__date__ ="Created:  7/10/158"

from phaseFix_lib import parseConfigFile, writeAndEcho, throwError, jobTypeSet, runPath

from gp_lib import genePool_submitter

from optparse import OptionParser

from os.path import abspath, join

from os import curdir, system

#==============================================================
def real_main():
    
    # Defining the program options
    usage  = "usage: %prog [options]"

    parser = OptionParser(usage)

    restartJob_def = 'align_reads'
    parser.add_option( "-s", \
                       "--restartJob", \
                       type    = 'str', \
                       help    = "Start at this point in the pipeline.  Default: %s"%restartJob_def, \
                       default = restartJob_def )

    # Parsing the arguments
    (options, args) = parser.parse_args()
    
    # Finding the base path
    basePath = abspath( curdir )
    
    # Setting the config file name
    configFile = join( basePath, 'phaseFix.config' )
    
    # Read the config file
    configDict, cmd_log = parseConfigFile( configFile, basePath )
    
    # Creating the initial BAM file to start the compuation off
    if ( options.restartJob not in jobTypeSet ):
        throwError( '%s is not a valid step in the SNP/INDEL calling pipeline'%options.restartJob )
    #####
    nextStep = options.restartJob
    
    # Building the command set
    cmdList = [ 'python %s/phaseFixing/phaseFixLauncher.py %s %s %s'%( runPath, basePath, configFile, nextStep) ]
    
    # Writing to the log file
    writeAndEcho( '\n-----------------------\n- Starting the computation:', cmd_log )
    for cmd in cmdList: writeAndEcho( cmd, cmd_log )
    
    # Running the job
    maxTime   = "05:00:00"
    maxMemory = '10G'
    sh_File   = join( basePath, '%s_phaseFixerJob.sh'%configDict['RUNID'] )
    genePool_submitter( cmdList, maxTime, maxMemory, sh_File, supressOutput=False, delete_sh_file=False, waitUntilDone=False, ncores=1 )

#==============================================================
if ( __name__ == '__main__' ):
    real_main()
    