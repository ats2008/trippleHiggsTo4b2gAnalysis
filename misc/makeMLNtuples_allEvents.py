#!/usr/bin/env python 
from __future__ import print_function
import os
import copy,json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-s',"--submit", help="Submit file to condor pool", action='store_true' )
parser.add_argument('-r',"--resubmit", help="Re-Submit file to condor pool", action='store_true' )
parser.add_argument('-t',"--test", help="Test Job", action='store_true' )
parser.add_argument('-p',"--printOnly", help="Only Print the commands", action='store_true' )
parser.add_argument('-n',"--njobs", help="Number of jobs to make",default='-6000')
parser.add_argument('-e',"--nevts", help="Number of events per job",default='-1')
parser.add_argument('-v',"--version", help="Vesion of the specific work",default='1p0')

args = parser.parse_args()

version=args.version
isTest=args.test
if isTest:
    args.njobs=2
    args.nevts=1000

njobs=int(args.njobs)
maxevents=int(args.nevts)
max_meterialize=250
jobsToProcess=None
submit2Condor=args.submit
resubmit2Condor=args.resubmit
onlyPrint=args.printOnly
print(" submit jobs ",submit2Condor)
print(" resubmit jobs ",resubmit2Condor)
print(" isTest ",isTest)
print(" printOnly ",onlyPrint)
print(" njobs ",njobs)
print(" maxEvt ",maxevents)

if submit2Condor or resubmit2Condor:
    choice=raw_input("Do you really want to submit the jobs to condor pool ? ")
    if 'y' not in choice.lower():
        print("Exiting ! ")
        exit(0)

jobDict={}
with open('misc/jsons/mlNtuplizer_allEvents.json') as f:
    jobDict=json.load(f)

fileListDict={}
with open('misc/jsons/fileListToUse.json') as f:
    fileListDict=json.load(f)

templateCMD="""
./misc/condorJobMakerGeneric.py 
       --exe             @@EXE 
       --fsrc            @@FILELIST 
       --runScript       @@SCRIPT_TPL 
       --cfg             @@CFG_TPL 
       --dest            @@DESTN 
       --jn              @@NJOBS 
       --fn              @@FILES_PER_JOB 
       --maxEvt          @@MAXEVENTS 
       --tag             @@TAG 
       --jobType         @@JOB_TYPE
       --maxMeterialize  @@MAX_METERIALIZE
"""

if onlyPrint :
    templateCMD= 'echo '+templateCMD 

allCondorSubFiles=[]
if jobsToProcess==None:
    jobsToProcess=list( jobDict.keys() )
jobsToProcess=['bkg','sig']
for jobTag in jobsToProcess:
    if njobs > 0:
        jobDict[jobTag]['njobs']=njobs
    if maxevents > 0:
        jobDict[jobTag]['maxevents']=maxevents
    cmd=templateCMD.replace("@@EXE",jobDict[jobTag]['exe'])
    cmd=cmd.replace("@@CFG_TPL"        ,jobDict[jobTag]['cfg'])
    cmd=cmd.replace("@@SCRIPT_TPL"     ,jobDict[jobTag]['script'])
    cmd=cmd.replace("@@NJOBS"          ,str(jobDict[jobTag]['njobs']))
    cmd=cmd.replace("@@JOB_TYPE"          ,str(jobDict[jobTag]['jobType']))
    cmd=cmd.replace("@@FILES_PER_JOB"        ,str(jobDict[jobTag]['files_per_job']))
    cmd=cmd.replace("@@MAXEVENTS"        , str(jobDict[jobTag]['maxevents']) )
    cmd=cmd.replace("@@MAX_METERIALIZE"        , str(jobDict[jobTag]['max_meterialize']) )
    for key in jobDict[jobTag]:
        if '@@' in key:
            cmd=cmd.replace( key ,jobDict[jobTag][ key ]  )

    for dset in jobDict[jobTag]['datasets']:
        dset=dset.encode('ascii','replace')
        if dset in fileListDict:
            flist=fileListDict[dset]
        else :
            print("File list for dataset : ",dset," not found ! skipping the datset ")
            continue
        tag=jobTag+'_'+dset+'_'+version
        cmdD=cmd.replace("@@TAG",tag)
        cmdD=cmdD.replace("@@FILELIST",flist)
        destination = jobDict[jobTag]['destn'] +'/'+tag+'/'
        cmdD=cmdD.replace("@@DESTN",destination)
        
        head='Condor/'+jobDict[jobTag]['jobType']+'/Jobs'+tag
        condorScriptName=head+'/job'+tag+'.sub'
        allCondorSubFiles.append(condorScriptName)
        if resubmit2Condor:
            continue
        print("")
        print("")
        print("="*40)
        print(cmdD)
        print(" - -"*10)
        print("")
        print("")
        cmdD=cmdD.replace("\n","")
        os.system( cmdD )
        if isTest:
            break
print("")
print("")
print("All condor submit files to be submitted ")
for fle in allCondorSubFiles:
    print('condor_submit '+fle)
    if submit2Condor or resubmit2Condor:
        os.system('condor_submit '+fle)
print("")
print("")
