#!/usr/bin/env python 
from __future__ import print_function
import os,glob
import copy,json
import argparse

condorScriptString="\
executable = $(filename)\n\
output = $Fp(filename)run.$(Cluster).stdout\n\
error = $Fp(filename)run.$(Cluster).stderr\n\
log = $Fp(filename)run.$(Cluster).log\n\
request_cpus = 3\n\
"
pwd=os.environ['PWD']
proxy_path=os.environ['X509_USER_PROXY']
HOME=os.environ['HOME']
NJOBS=20000
NEVENTS_PER_JOB = -1
ZERO_OFFSET=0
FILES_PER_JOB=1
maxMeterialize=100
offsetStep=10000
RESULT_BASE='/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/results'
JOB_TYPE='systExport'

parser = argparse.ArgumentParser()
parser.add_argument('-s',"--submit", help="Submit file to condor pool", action='store_true' )
parser.add_argument('-r',"--resubmit", help="Re-Submit file to condor pool", action='store_true' )
parser.add_argument("--resetSucess", help="Reset the sucessfull jobs as Sucessfull", action='store_true' )
parser.add_argument('-t',"--test", help="Test Job", action='store_true' )
parser.add_argument('-p',"--printOnly", help="Only Print the commands", action='store_true' )
parser.add_argument(     "--doNominal", help="Submit only Nominal Ones ! ", action='store_true' )
parser.add_argument('-n',"--njobs", help="Number of jobs to make",default='-6000')
parser.add_argument('-e',"--nevts", help="Number of events per job",default='-1')
parser.add_argument('-c',"--cats", help="Categories of jobs to process ",default=None)
parser.add_argument('-f',"--flist", help="Files to process",default=None)
parser.add_argument('-y',"--years", help="Year of jobs to process ",default="all")
parser.add_argument('-v',"--version", help="Vesion of the specific work",default='TEST')

args = parser.parse_args()

version=args.version

njobs=int(args.njobs)
maxevents=int(args.nevts)
max_meterialize=250
jobsToProcess=None
submit2Condor=args.submit
resubmit2Condor=args.resubmit
isTest=args.test
onlyPrint=args.printOnly
cats=args.cats
years=args.years.split(",")
job_hash=f'systExport_{args.version}'

procToProcess=[]
print(" submit jobs ",submit2Condor)
print(" resubmit jobs ",resubmit2Condor)
print(" isTest ",isTest)
print(" printOnly ",onlyPrint)
print(" njobs ",njobs)
print(" maxEvt ",maxevents)
print(" cats ",cats)
print(" years ",years)

if submit2Condor or resubmit2Condor:
    choice=input("Do you really want to submit the jobs to condor pool ? ")
    if 'y' not in choice.lower():
        print("Exiting ! ")
        exit(0)

jobDict={}
jobConfigJson='misc/jsons/systExtractor.json'
with open(jobConfigJson) as f:
    jobDict=json.load(f)


fileMapName='fileList/nanoFiles.jsob'
with open(fileMapName) as f:
    fileMap=json.load(f)

if isTest :
    njobs=2
    maxevents=1000
allCondorSubFiles=[]
jobsToProcess=list( jobDict.keys() )
#jobsToProcess=['sig']
#jobsToProcess=['data']
targetYear='run2'
#if cats:
#    jobsToProcess=cats.split(",")

print()
print(" Processing Job categories : ",jobsToProcess)
print()
i=0
njobs_total=0
for jobTag in jobsToProcess:
    htag=jobTag
    
    runScriptTemplate = jobDict[htag]['script'] 
    runScriptTxt=[]
    with open(runScriptTemplate,'r') as f:
        runScriptTxt=f.readlines()
    runScriptTxt=''.join(runScriptTxt)
    
    executable=jobDict[htag]["exe"]
    njobs=0
    print(htag,jobsToProcess)
    for dset in fileMap:
        
        i=0
        #if dset!="/ttHJetToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM":
        #if dset!="/VBFHToGG_M125_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM":   
        #if "TTHHToGGBB_TuneCP5" not in dset:    
        #    continue

        head=pwd+f'/Condor/{JOB_TYPE}/{job_hash}/{htag}/{dset}/'
        print(f"\rMaking Jobs for {htag}/{dset}\t\t\t\t",end="")
        for fls in fileMap[dset]:
            dirName =f'{head}/Job_{i}/'
            runScriptName=dirName+f'/{htag}_{i}_run.sh'
            destination=f'{RESULT_BASE}/{JOB_TYPE}/{job_hash}/{htag}/{dset}/'
            foutNmae=destination+f'output_f{i}.root'
            tmp=runScriptTxt.replace("@@DIRNAME",dirName)
            tmp=tmp.replace("@@IFNAME",fls)
            tmp=tmp.replace("@@IDX",str(i))
            tmp=tmp.replace("@@TAG",f"f{i}")
            tmp=tmp.replace("@@RUNSCRIPT",runScriptName)
            tmp=tmp.replace("@@EXECUTABLE",executable)
            tmp=tmp.replace("@@proxy_path",proxy_path)
            tmp=tmp.replace("@@HOME",HOME)
            tmp=tmp.replace("@@DESTINATION",destination)
            i+=1
            if args.resetSucess:
                if os.path.exists(foutNmae):
                    if not os.path.exists(f'{runScriptName+".sucess"}'):
                        cmd=f'mv {runScriptName} {runScriptName+".sucess"}'
                        os.system(cmd)
                else:
                    njobs_total+=1
                continue
            njobs+=1
            njobs_total+=1
            if args.resubmit:
                continue
            if not args.printOnly:
                if not os.path.exists(dirName):
                    os.system('mkdir -p '+dirName)
                if not os.path.exists(destination):
                    os.system('mkdir -p '+destination)
                runScript=open(runScriptName,'w')
                if os.path.exists(runScriptName+'.sucess'):
                   os.system('rm '+runScriptName+'.sucess')
                runScript.write(tmp)
                runScript.close()
                os.system('chmod +x '+runScriptName)
            if isTest:
                break
        condorScriptName=head+'/job.sub'
        allCondorSubFiles.append(condorScriptName)
        if args.resetSucess:
            continue
        if resubmit2Condor:
            continue
        if not onlyPrint:
            with open(condorScriptName,'w') as condorScript:
                condorScript.write(condorScriptString)
                condorScript.write("queue filename matching ("+head+"/*/*.sh)\n")
                print(f"Condor {njobs} Jobs made !\n\t submit file  : {condorScriptName}")
        print()
        print("Made ",i," jobs for ",dset)
        if isTest:
            break
print("")
print("")
if not args.printOnly:
    print("All condor submit files to be submitted ")
    for fle in allCondorSubFiles:
        if submit2Condor or resubmit2Condor:
            os.system('condor_submit '+fle)
        else:
            if args.resetSucess:
                pass
            else:
                print('condor_submit '+fle)
print("Total number of jobs ",njobs_total)
print("")
