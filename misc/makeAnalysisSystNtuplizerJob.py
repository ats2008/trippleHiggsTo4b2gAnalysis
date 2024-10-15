#!/usr/bin/env python3
from __future__ import print_function
import os
import copy,json,datetime
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-s',"--submit", help="Submit file to condor pool", action='store_true' )
parser.add_argument('-r',"--resubmit", help="Re-Submit file to condor pool", action='store_true' )
parser.add_argument('-t',"--test", help="Test Job", action='store_true' )
parser.add_argument('-p',"--printOnly", help="Only Print the commands", action='store_true' )
parser.add_argument("--dryRun", help="do Dry run", action='store_true' )
parser.add_argument("--doNominal", help="Only do nominal", action='store_true' )
parser.add_argument("--doOnlySyst", help="Only do systematics", action='store_true' )
parser.add_argument('-n',"--njobs", help="Number of jobs to make",default='-6000')
parser.add_argument('-e',"--nevts", help="Number of events per job",default='-1')
parser.add_argument('-j',"--nfiles_per_job", help="Number of files per job",default=-1,type=int)
parser.add_argument('-c',"--cats", help="Categories of jobs to process ",default=None)
parser.add_argument(     "--procs", help="procs jobs to process ",default=None)
parser.add_argument('-y',"--years", help="Year of jobs to process ",default="all")
parser.add_argument('-v',"--version", help="Vesion of the specific work",default='TEST')
parser.add_argument("--club_hash", help="club hash if to resubmit",default=None)

args = parser.parse_args()

version=args.version

njobs=int(args.njobs)
doOnlyNominal=int(args.doNominal)
maxevents=int(args.nevts)
max_meterialize=250
jobsToProcess=None
submit2Condor=args.submit
resubmit2Condor=args.resubmit
isTest=args.test
onlyPrint=args.printOnly
cats=args.cats
years=args.years.split(",")
club_hash= args.club_hash
if not club_hash:
    club_hash=datetime.datetime.now().strftime("%d%b%A%I%p%Mm%Ss")
if args.version!='TEST':
    club_hash=f'{args.version}_{club_hash}'
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
with open('misc/jsons/analysisWithSystNtuplizer.json') as f:
    jobDict=json.load(f)

fileListDict={}
#with open('results/mcDPMCleaner/mcSig.json') as f:
#with open('results/systSpliNtuples/systFileMap_v6p4.json') as f:
with open('results/systSpliNtuples/systFileMap_v6p4_all.json') as f:
    fileListDict_=json.load(f)
    #fileListDict={ dset.encode('ascii','replace') : fileListDict_[dset]  for dset in fileListDict_ }
    fileListDict=fileListDict_
#for ky in fileListDict:
#    print(ky,fileListDict[ky])
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
       --cfgExtras \"@@CFG_EXTRAS\"
"""

if isTest :
    njobs=2
    maxevents=1000
if onlyPrint:    
    templateCMD= 'echo '+templateCMD 
allCondorSubFiles=[]
if jobsToProcess==None:
    jobsToProcess=list( jobDict.keys() )

peakingProcs=['singleH','doubleH','sig','kappaScan']
jobsToProcess=['bkg','singleH','doubleH','sig','data','qcd','kappaScan','ttX']
procsToProcess=None
if cats:
    jobsToProcess=cats.split(",")
if args.procs:
    procsToProcess=args.procs.split(",")

print()
print(" Processing Job categories : ",jobsToProcess)
print(" Processing Job Prcess  : ",procsToProcess)
print()
nJobsTobeMade=0
for jobTag in jobsToProcess:
    #print(jobDict[jobTag])
    if njobs > 0:
        jobDict[jobTag]['njobs']=njobs
    if maxevents > 0:
        jobDict[jobTag]['maxevents']=maxevents
    if args.nfiles_per_job > 0:
        jobDict[jobTag]['files_per_job']=args.nfiles_per_job
    
    #print(jobDict[jobTag])
    cmd=templateCMD.replace("@@EXE",jobDict[jobTag]['exe'])
    cmd=cmd.replace(  "@@CFG_TPL"        ,jobDict[jobTag]['cfg'])
    cmd=cmd.replace(  "@@SCRIPT_TPL"     ,jobDict[jobTag]['script'])
    cmd=cmd.replace(  "@@NJOBS"          ,str(jobDict[jobTag]['njobs']))
    cmd=cmd.replace(  "@@JOB_TYPE"          ,str(jobDict[jobTag]['jobType'])+f'/{club_hash}')
    cmd=cmd.replace(  "@@FILES_PER_JOB"        ,str(jobDict[jobTag]['files_per_job']))
    cmd=cmd.replace(  "@@MAXEVENTS"        , str(jobDict[jobTag]['maxevents']) )
    cmd=cmd.replace(  "@@MAX_METERIALIZE"        , str(jobDict[jobTag]['max_meterialize']) )
    #cmd=cmd.replace(  "@@MAX_METERIALIZE"        , str(1000) )
    cfgExtra='""'
    for key in jobDict[jobTag]:
        if '@@' in key:
            cmd=cmd.replace( key ,jobDict[jobTag][ key ]  )

    for dset in jobDict[jobTag]['datasets']:
        if procsToProcess is not None:
            if len(procsToProcess):
                skp=True
                for prc in procsToProcess:
                    if prc in dset:
                        skp=False
                        break
                if skp:
                    #print("Skipping ",jobTag,dset)
                    continue
                   
        if "all" not in years:
            hasYear=False
            for yr in years:
                if yr in dset:
                    hasYear=True
                    break
            if not hasYear:
                #print("Skipping ",dset)
                continue
        #dset=dset.encode('ascii','replace')
        if doOnlyNominal:
            if jobTag in peakingProcs:
                #if 'nominal' not in dset:
                if 'sigma' in dset:
                    continue



        if args.doOnlySyst:
            if jobTag in peakingProcs:
                if 'nominal' in dset:
                    continue
            else:
                continue
        print("here",dset,years,procsToProcess)

        if dset in fileListDict:
            flist=fileListDict[dset]
        else :
            print("File list for dataset : ",dset," not found ! skipping the datset ")
            continue
        if args.dryRun:
            with open(flist) as ff:
                txt=ff.readlines()
                for k in txt:
                    if 'sig_kappaScan_mc_UL2016Post_c3_1_d4_0_v30/sig_kappaScan_mc_UL2016Post_c3_1_d4_0_v30_UL2016Post_0_0_FNUFEEDown01sigma_scaled_split.root' in k:
                        print(k)
                        print(flist)
                nf=len(txt)
                nn=int(nf / int(jobDict[jobTag]['files_per_job']) + 0.99999 )
            if nf < 2:
                print(flist)
            nJobsTobeMade+= nn
            print("Processing ",dset," njobs : ",nn," [ ", nf ,"/", jobDict[jobTag]['files_per_job']  ,"] "," Total jobs now  : ",nJobsTobeMade)
            print("    > filelist : ",flist)
            continue
        else:
            print("Processing ",dset)
        #if ('JEC' in dset) and ('2016' in dset):
        #    continue
        tag=jobTag+'_'+dset+'_'+version
        cmdD=cmd.replace("@@TAG",tag)
        print("\n\n\n\n")
        print(f"{flist}")
        cmdD=cmdD.replace("@@FILELIST",flist)
        destination = jobDict[jobTag]['destn'] +f'/{club_hash}/{version}/{jobTag}'+'/'+dset+'/'
        cmdD=cmdD.replace("@@DESTN",destination)
        cfgExtra_=cfgExtra
        if 'cfgExtras' in jobDict[jobTag]:
            if "*" in jobDict[jobTag]['cfgExtras'] :
                cfgExtra_= jobDict[jobTag]['cfgExtras']["*"]   
            if dset in jobDict[jobTag]['cfgExtras']:
                cfgExtra_=  jobDict[jobTag]['cfgExtras'][dset]  
        cmdD=cmdD.replace("@@CFG_EXTRAS",cfgExtra_)
        
        head='Condor/'+jobDict[jobTag]['jobType']+f'/{club_hash}'+'/Jobs'+tag
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

print("Total Number of Jobs ",nJobsTobeMade)
if not args.dryRun:
    print(f"All condor submit files to be submitted , {len(allCondorSubFiles)} jobs")
    for fle in allCondorSubFiles:
        print('condor_submit '+fle)
        if submit2Condor or resubmit2Condor:
            os.system('condor_submit '+fle)
print("== "*20)
print("  CLUB HASH ",club_hash)
print("== "*20)
print("")
