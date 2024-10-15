#!/usr/bin/env python 
from __future__ import print_function
import os,glob
import copy,json
import argparse
from correctionlib import _core


condorScriptString="\
executable = $(filename)\n\
output = $Fp(filename)run.$(Cluster).stdout\n\
error = $Fp(filename)run.$(Cluster).stderr\n\
log = $Fp(filename)run.$(Cluster).log\n\
request_cpus = 2\n\
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
JOB_TYPE='mva_tuple'

parser = argparse.ArgumentParser()
parser.add_argument('-s',"--submit", help="Submit file to condor pool", action='store_true' )
parser.add_argument('-r',"--resubmit", help="Re-Submit file to condor pool", action='store_true' )
parser.add_argument('-t',"--test", help="Test Job", action='store_true' )
parser.add_argument('-p',"--printOnly", help="Only Print the commands", action='store_true' )
parser.add_argument(     "--doNominal", help="Submit only Nominal Ones ! ", action='store_true' )
parser.add_argument(     "--doSystOnly", help="Submit only Systematics jobs Ones ! ", action='store_true' )
parser.add_argument('-n',"--njobs", help="Number of jobs to make",default='-6000')
parser.add_argument('-e',"--nevts", help="Number of events per job",default='-1')
parser.add_argument(     "--procs", help="Process to process ",default=None)
parser.add_argument(     "--systs", help="Systematics to process ",default=None)
parser.add_argument('-c',"--cats", help="Categories of jobs to process ",default=None)
parser.add_argument('-f',"--flist", help="Files to process",default=None)
parser.add_argument('-y',"--years", help="Year of jobs to process ",default="all")
parser.add_argument('-v',"--version", help="Vesion of the specific work",default='TEST')

args = parser.parse_args()

version=args.version
cats=args.cats
njobs=int(args.njobs)
maxevents=int(args.nevts)
max_meterialize=250
jobsToProcess=None
submit2Condor=args.submit
resubmit2Condor=args.resubmit
isTest=args.test
onlyPrint=args.printOnly

years=args.years.split(",")
job_hash=f'mvaTuple_{args.version}'

treeInMap={
    'kappaScan' : 'trees/ggHHH_125_13TeV' ,
    'ggHHH'  : 'trees/ggHHH_125_13TeV',
    'sig'  : 'trees/ggHHH_125_13TeV',
    'ggH'  : 'trees/ggHHH_125_13TeV',
    'vH'   : 'trees/ggHHH_125_13TeV',
    'vbfH' : 'trees/ggHHH_125_13TeV',
    'ttH'  : 'trees/ggHHH_125_13TeV',
    'ggHH' : 'trees/ggHHH_125_13TeV',
    'vbfHH': 'trees/ggHHH_125_13TeV',
    'ttWToQQHTo2G.root' : 'trees/ggHHH_125_13TeV',
    'ggZTo2BHHTo2B2G' : 'trees/ggHHH_125_13TeV',
    'WToQQHHTo2B2G'   : 'trees/ggHHH_125_13TeV',
    'ZToBBHHTo2B2G'   : 'trees/ggHHH_125_13TeV',
    'ttHHTo2B2G'      : 'trees/ggHHH_125_13TeV',
    'bkg' : 'trees/bkg_13TeV_TrippleHTag_0',
    'qcd' : 'trees/bkg_13TeV_TrippleHTag_0',
    'ttX' : 'trees/bkg_13TeV_TrippleHTag_0',
    'data' : 'trees/Data_13TeV_TrippleHTag_0'
}
treeOutMap={
    'kappaScan' : 'trees/ggHHH_125_13TeV' ,
    'ggHHH' : 'trees/ggHHH_125_13TeV' ,
    'sig' : 'trees/ggHHH_125_13TeV' ,
    'singleH' : 'trees/singleH_125_13TeV',
    'doubleH' : 'trees/doubleH_125_13TeV',
    'bkg' : 'trees/bkg_13TeV_TrippleHTag_0',
    'qcd' : 'trees/bkg_13TeV_TrippleHTag_0',
    'ttX' : 'trees/bkg_13TeV_TrippleHTag_0',
    'data' : 'trees/Data_13TeV_TrippleHTag_0'
}

procToProcess=[]
print(" submit jobs ",submit2Condor)
print(" resubmit jobs ",resubmit2Condor)
print(" isTest ",isTest)
print(" printOnly ",onlyPrint)
print(" njobs ",njobs)
print(" maxEvt ",maxevents)
print(" cats ",cats)
print(" years ",years)

if args.procs:
    procToProcess=args.procs.split(",")   
if submit2Condor or resubmit2Condor:
    choice=input("Do you really want to submit the jobs to condor pool ? ")
    if 'y' not in choice.lower():
        print("Exiting ! ")
        exit(0)

jobDict={}
jobConfigJson='misc/jsons/exportToMVATuple.json'
with open(jobConfigJson) as f:
    jobDict=json.load(f)


fbase='/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/results/analysisSystNtuples/v30p1_22MarFriday09PM59m33s/v30p1/doubleH/ggHH*/*/*/*.root'
fbase='/home/athachay/t3store3/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/results/analysisSystNtuples/v26p1/*kappa*/'
fbase='/home/athachay/t3store3/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/results/analysisSystNtuples/v26p*/*/*/*/*/*'
fbase=[]
fbase+=['/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/results/analysisSystNtuples/v29p2_07JanSunday05PM00m19s/v29p2/*/*/*/*/*.root'] # baseline
fbase+=['/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/results/analysisSystNtuples/v30_29FebThursday03PM28m34s/v30/sig/*/*/*/*.root'] # 2016Pre
fbase+=['/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/results/analysisSystNtuples/v30p1_08MarFriday06AM56m34s/v30p1/*/*/*/*/*.root'] # kappaScan Full UL
fbase+=['/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/results/analysisSystNtuples/v30p1_22MarFriday09PM59m33s/v30p1/*/*/*/*/*.root'] # ttHH , ggHH , bbH
fbase =['/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/results/analysisSystNtuples/v33p0_09AprTuesday04AM54m30s/v33p0/v33p0/*/*/*/*/*.root'] # nominal with btagScales
fbase =['/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/results/analysisSystNtuples/v33p1_09AprTuesday10PM23m58s/v33p1/v33p0/*/*/*/*/*.root']
fbase+=['/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/results/analysisSystNtuples/v33p0_09AprTuesday04AM54m30s/v33p0/v33p1/*/*/*/*/*.root']

### NEEDED To Correct the double Counting 
fbase+=['/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/results/analysisSystNtuples/36p0_06JunThursday07PM31m23s/36p0/v33p0/*/*/*/*/*.root']

## Missing c3-d4 (0,-1) files
fbase+=['/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/results/analysisSystNtuples/v40p0_17SepTuesday07PM40m20s/v40p0/v40p0/*/*/*/*/*.root']


print("Looking for file ! ")
allFiles=[]
if args.flist:
    with open(args.flist) as f:
        txt=f.readlines()
        for l in txt:
            allFiles.append(l[:-1])
else:        
    #allFiles=glob.glob(fstring)
    for fstring in fbase:
        print(f"Search string : {fstring}")
        fls_=glob.glob(fstring)
        print("  > number of files got : ",len(fls_))
        fls=[]
        for fl in fls_:
            if 'depricated' in fl:
                continue
            fls.append(fl)
        print("  > number of files after removing depricated files : ",len(fls))
        allFiles+=fls

print(f"Obtained : {len(allFiles)} Files")
fileMap={}
for fl in allFiles:
    items=fl.split("/")
   # print(items)
    tag=items[-5]
    proc=items[-4]
    syst=items[-3]
    era =items[-2]
    fileName=items[-1]
    # print(f"{tag=} {proc=} {syst=} {era=} {fileName=}")
    if 'all' not in years:
        if era not in years:
            continue
    if tag  not in fileMap:  
        fileMap[tag]={}  ; 
        print(f"Adding tag : {tag}")
    if proc not in fileMap[tag] : 
        fileMap[tag][proc]={} ; 
        print(f"Adding proc : {proc} for {tag}")
    if syst  not in fileMap[tag][proc] : 
        fileMap[tag][proc][syst]={} ;  
        #print(f"Adding syst : {syst} for {tag}/{proc}")
    if era not in fileMap[tag][proc][syst] : 
        fileMap[tag][proc][syst][era]=[] ;  
        #print(f"Adding era : {era} for {tag}/{proc}/{syst}") 
    #fileMap[tag][proc][syst][era].append(fl)
    fileMap[tag][proc][syst][era]=[fl]
for tag in fileMap:
    print(f"Got {tag} with {len(fileMap[tag])} procs")
    for proc in fileMap[tag]:
        n=0
        eras=set()
        for syst in fileMap[tag][proc]:
            n+=1
            for er in fileMap[tag][proc][syst]:
                eras.add(er)
        print(f"  > Got {proc} with {len(fileMap[tag][proc])} systs and {n} files total ! {eras} ")

with open("out.json",'w') as f:
    json.dump( fileMap ,f, indent=4  )


templateCMD="""
./misc/condorJobMakerJsonStyle.py 
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
#jobsToProcess=['sig']
#jobsToProcess=['data']
targetYear='run2'
if cats:
    jobsToProcess=cats.split(",")

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
    
    cfgScriptTemplate=jobDict[htag]['cfg']
    configurationTxt=[]
    with open(cfgScriptTemplate,'r') as f:
        configurationTxt=f.readlines()
    configurationTxt=''.join(configurationTxt)
    cutFile = jobDict[htag]['cuts']
    catFile = jobDict[htag]['cats']
    scaleFile = jobDict[htag]['scales']
    executable=jobDict[htag]["exe"]
    hasAJob=False
    for proc in fileMap[htag]:
        if len(procToProcess):
            if proc not in procToProcess:
                print("skipping  ",proc)
                continue

        head=pwd+f'/Condor/{JOB_TYPE}/{job_hash}/{htag}/{proc}/'
        print(f"Making Jobs for {htag}/{proc}")
        njobs=0
        treeIn  = 'bkg_13TeV_TrippleHTag_0'
        if htag in treeInMap:
            treeIn = treeInMap[htag]
        if proc in treeInMap:
            treeIn  = treeInMap[proc]
        treeOut = treeOutMap[htag]
        print(f"\t Setting {treeIn=}")
        print(f"\t Setting {treeOut=}")
        for syst in fileMap[htag][proc]:
            if args.doNominal and (syst!='nominal'):
                continue
            if args.doSystOnly and (syst=='nominal'):
                continue
            systname = syst
            if args.doSystOnly:
                pass
                #if ('JEC' not in syst):
                    #print("   > skipping ",syst)
                #    continue
            if args.systs:
                if syst not in args.systs.split(","):
                    print("   > skipping ",syst)
                    continue
            #if 'Y2016' in systname : systname=systname.replace('Y2016','Yearly')
            #if 'Y2017' in systname : systname=systname.replace('Y2017','Yearly')
            #if 'Y2018' in systname : systname=systname.replace('Y2018','Yearly')
            
            #if ('2016' not in systname) and ('2017' not in systname) and ( '2018' not in systname ) :
            #    continue

            #if '2016'  in systname : systname=systname.replace('2016' ,'Yearly')
            #if '2017'  in systname : systname=systname.replace('2017' ,'Yearly')
            #if '2018'  in systname : systname=systname.replace('2018' ,'Yearly')
            
            #print("Doing Syst ",syst,systname)


            print(f"\r\t Processing Syst : {syst}                                            ",end="")
            #print(f"\t Processing Syst : {syst}                                            ")
            for era in fileMap[htag][proc][syst]:
                i+=1
                njobs+=1
                if onlyPrint:
                    continue
                fileNames=','.join(fileMap[htag][proc][syst][era])
                dirName =f'{head}/Job_{era}_{syst}/'
                if not os.path.exists(dirName):
                    os.system('mkdir -p '+dirName)
                destination=f'{RESULT_BASE}/{JOB_TYPE}/{job_hash}/{htag}/{proc}/{syst}/{era}/'
                if resubmit2Condor:
                    continue
                if not os.path.exists(destination):
                    os.system('mkdir -p '+destination)
                cfgFileName='customize_'+str(i)+'.json'
                cfgFile=open(dirName+'/'+cfgFileName,'w')
                tmp=''
                tmp=configurationTxt.replace("@@FNAME",fileNames)
                tmp=tmp.replace("@@TREEIN" ,treeIn  )
                
                treeOutForSyst = treeOut
                if syst=='nominal':
                    tmp=tmp.replace("@@TREEOUT",treeOutForSyst )
                else:    
                    tmp=tmp.replace("@@TREEOUT",treeOutForSyst+f'_{systname}' )

                tmp=tmp.replace("@@CATFILE",catFile)
                tmp=tmp.replace("@@CUTFILE",cutFile)
                tmp=tmp.replace("@@SCALEFILE",scaleFile)
                tmp=tmp.replace("@@PTAG",htag)
                tmp=tmp.replace("@@IDX",str(i))
                tmp=tmp.replace("@@PROC" , str(proc))
                tmp=tmp.replace("@@ERA"  , str(era ))
                tmp=tmp.replace("@@TARGET_YEAR"  , str(targetYear))
                tmp=tmp.replace("@@SYST" , str(syst))
                tmp=tmp.replace("@@MAXEVENTS",str(NEVENTS_PER_JOB))
                tmp=tmp.replace("@@DESTINATION",destination)
                cfgFile.write(tmp)
                cfgFile.close()   
 
                runScriptName=dirName+f'/{htag}_{proc}_{syst}_{era}_run.sh'
                if os.path.exists(runScriptName+'.sucess'):
                   os.system('rm '+runScriptName+'.sucess')
                runScript=open(runScriptName,'w')
                tmp=runScriptTxt.replace("@@DIRNAME",dirName)
                tmp=tmp.replace("@@IDX",str(i))
                tmp=tmp.replace("@@CFGFILENAME",cfgFileName)
                tmp=tmp.replace("@@RUNSCRIPT",runScriptName)
                tmp=tmp.replace("@@EXECUTABLE",executable)
                tmp=tmp.replace("@@proxy_path",proxy_path)
                tmp=tmp.replace("@@HOME",HOME)
                tmp=tmp.replace("@@DESTINATION",destination)
                runScript.write(tmp)
                runScript.close()
                os.system('chmod +x '+runScriptName)
                hasAJob=True
            if isTest:
                break

        print()
        if hasAJob:
            condorScriptName=head+'/job'+htag+'.sub'
            if not onlyPrint:
                with open(condorScriptName,'w') as condorScript:
                    condorScript.write(condorScriptString)
                    if args.doNominal:
                        condorScript.write("queue filename matching ("+head+"/*nominal*/*.sh)\n")
                    else:
                        condorScript.write("queue filename matching ("+head+"/*/*.sh)\n")
                    print(f"Condor {njobs} Jobs made !\n\t submit file  : {condorScriptName}")
            allCondorSubFiles.append(condorScriptName)
            print("Made ",njobs," jobs for ",proc)
        njobs_total+=njobs
        if resubmit2Condor:
            continue
        if isTest:
            break
print("")
print("")
print("All condor submit files to be submitted ")
for fle in allCondorSubFiles:
    print('condor_submit '+fle)
    if submit2Condor or resubmit2Condor:
        os.system('condor_submit '+fle)
print("Total number of jobs ",njobs_total)
print("")
