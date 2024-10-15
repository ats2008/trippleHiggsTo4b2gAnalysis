#!/usr/bin/env python 
import os
import copy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-n',"--njobs", help="Number of jobs to make",default='10000')
parser.add_argument('-e',"--nevts", help="Number of events per job",default='-1')
parser.add_argument('-v',"--version", help="Vesion of the specific work",default='1p0')

args = parser.parse_args()

version=args.version

exe='/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntupleMaker_BDTScoreAdder.py'
destn='results/BDT_NTuples/'
cfgTemplate='misc/ntuplizerForML_Signal.tpl.cfg'
scriptTemplate='misc/runPython.tpl.sh'
njobs=int(args.njobs)
maxevents=args.nevts
max_meterialize=250
jobsToProcess=None

jobDict={}
jobDict['bdtScoreAdder']={}
jobDict['bdtScoreAdder']['exe']=exe
jobDict['bdtScoreAdder']['cfg']='misc/cfg/ntuplizer_BDTscoreAdder_sig.tpl.cfg'
jobDict['bdtScoreAdder']['script']=scriptTemplate
jobDict['bdtScoreAdder']['destn']=destn
jobDict['bdtScoreAdder']['njobs']=njobs
jobDict['bdtScoreAdder']['files_per_job']=1
jobDict['bdtScoreAdder']['maxevents']=maxevents
jobDict['bdtScoreAdder']['max_meterialize']=max_meterialize
jobDict['bdtScoreAdder']['datasets']={}
jobDict['bdtScoreAdder']['datasets']['ggHHH2018'] = 'fileList/ntuples/ggHHH_UL2018_mlScoreUpdated.fls'

jobDict['bdtScoreAdderBkg']=copy.copy( jobDict['bdtScoreAdder'] )
jobDict['bdtScoreAdderBkg']['cfg']='misc/cfg/ntuplizer_BDTscoreAdder_bkg.tpl.cfg'
jobDict['bdtScoreAdderBkg']['datasets']={}
jobDict['bdtScoreAdderBkg']['datasets']['ggJBox1B']  = 'fileList/ntuples/DiPhotonJetsBox1BJet_MGG-80toInf_13TeV_mlScoreUpdated.fls'
jobDict['bdtScoreAdderBkg']['datasets']['ggJBox2B']  = 'fileList/ntuples/DiPhotonJetsBox2BJets_MGG-80toInf_13TeV_mlScoreUpdated.fls'
jobDict['bdtScoreAdderBkg']['datasets']['ggJBox']    = 'fileList/ntuples/DiPhotonJetsBox_MGG-80toInf_13TeV_mlScoreUpdated.fls'

jobDict['bdtScoreAdderData']=copy.copy( jobDict['bdtScoreAdder'] )
jobDict['bdtScoreAdderData']['cfg']='misc/cfg/ntuplizer_BDTscoreAdder_data.tpl.cfg'
jobDict['bdtScoreAdderData']['datasets']={}
jobDict['bdtScoreAdderData']['datasets']['2018A']  = 'fileList/ntuples/data2018A_v2_mlupd.fls'
jobDict['bdtScoreAdderData']['datasets']['2018B']  = 'fileList/ntuples/data2018B_v2_mlupd.fls'
jobDict['bdtScoreAdderData']['datasets']['2018C']  = 'fileList/ntuples/data2018C_v2_mlupd.fls'
jobDict['bdtScoreAdderData']['datasets']['2018D']  = 'fileList/ntuples/data2018D_v2_mlupd.fls'

jobsToProcess=['bdtScoreAdderData']

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
       --maxMeterialize  @@MAX_METERIALIZE
"""

allCondorSubFiles=[]
if jobsToProcess==None:
    jobsToProcess=list( jobDict.keys() )
for jobTag in jobsToProcess:
    cmd=templateCMD.replace("@@EXE",jobDict[jobTag]['exe'])
    cmd=cmd.replace("@@CFG_TPL"        ,jobDict[jobTag]['cfg'])
    cmd=cmd.replace("@@SCRIPT_TPL"     ,jobDict[jobTag]['script'])
    cmd=cmd.replace("@@NJOBS"          ,str(jobDict[jobTag]['njobs']))
    cmd=cmd.replace("@@FILES_PER_JOB"        ,str(jobDict[jobTag]['files_per_job']))
    cmd=cmd.replace("@@MAXEVENTS"        , str(jobDict[jobTag]['maxevents']) )
    cmd=cmd.replace("@@MAX_METERIALIZE"        , str(jobDict[jobTag]['max_meterialize']) )
    for key in jobDict[jobTag]:
        if '@@' in key:
            cmd=cmd.replace( key ,jobDict[jobTag][ key ]  )

    for dset in jobDict[jobTag]['datasets']:
        tag=jobTag+'_'+dset+'_'+version
        cmdD=cmd.replace("@@TAG",tag)
        cmdD=cmdD.replace("@@FILELIST",jobDict[jobTag]['datasets'][dset])
        destination = jobDict[jobTag]['destn'] +'/'+tag+'/'
        cmdD=cmdD.replace("@@DESTN",destination)
        print("")
        print("")
        print("="*40)
        print(cmdD)
        print(" - -"*10)
        print("")
        print("")
        cmdD=cmdD.replace("\n","")
        os.system( cmdD )
        head='Condor/Jobs'+tag
        condorScriptName=head+'/job'+tag+'.sub'
        allCondorSubFiles.append(condorScriptName)
print("")
print("")
print("All condor submit files to be submitted ")
for fle in allCondorSubFiles:
    print('csub '+fle)
print("")
print("")
