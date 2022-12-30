#!/usr/bin/env python 
import os
import copy

version='1p1'

exe='/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntuple_forBDT.py'
destn='results/BDT_NTuples/'
cfgTemplate='misc/ntuplizerForML_Signal.tpl.cfg'
scriptTemplate='misc/runPython.tpl.sh'
njobs=1000
max_meterialize=250

jobDict={}
jobDict['bdtNtuplizer']={}
jobDict['bdtNtuplizer']['exe']=exe
jobDict['bdtNtuplizer']['cfg']='misc/ntuplizer_BDT_sig.tpl.cfg'
jobDict['bdtNtuplizer']['script']=scriptTemplate
jobDict['bdtNtuplizer']['destn']=destn
jobDict['bdtNtuplizer']['njobs']=njobs
jobDict['bdtNtuplizer']['files_per_job']=1
jobDict['bdtNtuplizer']['maxevents']=-1
jobDict['bdtNtuplizer']['max_meterialize']=max_meterialize
jobDict['bdtNtuplizer']['datasets']={}
jobDict['bdtNtuplizer']['datasets']['ggHHH2018'] = 'fileList/ntuples/ggHHH_UL2018_mlScoreUpdated.fls'

jobDict['bdtNtuplizerBkg']=copy.copy( jobDict['bdtNtuplizer'] )
jobDict['bdtNtuplizerBkg']['cfg']='misc/ntuplizer_BDT.tpl.cfg'
jobDict['bdtNtuplizerBkg']['datasets']={}
jobDict['bdtNtuplizerBkg']['datasets']['ggJBox1B']  = 'fileList/ntuples/DiPhotonJetsBox1BJet_MGG-80toInf_13TeV_mlScoreUpdated.fls'
jobDict['bdtNtuplizerBkg']['datasets']['ggJBox2B']  = 'fileList/ntuples/DiPhotonJetsBox2BJets_MGG-80toInf_13TeV_mlScoreUpdated.fls'
jobDict['bdtNtuplizerBkg']['datasets']['ggJBox']    = 'fileList/ntuples/DiPhotonJetsBox_MGG-80toInf_13TeV_mlScoreUpdated.fls'

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
for jobTag in jobDict:
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
