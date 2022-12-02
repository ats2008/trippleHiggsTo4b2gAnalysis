#!/usr/bin/env python3
import os
import sys
import argparse

version='v1'

useStr="\
Usage\n\
    ./jobMakerForCmsRuns.py <executable> <InputFileListFname> <destination> <NJOBS> <FILES_PER_JOB> <jobPrefix>\n\
\n\
"

executable='recoStepUL2018.py'

NJOBS=20000
NEVENTS_PER_JOB = -1
ZERO_OFFSET=0
FILES_PER_JOB=1
maxMeterialize=100
offsetStep=10000

pwd=os.environ['PWD']
proxy_path=os.environ['X509_USER_PROXY']
HOME=os.environ['HOME']
xrdRedirector="root://cms-xrd-global.cern.ch/"

FileSource ="bmm5FileList.txt"
destination='/grid_mnt/t3storage3/athachay/bs2mumug/run2studies/CMSSW_10_6_19_patch2/src/BsMMGAnalysis/MergeWithBMMNtuples/RunLumiEventFileMaker/runLumiList/'
tag=""



parser = argparse.ArgumentParser()
parser.add_argument("--exe", help="Executable")
parser.add_argument("--fsrc", help="file source")
parser.add_argument("--runScript", help="template runscript")
parser.add_argument("--cfg", help="Configuration template file")
parser.add_argument("--dest", help="destination",default='./')
parser.add_argument("--jn", help="Max number of Jobs",default='10')
parser.add_argument("--fn", help="File per Job",default=1)
parser.add_argument("--maxEvt", help="Mac Events per Job , -1 for all",default=100)
parser.add_argument("--tag", help="Job Tag",default="")
parser.add_argument("--offsetStep", help="offsetStep for event number re-defenition",default="5000")
parser.add_argument("--maxMeterialize", help="Executable",default=None)

args = parser.parse_args()

executable=args.exe
FileSource=args.fsrc
runScriptTemplate=args.runScript
cfgScriptTemplate=args.cfg
destination=args.dest
NJOBS=int(args.jn)
FILES_PER_JOB=int(args.fn)
NEVENTS_PER_JOB=int(args.maxEvt)
tag=args.tag
maxMeterialize=int(args.maxMeterialize)
offsetStep=int(args.offsetStep)

#if len(sys.argv) > 2:
#    executable=sys.argv[1]
#    FileSource=sys.argv[2]  
#    runScriptTemplate=sys.argv[3]  
#    cfgScriptTemplate=sys.argv[4]  
#else:
#    print(useStr)
#    sys.exit(1)
#argC=5
#if len(sys.argv) > argC:
#    destination=sys.argv[argC]  
#argC+=1
#if len(sys.argv) > argC :
#    NJOBS=int(sys.argv[argC])  
#argC+=1
#if len(sys.argv) > argC :
#    FILES_PER_JOB=int(sys.argv[argC])  
#argC+=1
#if len(sys.argv) > argC :
#    NEVENTS_PER_JOB=int(sys.argv[argC])  
#argC+=1
#if len(sys.argv) > argC :
#    tag=sys.argv[argC]  
#argC+=1
#if len(sys.argv) > argC :
#    maxMeterialize=int(sys.argv[argC])  
#argC+=1
#if len(sys.argv) > argC :
#    offsetStep=int(sys.argv[argC])  

if(not os.path.exists(destination)):
    os.system("mkdir -p "+destination)

destination=os.path.abspath(destination)


print(" Executable ",executable)
print("Source file list ",FileSource)
print("Run Script Template ",runScriptTemplate)
print("Configuration ",cfgScriptTemplate )
print("destination : ",destination)
print("NJOBS : ",NJOBS)
print("FILES_PER_JOB : ",FILES_PER_JOB)
print("NEVENTS_PER_JOB : ",NEVENTS_PER_JOB)
print("maxMeterialize : ",maxMeterialize)
print("tag : ",tag)
Fnames=open(FileSource,'r')
sourceFileList=Fnames.readlines()
Fnames.close()
print("Number avilable files = ",len(sourceFileList))


condorScriptString="\
executable = $(filename)\n\
output = $Fp(filename)run.$(Cluster).stdout\n\
error = $Fp(filename)run.$(Cluster).stderr\n\
log = $Fp(filename)run.$(Cluster).log\n\
+JobFlavour = \"longlunch\"\n\
"
if maxMeterialize >0 :
    condorScriptString+="max_materialize="+str(maxMeterialize)+"\n"

f=open(runScriptTemplate,'r')
runScriptTxt=f.readlines()
f.close()
runScriptTxt=''.join(runScriptTxt)
runScriptTxt=runScriptTxt.replace("@@proxy_path",proxy_path)
runScriptTxt=runScriptTxt.replace("@@HOME",HOME)
runScriptTxt=runScriptTxt.replace("@@DESTINATION",destination)
#print(runScriptTxt)

f=open(cfgScriptTemplate,'r')
configurationTxt=f.readlines()
f.close()
configurationTxt=''.join(configurationTxt)
#print(configurationTxt)


head='Condor/Jobs'+tag
if not os.path.exists(head ):
    os.system('mkdir -p '+head)

condorScriptName=head+'/job'+tag+'.sub'
condorScript=open(condorScriptName,'w')
condorScript.write(condorScriptString)

n=int(len(sourceFileList)/FILES_PER_JOB) + 1
if n < NJOBS:
    NJOBS=n
print("Making ",NJOBS," Jobs ")

njobs=0
for ii in range(NJOBS):
    i=ii+ZERO_OFFSET
    
    if len(sourceFileList)<FILES_PER_JOB:
       print("\nfname count less than required .. stoping ")
       FILES_PER_JOB=len(sourceFileList)
    
    if len(sourceFileList) ==0:
       break 

    dirName= pwd+'/'+head+'/Job_'+str(i)
    
    if(ii%10==0) : print("\nJob Made : ",end = " " )
    print(ii,end =" ")

    if os.path.exists(dirName):
        k=True
    else:
        os.system('mkdir '+dirName)
    
    cfgFileName='customize_'+str(i)+'.cfg'
    cfgFile=open(dirName+'/'+cfgFileName,'w')
    tmp=''
    for j in range(FILES_PER_JOB):
      tmp+=sourceFileList.pop(0)[:-1]+"\n"
    tmp=configurationTxt.replace("@@FNAMES",tmp[:-1])
    tmp=tmp.replace("@@TAG",tag)
    tmp=tmp.replace("@@IDX",str(i))
    tmp=tmp.replace("@@OFFSET_IDX",str(i*offsetStep))
    tmp=tmp.replace("@@MAXEVENTS",str(NEVENTS_PER_JOB))
    cfgFile.write(tmp)
    cfgFile.close()   
    
    runScriptName=dirName+'/'+tag+'run'+str(i)+'.sh'
    if os.path.exists(runScriptName+'.sucess'):
       os.system('rm '+runScriptName+'.sucess')
    runScript=open(runScriptName,'w')
    tmp=runScriptTxt.replace("@@DIRNAME",dirName)
    tmp=tmp.replace("@@IDX",str(i))
    tmp=tmp.replace("@@CFGFILENAME",cfgFileName)
    tmp=tmp.replace("@@RUNSCRIPT",runScriptName)
    tmp=tmp.replace("@@EXECUTABLE",executable)
    runScript.write(tmp)
    runScript.close()
    os.system('chmod +x '+runScriptName)
    if maxMeterialize < 0:
        condorScript.write("queue filename matching ("+runScriptName+")\n")
    njobs+=1
if maxMeterialize >=0:
    condorScript.write("queue filename matching ("+head+"/*/*.sh)\n")
    
print()
print(" Number of jobs made : ", njobs)
print(" Number of files left : ", len(sourceFileList) )
print(" Condor submit file  : ", condorScriptName)
condorScript.close()
