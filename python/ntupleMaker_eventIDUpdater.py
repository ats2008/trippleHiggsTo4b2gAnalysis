from __future__ import print_function
from collections import OrderedDict
import ROOT 
import numpy as np
from trippleHiggsUtils import *
from Util  import *
from branches import *
from array import array
from trippleHiggsSelector import *
from TMVA_Model import *

import  trippleHiggsMLInterface as hhhMLI

import os,sys

NJETS=8

cfgFileName=''
if len(sys.argv) <2:
    print("Usage\n\t ~$python recoAnalyzer.py <configFile>\n")
    exit(1)
else:
    cfgFileName=sys.argv[1]
maxEvtsSuperSeeder=-1
if len(sys.argv) >2:
    maxEvtsSuperSeeder=int(sys.argv[2])
print("Loading cfg file ",cfgFileName)
f=open(cfgFileName,'r')
cfgTxt=f.readlines()
f.close()

headers=getListOfStringsFromConfigs(cfgTxt,"#HEADER_BEG","#HEADER_END")
for header in headers:
    print("Loading cfg file ",header)
    f=open(header,'r')
    tmp=f.readlines()
    f.close()
    for l in tmp:
        cfgTxt.append(l)

def insideTheEllipse( x,y,x1,y1,x2,y2,a):
    return np.sqrt( (x1-x)*(x1-x)+(y1-y)*(y1-y) ) + np.sqrt( (x2-x)*(x2-x)+(y2-y)*(y2-y) ) < a           

allFnames=getListOfStringsFromConfigs(cfgTxt,"#FNAMES_BEG","#FNAMES_END")
foutName=getValueFromConfigs(cfgTxt,"OutpuFileName","fggHists.root")
processID=getValueFromConfigs(cfgTxt,"processID",default="DATA")
treeName=getValueFromConfigs(cfgTxt,"treeName",default="tagsDumper/trees/Data_13TeV_TrippleHTag_0")
outTreeName=getValueFromConfigs(cfgTxt,"outTreeName",default="Data_13TeV_TrippleHTag_0")
eventIdxOffset=int(getValueFromConfigs(cfgTxt,"EventIdxOffset",default=0))

etaMax=2.5
pTMin=25.0

print("allFnames   :  ",              allFnames)
print("foutName   :  ",               foutName)
print("processID   :  ",              processID)
print("treeName   :  ",               treeName)
print("eventIdxOffset   :  ", eventIdxOffset)

maxEvents=-1
tmp_=getValueFromConfigs(cfgTxt,"MaxEvents")
if tmp_!='':
    maxEvents=int(tmp_)

for i in allFnames:
    print(" file : ",i)
if(maxEvtsSuperSeeder > 0):
    maxEvents=maxEvtsSuperSeeder
print("maxevents : ",maxEvents)


nNoHtoGammaGamma=0
nNoHHto4B=0
totalEvents=0
nMiss=0
nProcessed=0
fname=allFnames[0]
simFile = ROOT.TFile(fname,'READ')
eTree=simFile.Get(treeName)
if not eTree:
   treeName='tagsDumper/trees/EGamma_13TeV_TrippleHTag_0'
   print(" Trying tree name  : ",treeName )
   eTree=simFile.Get(treeName)
allBranchesFromFile=[]
for ky in eTree.GetListOfBranches():
    allBranchesFromFile.append(ky.GetName())

branchesToFill=allBranchesFromFile
branchesToFill.append("eventIdx")

branches=np.unique(branchesToFill)

filemode="RECREATE"
fout = ROOT.TFile(foutName, filemode)
dir_ = fout.mkdir('tagsDumper')
dir_ = dir_.mkdir('trees')

tofill = OrderedDict(zip(branches, [np.nan]*len(branches)))
outputDataDict=OrderedDict(zip(branches, [np.nan]*len(branches))) 

ntuple={}
ntuple['all'] = ROOT.TNtuple(outTreeName, outTreeName, ':'.join(branches))


kyList = [ ky for ky in tofill ]
print("len(branches) : " , len(branches))




sumWeights=ROOT.TH1F("sumEvts","sumEvts",1,0.0,1.0)
sumWeights.SetCanExtend(ROOT.TH1.kAllAxes)
th1Store={}

eventIdx=eventIdxOffset
for fname in allFnames:
    print("Opening file : ",fname)
    simFile = ROOT.TFile(fname,'READ')
    eTree=simFile.Get(treeName)
    print(" NEntries = ", eTree.GetEntries() )
    if not eTree:
        treeName='tagsDumper/trees/EGamma_13TeV_TrippleHTag_0'
        print(" Trying tree name  : ",treeName )
        eTree=simFile.Get(treeName)
    maxEvents_ = eTree.GetEntries()
    if(maxEvents >0  and (totalEvents+maxEvents_) > maxEvents):
        maxEvents_= (maxEvents - totalEvents)
    totalEvents+=maxEvents_
    allBranches=[]
    for ky in eTree.GetListOfBranches():
        allBranches.append(ky.GetName())
    for evt in range(maxEvents_):
        eTree.GetEntry(evt)
        sumWeights.Fill('Events',1)
        sumWeights.Fill('weights',eTree.weight)
        if(evt%1000==0):
            print("   Doing i = ",evt," / ",maxEvents_," n Found : ", nProcessed)
            print("     gg mass , h1MassVsh2Mass  : ",eTree.CMS_hgg_mass )
            
        wei=eTree.weight
        for ky in tofill:
            if ky in allBranches:
                tofill[ky]=getattr(eTree,ky)
            else:
                tofill[ky]=-999e3
        
        tofill['eventIdx']=eTree.event
        tofill['event']=eventIdx
        eventIdx+=1
       
        for i in outputDataDict:
            outputDataDict[i]=tofill[i]
        if(len(branches)!=len(outputDataDict)):
            print("len(branches)!=len(outputDataDict) !! ",len(branches), " != ",len(outputDataDict))
            for i in branches:
                if i not in outputDataDict:
                    print(i,"in branches not in outputDataDict")
            for i in outputDataDict:
                if i not in branches:
                    print(i,"in outputDataDict not in branches")
            print("EXITING") 
            exit(1)
        cat='cat0'
        ntuple['all'].Fill(array('f', outputDataDict.values()))
        nProcessed+=1
    simFile.Close()           
    print("Closing file : ",fname)
print("Number of evets Processed ",nProcessed)

dir_.cd()    
sumWeights.Write()

for cat in ntuple:
    ntuple[cat].Write()
fout.Close()
print(" File written out  : ",foutName)

