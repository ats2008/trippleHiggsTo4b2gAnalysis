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
nMatchFound=0
fname=allFnames[0]
simFile = ROOT.TFile(fname,'READ')
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

ntuple={}
outputDataDict=OrderedDict(zip(branches, [np.nan]*len(branches))) 
ntuple['all'] = ROOT.TNtuple(outTreeName, outTreeName, ':'.join(branches))


kyList = [ ky for ky in tofill ]
print("len(branches) : " , len(branches))


m1m2 = ROOT.TH2F("m1jj_m2jj","H1bb , H2bb mass",300,0.0,300.0,300,0.0,300. )


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
        eTree=simFile.Get('tagsDumper/trees/Data_13TeV_TrippleHTag_0')
    maxEvents_ = eTree.GetEntries()
    if(maxEvents >0  and (totalEvents+maxEvents_) > maxEvents):
        maxEvents_= (maxEvents - totalEvents)
    totalEvents+=maxEvents_
    allBranches=[]
    for ky in eTree.GetListOfBranches():
        allBranches.append(ky.GetName())
    for evt in range(maxEvents_):
        eTree.GetEntry(evt)
        if(evt%100==0):
            print("   Doing i = ",evt," / ",maxEvents_,
                  " n Found : ", nMatchFound,
                  " n Miss  : ", nMiss)
            print(" gg mass , h1MassVsh2Mass  : ",eTree.CMS_hgg_mass ,"( ",eTree.M1jj , eTree.M2jj,")" )
            
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
        if(len(branchesToFill)!=len(outputDataDict)):
            print("len(branchesToFill)!=len(outputDataDict) !! ")
            exit(1)
        cat='cat0'
        ntuple['all'].Fill(array('f', outputDataDict.values()))

    simFile.Close()           
    print("Closing file : ",fname)
print("Number of evetsProcess ",nMatchFound)

dir_.cd()    
sumWeights.Write()
m1m2.Write()

for cat in ntuple:
    ntuple[cat].Write()
fout.Close()
print(" File written out  : ",foutName)

