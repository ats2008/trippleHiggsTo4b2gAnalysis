from __future__ import print_function
from collections import OrderedDict
import ROOT 
import numpy as np
from array import array

import Util as utl

import trippleHiggsMLInterface as hhhMLI
import trippleHiggsUtils as hhhUtil
import trippleHiggsSelector as hhhSelector

import os,sys
import datetime

"""
    Tries to save all the required varables into  histograms for data and MC
"""


NJETS=8

cfgFileName=''
if len(sys.argv) <2:
    print("Usage\n\t ~$python SCRPPT.py <configFile>\n")
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

headers=utl.getListOfStringsFromConfigs(cfgTxt,"#HEADER_BEG","#HEADER_END")
for header in headers:
    print("Loading cfg file ",header)
    f=open(header,'r')
    tmp=f.readlines()
    f.close()
    for l in tmp:
        cfgTxt.append(l)

def insideTheEllipse( x,y,x1,y1,x2,y2,a):
    return np.sqrt( (x1-x)*(x1-x)+(y1-y)*(y1-y) ) + np.sqrt( (x2-x)*(x2-x)+(y2-y)*(y2-y) ) < a           

allFnames =utl.getListOfStringsFromConfigs(cfgTxt,"#FNAMES_BEG","#FNAMES_END")
foutName  =utl.getValueFromConfigs(cfgTxt,"OutpuFileName","fggHists.root")
processID =utl.getValueFromConfigs(cfgTxt,"processID",default="DATA")
treeName  =utl.getValueFromConfigs(cfgTxt,"treeName",default="tagsDumper/trees/Data_13TeV_TrippleHTag_0")
doOverlapRemoval =utl.getValueFromConfigs(cfgTxt,"doOverlapRemoval",default="1") ; doOverlapRemoval = int(doOverlapRemoval) > 0.5
isData =utl.getValueFromConfigs(cfgTxt,"isData",default="1") ; isData = int(isData) > 0.5
etaMax =float(utl.getValueFromConfigs(cfgTxt,"etaMax",default="2.5"))
pTMin =float(utl.getValueFromConfigs(cfgTxt,"pTMin",default="25.0"))
overlapRemovalDRMax =float(utl.getValueFromConfigs(cfgTxt,"overlapRemovalDRMax",default="0.4"))

etaMax=2.5
pTMin=25.0

print("allFnames   :  ",              allFnames)
print("foutName   :  ",               foutName)
print("processID   :  ",              processID)
print("pTMin   :  ",              pTMin)
print("etaMax   :  ",             etaMax)
print("doOverlapRemoval   :  ", doOverlapRemoval)

maxEvents=-1
tmp_=utl.getValueFromConfigs(cfgTxt,"MaxEvents")
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

filemode="RECREATE"
fout = ROOT.TFile(foutName, filemode)
dir_ = fout


sumWeights=ROOT.TH1F("sumEvts","sumEvts",1,0.0,1.0)
sumWeights.SetCanExtend(ROOT.TH1.kAllAxes)
isDataHist=ROOT.TH1F("isData","isData",1,0.0,1.0)
isDataHist.SetCanExtend(ROOT.TH1.kAllAxes)

histStore = hhhUtil.getRecoHistos()
hhhUtil.addKinematicVars(histStore)
hhhUtil.addFlasggVars(histStore)
hhhUtil.addOtherDerivedVariables(histStore)
beg=datetime.datetime.now()

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
    
    isData_= isData or 'EGamma' in fname.split('/')[-1] 
    isData_= isData_ or ('DoubleEG' in fname.split('/')[-1])
    print(" isData : ",isData_)
    isDataHist.Fill('isData',isData_)

    for ky in eTree.GetListOfBranches():
        allBranches.append(ky.GetName())
    for evt in range(maxEvents_):
        eTree.GetEntry(evt)

        isMasked = isData_ and  (eTree.CMS_hgg_mass > 115.0) and (eTree.CMS_hgg_mass < 135.0)
        if isMasked:
            continue
        wei=eTree.weight
        histStore["events"]['nEvents' ].Fill('total',1)

        sumWeights.Fill('Events',1)
        sumWeights.Fill('weights',wei)
        if(evt%1000==0):
            i=evt
            now=datetime.datetime.now()
            timeSpendSec=np.round((now-beg).total_seconds() , 2)
            timeLeftSec =np.round(1.0*(maxEvents_-i)*timeSpendSec/( i +1e-3),2)
            print("   Doing i = ",i," / ",maxEvents_)
            print("      time left : ", str(datetime.timedelta(seconds= timeLeftSec)),
                    " [ time elapsed : ",datetime.timedelta(seconds= timeSpendSec), " s ]")
            print(" gg mass   : ",eTree.CMS_hgg_mass)
        jetMask=hhhSelector.getSelectedJetCollectionMaskEta(eTree,etaMax=etaMax)
        if sum(jetMask) < 4 :
            histStore["events"]['allFailIdsEvents' ].Fill('nJetPreselectionEta',1)
            histStore["events"]['allFailIdsWeights'].Fill('nJetPreselectionEta',wei)
            continue

        jetMask=hhhSelector.getSelectedJetCollectionMaskPt(eTree,jetMask=jetMask,pTMin=pTMin)
        if sum(jetMask) < 4 :
            histStore["events"]['allFailIdsEvents' ].Fill('nJetPreselectionPt',1)
            histStore["events"]['allFailIdsWeights'].Fill('nJetPreselectionPt',wei)
            continue
        if doOverlapRemoval:
            jetMask=hhhSelector.getSelectedJetCollectionMaskOverLap(eTree,jetMask=jetMask,overlapRemovalDRMax=overlapRemovalDRMax)
            if sum(jetMask) < 4 :
                histStore["events"]['allFailIdsEvents' ].Fill('nJetPreselectionOR',1)
                histStore["events"]['allFailIdsWeights'].Fill('nJetPreselectionOR',wei)
                continue
    
        allQuads=hhhSelector.getBJetParisFGG(eTree,mask=jetMask)
        if not allQuads['isValid']:
            histStore["events"]['allFailIdsEvents'].Fill(allQuads['fail'],1)
            histStore["events"]['allFailIdsWeights'].Fill(allQuads['fail'],wei)
            continue

        quad=allQuads['bJetQuad']
        
        LVStore = hhhUtil.getLVStoreFromTreeAndQuad(eTree,quad)
        
        hhhUtil.fillTrippleHRecoVariables( eTree,histStore,LVStore,quad  )
        hhhUtil.fillKinematicVarsFromLV(LVStore,histStore["kinematicVars"],wei=wei)
        hhhUtil.fillFlashggVars(histStore,eTree,LVStore=LVStore,wei=wei)    
        varDict=hhhUtil.getOtherDerivedVariables(eTree,LVStore,quad)
        hhhUtil.fillOtherDerivedVariables(histStore,varDict,wei)

        nProcessed+=1
        histStore["events"]['nEvents' ].Fill('finalCount',1)
    simFile.Close()           
    print("Closing file : ",fname)
print("Number of evets Processed ",nProcessed)

dir_.cd()    
sumWeights.Write()
isDataHist.Write()
utl.saveToFile(histStore,fout)
fout.Close()
print(" File written out  : ",foutName)

