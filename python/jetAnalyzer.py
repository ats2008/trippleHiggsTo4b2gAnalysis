from __future__ import print_function
import ROOT 
import numpy as np
from trippleHiggsSelector import *
from Util  import *

import os,sys

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

allFnames=getListOfStringsFromConfigs(cfgTxt,"#FNAMES_BEG","#FNAMES_END")
foutName=getValueFromConfigs(cfgTxt,"OutpuFileName")

maxEvents=-1
tmp_=getValueFromConfigs(cfgTxt,"MaxEvents")
if tmp_!='':
    maxEvents=int(tmp_)

for i in allFnames:
    print(" file : ",i)
if(maxEvtsSuperSeeder > 0):
    maxEvents=maxEvtsSuperSeeder
print("maxevents : ",maxEvents)

histStore = getRecoHistos()
histStore['jetIDStudy']= getHHHJethistos()
histStore['jetIDSelectedQuad']= getHHHJethistos()

nNoHtoGammaGamma=0
nNoHHto4B=0
totalEvents=0

isMC = True
isMC = False

for fname in allFnames:
    
    simFile = ROOT.TFile(fname,'READ')
    eTree=simFile.Get('ntuplizer/tree')
    if not eTree:
        eTree=simFile.Get('Ntuple/tree')
    maxEvents_ = eTree.GetEntries()
    if(maxEvents >0  and (totalEvents+maxEvents_) > maxEvents):
        maxEvents_= (maxEvents - totalEvents)
    totalEvents+=maxEvents_
    
    for i in range(maxEvents_):

        eTree.GetEntry(i)
        if(i%100==0):
            print("Doing i = ",i," / ",maxEvents_,
                  " n No HiggsCand : ", nNoHtoGammaGamma,
                  " n No 4BCands  : ", nNoHHto4B)
            
        histStore['events']['nEvents'].Fill("nEvents",1);
        
        if( not triggerSelection(eTree)):        continue
        histStore['events']['nEvents'].Fill("nTriggered",1);

        diPhotons_=getDiPhotons(eTree)
        
        if(not diPhotons_['valid']):
            nNoHtoGammaGamma+=1
            continue
        
        histStore['events']['nEvents'].Fill("nDiPhotonSelected",1);
        idx=diPhotons_['diPhotonCand']

        if idx < 0:
            continue

        histStore['events']['nEvents'].Fill("nDiPhotonCands",1);
        histStore['events']['nEvents'].Fill("nDiPhotonCandsInWindow120to130",1);
        bJets=getBJetParis(eTree,diPhotons_['diPhotons'][idx]['index'])
        
        idxs=getBestGetMatchesGlobal(eTree)

        if(not bJets['valid']): 
            nNoHHto4B+=1
        
        quad=bJets['bJetQuad']
        histStore['events']['nEvents'].Fill("nBJetQuads",1);
        if( bJets['valid']): 
            if quad['r_HH'] < 25.0 :
                histStore['events']['nEvents'].Fill("nDiPhotonCandEventsIn25GevMbbSignalWindow",1);
   
        
        jetPts=np.array(eTree.jets_pt)
        mask= jetPts>25
        jetEta=np.array(eTree.jets_eta)
        mask=np.logical_and(mask , abs(jetEta) < 2.5 )
        nGoodJets=sum(mask)
        histStore['hhhTo4b2gamma']['eventJetMultiplicity'].Fill(nGoodJets)
        idxs=[]

        idxs=getBestGetMatchesGlobal(eTree)
        if(bJets['valid'] and quad['r_HH'] < 25.0):
            isValid,hasHH,hasLeadHiggs,hasSubLeadHiggs=checkGenMatches(eTree,quad['idxs'])
            hasRecoedH1=True
            hasRecoedH2=True
            hasRecoedHH=True
            if idxs[0] < 0 or idxs[1] <0 :
                hasRecoedH1=False
            if idxs[2] < 0 or idxs[3] <0 :
                hasRecoedH2=False
            hasRecoedHH = hasRecoedH1 and hasRecoedH2
            if hasRecoedH1:
                histStore["events"]["nHiggsMatch"].Fill("hasAllRecoedH1",1)      
            if hasRecoedH2:
                histStore["events"]["nHiggsMatch"].Fill("hasAllRecoedH2",1)      
            if hasRecoedHH:
                histStore["events"]["nHiggsMatch"].Fill("hasAllRecoedHH",1)      
            if isValid :
                histStore["events"]["nHiggsMatch"].Fill("AllInSR",1)      
            if isValid and hasHH:
                histStore["events"]["nHiggsMatch"].Fill("HH",1)      
            if isValid and hasLeadHiggs:
                histStore["events"]["nHiggsMatch"].Fill("H1",1)      
            if isValid and hasSubLeadHiggs:
                histStore["events"]["nHiggsMatch"].Fill("H2",1)      
        idx=idxs
        if idx[0] > 0:
            histStore["events"]["nHiggsMatch"].Fill("h1j1",1) 
            fillHHHJetHisto( histStore['jetIDStudy']['jets']['h1j1'] , eTree , idx[0] )
        if idx[1] > 0:
            histStore["events"]["nHiggsMatch"].Fill("h1j2",1)      
            fillHHHJetHisto( histStore['jetIDStudy']['jets']['h1j2'] , eTree , idx[1] )
        if idx[2] > 0:
            histStore["events"]["nHiggsMatch"].Fill("h2j1",1)      
            fillHHHJetHisto( histStore['jetIDStudy']['jets']['h2j1'] , eTree , idx[2] )
        if idx[3] > 0:
            histStore["events"]["nHiggsMatch"].Fill("h2j2",1)      
            fillHHHJetHisto( histStore['jetIDStudy']['jets']['h2j2'] , eTree , idx[3] )
        
        if bJets['valid']:
            if(quad['r_HH'] < 25.0):
                i1,i2,i3,i4=quad['idxs']
                fillHHHJetHisto( histStore['jetIDSelectedQuad']['jets']['h1j1'] , eTree , i1 )
                fillHHHJetHisto( histStore['jetIDSelectedQuad']['jets']['h1j2'] , eTree , i2 )
                fillHHHJetHisto( histStore['jetIDSelectedQuad']['jets']['h2j1'] , eTree , i3 )
                fillHHHJetHisto( histStore['jetIDSelectedQuad']['jets']['h2j2'] , eTree , i4 )


    simFile.Close()           
saveTheDictionary(histStore,foutName)
print(" File written out  : ",foutName)

