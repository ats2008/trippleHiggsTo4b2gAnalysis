from __future__ import print_function
import ROOT 
import numpy as np
from trippleHiggsGenAnlayzer import *
import trippleHiggsUtils as hhhUtil

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
treeName=getValueFromConfigs(cfgTxt,"treeName",default="ntuplizer/tree")

maxEvents=-1
tmp_=getValueFromConfigs(cfgTxt,"MaxEvents")
if tmp_!='':
    maxEvents=int(tmp_)

for i in allFnames:
    print(" file : ",i)
if(maxEvtsSuperSeeder > 0):
    maxEvents=maxEvtsSuperSeeder
print("maxevents : ",maxEvents)
print("treeName : ",treeName)

histStore = getGenHistos()

print(histStore.keys())
print(histStore['kinematicVars'].keys())
hhhUtil.addObjectDeltaRValues(histStore)

nNoHtoGammaGamma=0
nNoHHto4B=0
totalEvents=0

isMC = True

bLV=[]
for i in range(4):
    bLV.append(ROOT.TLorentzVector())
gammaLV=[]
for i in range(4):
    gammaLV.append(ROOT.TLorentzVector())

for fname in allFnames:
    
    simFile = ROOT.TFile(fname,'READ')
    eTree=simFile.Get(treeName)
    if not eTree:
        eTree=simFile.Get(treeName)
        
    maxEvents_ = eTree.GetEntries()
    if(maxEvents >0  and (totalEvents+maxEvents_) > maxEvents):
        maxEvents_= (maxEvents - totalEvents)
    totalEvents+=maxEvents_

    for i in range(maxEvents_):
        eTree.GetEntry(i)
        bIdx=0
        gammaIdx=0        
        
        if(i%100==0):
            print("Doing i = ",i," / ",maxEvents_)
        

        bDaus=hhhUtil.getHiggsDauP4s(eTree,5)       
        gammaDaus=hhhUtil.getHiggsDauP4s(eTree,22)       
        bIdx=len(bDaus)
        gammaIdx=len(gammaDaus)
        
        bLV=bDaus
        gammaLV=gammaDaus
        if(gammaIdx!=2 or bIdx!=4):
            print("Problem !!! ",gammaIdx,bIdx)
        
        nHiggsAbove300=0
        
        if ( bLV[0] +bLV[1] ).Pt() > 300: nHiggsAbove300+=1
        if ( bLV[2] +bLV[3] ).Pt() > 300: nHiggsAbove300+=1
        if ( gammaLV[0] +gammaLV[1] ).Pt() > 300: nHiggsAbove300+=1

        histStore['vars']['nHiggsAbove300'].Fill(str(nHiggsAbove300),1)

        histStore['vars']['mass_4b'].Fill( (bLV[0]+bLV[1]+bLV[2]+bLV[3]).M())
        histStore['vars']['mass_2gamma'].Fill( (gammaLV[0]+gammaLV[1]).M())
        histStore['vars']['mass_2bb'].Fill( (bLV[0]+bLV[1]).M())
        histStore['vars']['mass_2bb'].Fill( (bLV[2]+bLV[3]).M())
        histStore['vars']['mass_4b2gamma'].Fill( (bLV[0]+bLV[1]+bLV[2]+bLV[3] + gammaLV[0]+ gammaLV[1]).M())
        
        allDauP4=bLV[0]+bLV[1]+bLV[2]+bLV[3] + gammaLV[0]+ gammaLV[1]
        diPhotonLV=gammaLV[0]+ gammaLV[1]

        quad={}
        if (bLV[0]+bLV[1]).Pt() > (bLV[2]+bLV[3]).Pt() :
            quad['p4_h1']=bLV[0]+bLV[1]
            quad['p4_h2']=bLV[2]+bLV[3]
        else:
            quad['p4_h2']=bLV[0]+bLV[1]
            quad['p4_h1']=bLV[2]+bLV[3]
            
        mx0 =  allDauP4.M() - ( quad['p4_h1'].M()+quad['p4_h2'].M()+diPhotonLV.M() - 375.0 )
        mx1 = (quad['p4_h1']+quad['p4_h2']).M() - ( quad['p4_h1'].M()+quad['p4_h2'].M() - 250.0 ) 
        mx2 = (quad['p4_h1']+diPhotonLV).M()-( quad['p4_h1'].M()+diPhotonLV.M() - 250.0 ) 
        mx3 = (quad['p4_h2']+diPhotonLV).M()-( quad['p4_h2'].M()+diPhotonLV.M() - 250.0 ) 
        


        histStore['vars']['mass_X0'].Fill( mx0)
        histStore['vars']['mass_X1'].Fill( mx1)
        histStore['vars']['mass_X2'].Fill( mx2)
        histStore['vars']['mass_X3'].Fill( mx3)
    
    
        histStore['vars']['H1Pt'].Fill(eTree.gen_H1_pt)
        histStore['vars']['H2Pt'].Fill(eTree.gen_H2_pt)
        histStore['vars']['H3Pt'].Fill(eTree.gen_H3_pt)
        
        histStore['vars']['H1Eta'].Fill(eTree.gen_H1_eta)
        histStore['vars']['H2Eta'].Fill(eTree.gen_H2_eta)
        histStore['vars']['H3Eta'].Fill(eTree.gen_H3_eta)
        
        histStore['vars']['H1y'].Fill(eTree.gen_H1_y)
        histStore['vars']['H2y'].Fill(eTree.gen_H2_y)
        histStore['vars']['H3y'].Fill(eTree.gen_H3_y)
        
        histStore['vars']['H1Phi'].Fill(eTree.gen_H1_phi)
        histStore['vars']['H2Phi'].Fill(eTree.gen_H1_phi)
        histStore['vars']['H3Phi'].Fill(eTree.gen_H1_phi)
        
        histStore['vars']['H1H2DeltaR'].Fill(deltaR(eTree.gen_H1_eta,eTree.gen_H1_phi,eTree.gen_H2_eta,eTree.gen_H2_phi))
        histStore['vars']['H2H3DeltaR'].Fill(deltaR(eTree.gen_H2_eta,eTree.gen_H2_phi,eTree.gen_H3_eta,eTree.gen_H3_phi))
        histStore['vars']['H3H1DeltaR'].Fill(deltaR(eTree.gen_H3_eta,eTree.gen_H3_phi,eTree.gen_H1_eta,eTree.gen_H1_phi))
    
    
        histStore['vars']['H1H2DeltaEta'].Fill(abs(eTree.gen_H1_eta - eTree.gen_H2_eta))
        histStore['vars']['H2H3DeltaEta'].Fill(abs(eTree.gen_H2_eta - eTree.gen_H3_eta))
        histStore['vars']['H3H1DeltaEta'].Fill(abs(eTree.gen_H3_eta - eTree.gen_H1_eta))
    
    
        histStore['vars']['H1H2DeltaPhi'].Fill(abs(eTree.gen_H1_phi - eTree.gen_H2_phi))
        histStore['vars']['H2H3DeltaPhi'].Fill(abs(eTree.gen_H2_phi - eTree.gen_H3_phi))
        histStore['vars']['H3H1DeltaPhi'].Fill(abs(eTree.gen_H3_phi - eTree.gen_H1_phi))
        for ii in range(4):
            #print(ii," ) ",np.round(bDaus[ii].Pt(),2 ), " | ",end="")
            histStore['vars']["B"+str(ii+1)+"Pt"].Fill(bDaus[ii].Pt())
            histStore['vars']["B"+str(ii+1)+"Eta"].Fill(bDaus[ii].Eta())
            histStore['vars']["B"+str(ii+1)+"Phi"].Fill(bDaus[ii].Phi())
            histStore['vars']["B"+str(ii+1)+"Y"].Fill(bDaus[ii].Rapidity())

        for ii in range(2):
            histStore['vars']["gamma"+str(ii+1)+"Pt"].Fill(gammaDaus[ii].Pt())
            histStore['vars']["gamma"+str(ii+1)+"Eta"].Fill(gammaDaus[ii].Eta())
            histStore['vars']["gamma"+str(ii+1)+"Phi"].Fill(gammaDaus[ii].Phi())
        
        LVStore=getLVStoreFromGen(eTree)
        fillTrippleHGenVariablesFromLVStore(histStore,LVStore)
        hhhUtil.fillObjectDeltaRValuesFromLVStore(histStore, LVStore)

    simFile.Close()           

saveTheDictionary(histStore,foutName)
print(" File written out  : ",foutName)

