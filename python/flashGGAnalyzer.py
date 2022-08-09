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
foutName=getValueFromConfigs(cfgTxt,"OutpuFileName","fggHists.root")
processID=getValueFromConfigs(cfgTxt,"processID",default="DATA")
treeName=getValueFromConfigs(cfgTxt,"treeName",default="tagsDumper/trees/Data_13TeV_TrippleHTag_0")
maskSignalMgg=getBoolFromConfigs(cfgTxt,"maskSignalMgg",default=True)
doBjetCounting=getBoolFromConfigs(cfgTxt,"doBjetCounting",default=False)
minNBjetsFromMC=int(getValueFromConfigs(cfgTxt,"minNBjetsFromMC",default="0"))
maxNBjetsFromMC=int(getValueFromConfigs(cfgTxt,"maxNBjetsFromMC",default="1000"))
weightScale=float(getValueFromConfigs(cfgTxt,"WeightScale",default="1.0"))
doPtReWeighting=getBoolFromConfigs(cfgTxt,"doPtReWeighting",default=False)
pTReweitingFile=getValueFromConfigs(cfgTxt,"pTReweitingFile",default="")
pTReweitingHistName=getValueFromConfigs(cfgTxt,"pTReweitingHistName",default="")
pTReweitingValFile=getValueFromConfigs(cfgTxt,"pTReweitingValFile",default="")
pTReweitingHistValName=getValueFromConfigs(cfgTxt,"pTReweitingHistValName",default="")
pTReweitingHistCatagories=getValueFromConfigs(cfgTxt,"pTReweitingHistCatagories",default="")
pTReweitingHistValName=pTReweitingHistValName.split(',')
pTReweitingHistCatagories=pTReweitingHistCatagories.split(',')
reweighter={}
for i,j in zip(pTReweitingHistCatagories,pTReweitingHistValName):
    reweighter[i]=j
    

print("allFnames   :  ",              allFnames)
print("foutName   :  ",               foutName)
print("processID   :  ",              processID)
print("treeName   :  ",               treeName)
print("maskSignalMgg   :  ",          maskSignalMgg)
print("doBjetCounting   :  ",         doBjetCounting)
print("minNBjetsFromMC   :  ",        minNBjetsFromMC)
print("maxNBjetsFromMC   :  ",        maxNBjetsFromMC)
print("weightScale   :  ",            weightScale)
print("doPtReWeighting   :  ",        doPtReWeighting)
print("pTReweitingFile   :  ",        pTReweitingFile)
print("pTReweitingHistName   :  ",    pTReweitingHistName)
print("pTReweitingValFile   :  ",     pTReweitingValFile)
print("pTReweitingHistValName   :  ", pTReweitingHistValName)


scaler=mcScaler()
scalerVal=mcScaler()
if doPtReWeighting:
    scaler.setSFHistFromFile(pTReweitingFile,reweighter)
#    scalerVal.setSFHistFromFile(pTReweitingValFile,pTReweitingHistValName)

print("")

maxEvents=-1
tmp_=getValueFromConfigs(cfgTxt,"MaxEvents")
if tmp_!='':
    maxEvents=int(tmp_)

for i in allFnames:
    print(" file : ",i)
if(maxEvtsSuperSeeder > 0):
    maxEvents=maxEvtsSuperSeeder
print("maxevents : ",maxEvents)

histStore = {}

addCandidateVars(histStore,['allCands' ,'trippleHCands','controlCands','validation_CR','validation_SR'])

nNoHtoGammaGamma=0
nNoHHto4B=0
totalEvents=0

isMC = True
isMC = False

for fname in allFnames:
    print("Opening file : ",fname)
    simFile = ROOT.TFile(fname,'READ')
    eTree=simFile.Get(treeName)
    print(" NEntries = ", eTree.GetEntries())
    if not eTree:
        eTree=simFile.Get('tagsDumper/trees/Data_13TeV_TrippleHTag_0')
    maxEvents_ = eTree.GetEntries()
    if(maxEvents >0  and (totalEvents+maxEvents_) > maxEvents):
        maxEvents_= (maxEvents - totalEvents)
    totalEvents+=maxEvents_
    
    for i in range(maxEvents_):
        eTree.GetEntry(i)
        if(i%500==0):
            print("   Doing i = ",i," / ",maxEvents_,
                  " n No HiggsCand : ", nNoHtoGammaGamma,
                  " n No 4BCands  : ", nNoHHto4B)
            print(" gg mass , h1MassVsh2Mass  : ",eTree.CMS_hgg_mass ,"( ",eTree.M1jj , eTree.M2jj,")" )
            

        wei=eTree.weight
        if processID=="MC":
            wei*=weightScale
            if doBjetCounting:
                nBs=getNBsFromCand(eTree)
                if nBs > maxNBjetsFromMC:
                    continue
                if nBs < minNBjetsFromMC:
                    continue
        elif processID=="DATA":
            pass

        weiVal=wei
        if doPtReWeighting:
            cat='BB'
            if abs(eTree.leadingPhoton_eta) < 1.44 and abs(eTree.subleadingPhoton_eta) > 1.567:
                cat='BE'
            if abs(eTree.leadingPhoton_eta) > 1.567 and abs(eTree.subleadingPhoton_eta) < 1.44:
                cat='BE'
            if abs(eTree.leadingPhoton_eta) > 1.567 and abs(eTree.subleadingPhoton_eta) < 1.567:
                cat='EE'
            pT=eTree.diphoton_pt
            scaleFactor=scaler.getSFForX(pT,cat)
            wei*=scaleFactor
            #print("Setting SF : ",scaleFactor," | for pT : ",pT," , cat : ",cat)
            #scaleFactor=scalerVal.getSFForX(pT,cat)
            #weiVal*=scaleFactor
        
        isBTaggedAllLoose = (eTree.h1LeadingJet_DeepFlavour > 0.0490 ) and ( eTree.h1SubleadingJet_DeepFlavour > 0.0490 ) and ( eTree.h2LeadingJet_DeepFlavour > 0.0490 ) and ( eTree.h2SubleadingJet_DeepFlavour > 0.0490 )
        

        #if not isBTaggedAllLoose:
        #    continue

        ## Control and signal region
        dr = np.sqrt((eTree.M1jj-125.0)**2 + (eTree.M2jj - 125.0)**2)
        
        isMasked =  (dr < 25.0) and maskSignalMgg and processID=="DATA" and ( eTree.CMS_hgg_mass > 115.0) and (eTree.CMS_hgg_mass < 135.0)
        
        if not isMasked:
            fillCandidateHistograms(histStore["allCands"],eTree,wei)

        if( dr < 25.0 ):
            if isMasked:
                continue
            fillCandidateHistograms(histStore["trippleHCands"],eTree,wei)     
        elif ( dr < 50.0 ):
            fillCandidateHistograms(histStore["controlCands"],eTree,wei)     
        
        dr = np.sqrt((eTree.M1jj-200.0)**2 + (eTree.M2jj - 200.0)**2)
        if( dr < 25.0 ):
            fillCandidateHistograms(histStore["validation_SR"],eTree,weiVal)     
        elif ( dr < 50.0 ):
            fillCandidateHistograms(histStore["validation_CR"],eTree,weiVal)     
        
        continue

        if(eTree.h1LeadingJet_DeepFlavour > 0.7100 and eTree.h1SubleadingJet_DeepFlavour > 0.7100 
            and  eTree.h2LeadingJet_DeepFlavour > 0.7100 and eTree.h2SubleadingJet_DeepFlavour > 0.7100 ):
            histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_tightID"].Fill(eTree.M1jj , eTree.M2jj )
            histStore["h1Tobb"]["mass_tightID" ].Fill(eTree.M1jj)
            histStore["h2Tobb"]["mass_tightID" ].Fill(eTree.M2jj)

        if(eTree.h1LeadingJet_DeepFlavour > 0.2783 and eTree.h1SubleadingJet_DeepFlavour > 0.2783 
            and  eTree.h2LeadingJet_DeepFlavour > 0.2783 and eTree.h2SubleadingJet_DeepFlavour > 0.2783 ):
            histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_mediumID"].Fill(eTree.M1jj , eTree.M2jj )
            histStore["h1Tobb"]["mass_mediumID" ].Fill(eTree.M1jj)
            histStore["h2Tobb"]["mass_mediumID" ].Fill(eTree.M2jj)

        
        if(eTree.h1LeadingJet_DeepFlavour > 0.0490 and eTree.h1SubleadingJet_DeepFlavour > 0.0490 
            and  eTree.h2LeadingJet_DeepFlavour > 0.0490 and eTree.h2SubleadingJet_DeepFlavour > 0.0490 ):
            histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_looseID"].Fill(eTree.M1jj , eTree.M2jj )
            histStore["h1Tobb"]["mass_looseID" ].Fill(eTree.M1jj)
            histStore["h2Tobb"]["mass_looseID" ].Fill(eTree.M2jj)

        
    simFile.Close()           
    print("Closing file : ",fname)
saveTheDictionary(histStore,foutName)
print(" File written out  : ",foutName)

