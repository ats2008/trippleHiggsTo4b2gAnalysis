from __future__ import print_function
from collections import OrderedDict
import ROOT 
import numpy as np
from trippleHiggsUtils import *
from Util  import *
from branches import *
from array import array
from trippleHiggsSelector import *

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
    scaler.setSFHistFromFile(pTReweitingFile,pTReweitingHistName)
    scalerVal.setSFHistFromFile(pTReweitingValFile,pTReweitingHistValName)
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
filemode="RECREATE"
fout = ROOT.TFile(foutName, filemode)
if filemode=='update':
    ntuple = fout.Get('tree')
else:
    ntuple = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
tofill = OrderedDict(zip(branches, [np.nan]*len(branches)))


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
    allBranches=[]
    for ky in eTree.GetListOfBranches():
        allBranches.append(ky.GetName())
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
            pT=eTree.diphoton_pt
            scaleFactor=scaler.getSFForX(pT)
            wei*=scaleFactor
            scaleFactor=scalerVal.getSFForX(pT)
            weiVal*=scaleFactor
        
        isBTaggedAllLoose = (eTree.h1LeadingJet_DeepFlavour > 0.0490 ) and ( eTree.h1SubleadingJet_DeepFlavour > 0.0490 ) and ( eTree.h2LeadingJet_DeepFlavour > 0.0490 ) and ( eTree.h2SubleadingJet_DeepFlavour > 0.0490 )
        #if not isBTaggedAllLoose:
        #    continue

        ## Control and signal region
        dr = np.sqrt((eTree.M1jj-125.0)**2 + (eTree.M2jj - 125.0)**2)
        
        isMasked =  (dr < 25.0) and maskSignalMgg and processID=="DATA" and ( eTree.CMS_hgg_mass > 115.0) and (eTree.CMS_hgg_mass < 135.0)
        
        if isMasked:
            continue

        for ky in tofill:
            if ky in allBranches:
                tofill[ky]=getattr(eTree,ky)
     #           print(ky,"--> ",tofill[ky])
            else:
     #           print("\t Setting ",ky," -1.111e3")
                tofill[ky]=-1.111e3
        LVStore = getLVStore(eTree)
        j1CosTheta,k1CosTheta,ggCostheta,drMin,drOther=getCosthetaVars(eTree,LVStore)
        
        tofill['weight'] =  wei
        tofill['hhhCosThetaH1'] = eTree.absCosThetaStar_CS 
        tofill["hh4CosThetaLeadJet"] = eTree.absCosTheta_bb
        tofill["h1bbCosThetaLeadJet"]= abs(j1CosTheta)
        tofill["h2bbCosThetaLeadJet"]= abs(k1CosTheta)
        tofill["h2bbCosThetaLeadJet"]= abs(k1CosTheta)
        tofill["PhoJetMinDr"]= drMin
        tofill["PhoJetOtherDr"]= drOther
        tofill['pTleadG_overMgg'] =  eTree.leadingPhoton_pt / eTree.CMS_hgg_mass 
        tofill['pTh1leadJ_overMh1'] =  eTree.h1LeadingJet_pt / eTree.M1jj 
        tofill['pTh2leadJ_overMh2'] =  eTree.h2LeadingJet_pt / eTree.M2jj 

        tofill['pTsubleadG_overMgg'] =  eTree.subleadingPhoton_pt / eTree.CMS_hgg_mass 
        tofill['pTh1subleadJ_overMh1'] =  eTree.h1SubleadingJet_pt / eTree.M1jj 
        tofill['pTh2subleadJ_overMh2'] =  eTree.h2SubleadingJet_pt / eTree.M2jj 

#        if( dr < 25.0 ):
#            if isMasked:
#                continue
#            fillCandidateHistograms(histStore["trippleHCands"],eTree,wei)     
#        elif ( dr < 50.0 ):
#            fillCandidateHistograms(histStore["controlCands"],eTree,wei)     
#        
#        dr = np.sqrt((eTree.M1jj-200.0)**2 + (eTree.M2jj - 200.0)**2)
#        if( dr < 25.0 ):
#            fillCandidateHistograms(histStore["validation_SR"],eTree,weiVal)     
#        elif ( dr < 50.0 ):
#            fillCandidateHistograms(histStore["validation_CR"],eTree,weiVal)     
        
        ntuple.Fill(array('f', tofill.values()))
        
    simFile.Close()           
    print("Closing file : ",fname)
fout.cd()    
ntuple.Write()
fout.Close()
print(" File written out  : ",foutName)

