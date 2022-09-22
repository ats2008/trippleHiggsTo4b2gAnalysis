from __future__ import print_function
from collections import OrderedDict
import ROOT
from array import array
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
outTreeName=getValueFromConfigs(cfgTxt,"outTreeName",default="Data_13TeV_TrippleHTag_0")
maskSignalMgg=getBoolFromConfigs(cfgTxt,"maskSignalMgg",default=True)
doBjetCounting=getBoolFromConfigs(cfgTxt,"doBjetCounting",default=False)
minNBjetsFromMC=int(getValueFromConfigs(cfgTxt,"minNBjetsFromMC",default="0"))
maxNBjetsFromMC=int(getValueFromConfigs(cfgTxt,"maxNBjetsFromMC",default="1000"))
weightScale=float(getValueFromConfigs(cfgTxt,"WeightScale",default="1.0"))
resetWeight=float(getValueFromConfigs(cfgTxt,"resetWeight",default=-1e5))

doPtReWeighting=getBoolFromConfigs(cfgTxt,"doPtReWeighting",default=False)
pTReweitingFile=getValueFromConfigs(cfgTxt,"pTReweitingFile",default="")
pTReweitingHistName=getValueFromConfigs(cfgTxt,"pTReweitingHistName",default="")
pTReweitingHistCatagories=getValueFromConfigs(cfgTxt,"pTReweitingHistCatagories",default="")

pTReweitingHistName=pTReweitingHistName.split(',')
pTReweitingHistCatagories=pTReweitingHistCatagories.split(',')

reweighter={}
for i,j in zip(pTReweitingHistCatagories,pTReweitingHistName):
    print(i,j)
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
#print("pTReweitingValFile   :  ",     pTReweitingValFile)
#print("pTReweitingHistValName   :  ", pTReweitingHistValName)


scaler=mcScaler()
scalerVal=mcScaler()
if doPtReWeighting:
    print(reweighter)
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

#addCandidateVars(histStore,['allCands'])# ,'trippleHCands','controlCands','validation_CR','validation_SR'])
#histStore["hhhTo4b2gamma"]={}
#histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_tightID"]= histStore["allCands"]["h1MassVsh2Mass"].Clone()
#histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_tightID"].SetName("h1MassVsh2Mass_tightID")
#histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_mediumID"]= histStore["allCands"]["h1MassVsh2Mass"].Clone()
#histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_mediumID"].SetName("h1MassVsh2Mass_mediumID")
#histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_looseID"]= histStore["allCands"]["h1MassVsh2Mass"].Clone()
#histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_looseID"].SetName("h1MassVsh2Mass_looseID")
#histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_looseID"]= histStore["allCands"]["h1MassVsh2Mass"].Clone()
#histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_looseID"].SetName("h1MassVsh2Mass_looseID")
#histStore["h1Tobb"]={}
#histStore["h2Tobb"]={}
#histStore["h1Tobb"]["mass_tightID" ] =ROOT.TH1F("mass_tightID","",2000,0.0,2000)
#histStore["h1Tobb"]["mass_mediumID" ]=ROOT.TH1F("mass_mediumID","",2000,0.0,2000)
#histStore["h1Tobb"]["mass_looseID" ] =ROOT.TH1F("mass_looseID","",2000,0.0,2000)
#histStore["h2Tobb"]["mass_tightID" ] =ROOT.TH1F("mass_tightID","",2000,0.0,2000)
#histStore["h2Tobb"]["mass_mediumID" ]=ROOT.TH1F("mass_mediumID","",2000,0.0,2000)
#histStore["h2Tobb"]["mass_looseID" ] =ROOT.TH1F("mass_looseID","",2000,0.0,2000)

nNoHtoGammaGamma=0
nNoHHto4B=0
totalEvents=0



########################################     Ntuple defenition


filemode="RECREATE"
fout = ROOT.TFile(foutName, filemode)
dir_ = fout.mkdir('tagsDumper')
dir_ = dir_.mkdir('trees')
dir_.cd()
branches_skimmed=[
                  'M1jj','M2jj','CMS_hgg_mass','dZ','weight','triHiggs_mass'
                  'h1LeadingJet_pt','h1SubleadingJet_pt','h2LeadingJet_pt','h2SubleadingJet_pt',
                  'h1LeadingJet_eta','h1SubleadingJet_eta','h2LeadingJet_eta','h2SubleadingJet_eta',
                  'h1LeadingJet_puJetIdMVA','h1SubleadingJet_puJetIdMVA','h2LeadingJet_puJetIdMVA','h2SubleadingJet_puJetIdMVA',
                  'h1LeadingJet_DeepFlavour','h1SubleadingJet_DeepFlavour','h2LeadingJet_DeepFlavour','h2SubleadingJet_DeepFlavour',
                 ]

ntuple={}
ntuple['all'] = ROOT.TNtuple(outTreeName+'_NOTAG', outTreeName, ':'.join(branches_skimmed))
tofill = OrderedDict(zip(branches_skimmed, [np.nan]*len(branches_skimmed)))

##                           <------------------------------------->

def cutBasedSelectionFGG( 
                         eTree,
                         cutDetails={
                                     'minPt'  : 25.0 ,
                                     'etaMax' : 2.5 ,
                                     'puID'   : 'medium',
                                     'm1jj'   : 70.0,
                                     'm2jj'   : 70.0,
                                     'mHHH_min' : 300.0
                                     },
                         
                       ):


    if eTree.h1LeadingJet_pt    <  cutDetails['minPt'] : return False
    if eTree.h1SubleadingJet_pt <  cutDetails['minPt'] : return False
    if eTree.h2LeadingJet_pt    <  cutDetails['minPt'] : return False
    if eTree.h2SubleadingJet_pt <  cutDetails['minPt'] : return False

    if abs(eTree.h1LeadingJet_eta)    > cutDetails['etaMax']    : return False
    if abs(eTree.h1SubleadingJet_eta) > cutDetails['etaMax']    : return False
    if abs(eTree.h2LeadingJet_eta)    > cutDetails['etaMax']    : return False
    if abs(eTree.h2SubleadingJet_eta) > cutDetails['etaMax']    : return False

    if not puJetID( eTree.h1LeadingJet_pt    , eTree.h1LeadingJet_puJetIdMVA , cutDetails['puID'] )  : return False
    if not puJetID( eTree.h1SubleadingJet_pt , eTree.h1SubleadingJet_pt      , cutDetails['puID'] )  : return False
    if not puJetID( eTree.h2LeadingJet_pt    , eTree.h2LeadingJet_pt         , cutDetails['puID'] )  : return False
    if not puJetID( eTree.h2SubleadingJet_pt , eTree.h2SubleadingJet_pt      , cutDetails['puID'] )  : return False
    
    if eTree.M1jj < cutDetails['m1jj'] : return False
    if eTree.M2jj < cutDetails['m2jj'] : return False
    if eTree.triHiggs_mass < cutDetails['mHHH_min']     : return False
    
    return True



##    < ------------------- Category defenition ----------------------------- >

tagList=['']
catList=['cat'+str(i) for i in range(6) ]
catList.append('X')
if processID=='DATA':
    tagList.append('_fitRange')
for tag in tagList:
    if filemode=='update':
        for cat in catList:
            ntuple[cat+tag] = fout.Get(outTreeName+'_'+cat)
    else:
        for cat in catList:
            ntuple[cat+tag] = ROOT.TNtuple(outTreeName+tag+'_'+cat, outTreeName+tag+'_'+cat, ':'.join(branches_skimmed))

def getCategory(eTree):
    
    if btagID(eTree.h1LeadingJet_DeepFlavour,'tight') and btagID(eTree.h1SubleadingJet_DeepFlavour,'tight') \
        and btagID(eTree.h2LeadingJet_DeepFlavour,'tight') and btagID(eTree.h2SubleadingJet_DeepFlavour,'tight'):
        
        return 'cat0'
    
    if btagID(eTree.h1LeadingJet_DeepFlavour,'tight') and btagID(eTree.h1SubleadingJet_DeepFlavour,'tight') \
        and btagID(eTree.h2LeadingJet_DeepFlavour,'medium') and btagID(eTree.h2SubleadingJet_DeepFlavour,'medium'):
        
        return 'cat1'
    
    if btagID(eTree.h1LeadingJet_DeepFlavour,'tight') and btagID(eTree.h1SubleadingJet_DeepFlavour,'tight') \
        and btagID(eTree.h2LeadingJet_DeepFlavour,'loose') and btagID(eTree.h2SubleadingJet_DeepFlavour,'loose'):
        
        return 'cat2'
    
    if btagID(eTree.h1LeadingJet_DeepFlavour,'medium') and btagID(eTree.h1SubleadingJet_DeepFlavour,'medium') \
        and btagID(eTree.h2LeadingJet_DeepFlavour,'loose') and btagID(eTree.h2SubleadingJet_DeepFlavour,'loose'):
        
        return 'cat3'
    
    if btagID(eTree.h1LeadingJet_DeepFlavour,'loose') and btagID(eTree.h1SubleadingJet_DeepFlavour,'loose') \
        and btagID(eTree.h2LeadingJet_DeepFlavour,'loose') and btagID(eTree.h2SubleadingJet_DeepFlavour,'loose'):
        
        return 'cat4'
    
    if btagID(eTree.h1LeadingJet_DeepFlavour,'loose') and btagID(eTree.h1SubleadingJet_DeepFlavour,'loose') :
        return 'cat5'
    
    return 'X'


########################################     <<-------------------


isMC = True
isMC = False

sumEvts=ROOT.TH1F("sumEvts","sumEvts",1,0.0,1.0)
sumEvts.SetCanExtend(ROOT.TH1.kAllAxes)
sumWeights=ROOT.TH1F("sumWeights","sumWeights",1,0.0,1.0)
sumWeights.SetCanExtend(ROOT.TH1.kAllAxes)

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
        if resetWeight > -1e4:
            wei=resetWeight
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
        ## ----->  Event Selection <--------
        
        if not cutBasedSelectionFGG(eTree):
            continue
        ## ----->                  <--------


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

        weiVal=wei
        
        isBTaggedAllLoose = (eTree.h1LeadingJet_DeepFlavour > 0.0490 ) and ( eTree.h1SubleadingJet_DeepFlavour > 0.0490 ) and ( eTree.h2LeadingJet_DeepFlavour > 0.0490 ) and ( eTree.h2SubleadingJet_DeepFlavour > 0.0490 )
        
        dr = np.sqrt((eTree.M1jj-125.0)**2 + (eTree.M2jj - 125.0)**2)
        isMasked =  maskSignalMgg and processID=="DATA" and ( eTree.CMS_hgg_mass > 115.0) and (eTree.CMS_hgg_mass < 135.0)

        for ky in tofill:
            if ky in allBranches:
                tofill[ky]=getattr(eTree,ky)
            else:
                tofill[ky]=-1.111e3

        cat=getCategory(eTree)

        sumEvts.Fill('total_'+cat, 1)
        sumWeights.Fill('total_wei'+cat, wei)
        tofill['weight']=wei
        ntuple[cat].Fill(array('f', tofill.values()))
        if tag!='':
            ntuple[cat+tag].Fill(array('f', tofill.values()))
        ntuple['all'].Fill(array('f', tofill.values()))

        
    simFile.Close()           
    print("Closing file : ",fname)
dir_.cd()
for cat in ntuple:
    ntuple[cat].Write()
#saveTheDictionary(histStore,None,fout)

sumEvts.Write()
sumWeights.Write()
fout.Close()

print(" File written out  : ",foutName)

