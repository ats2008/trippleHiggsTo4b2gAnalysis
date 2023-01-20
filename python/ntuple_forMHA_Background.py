from __future__ import print_function
import ROOT 
from collections import OrderedDict
import numpy as np
import trippleHiggsUtils as hhhUtil
from Util  import *
import branches as brList
from array import array
import trippleHiggsSelector as hhhSelector
from TMVA_Model import *

import datetime

import os,sys

cfgFileName=''
if len(sys.argv) <2:
    print("Usage\n\t ~$python recoAnalyzer.py <configFile>\n")
    exit(1)
else:
    cfgFileName=sys.argv[1]
maxEvtsSuperSeeder=-1
etaMax=2.5
pTMin=25.0
drMax=0.4
ext=''
if len(sys.argv) >2:
    maxEvtsSuperSeeder=int(sys.argv[2])
if len(sys.argv) >3:
    etaMax=float(sys.argv[3])
if len(sys.argv) >4:
    pTMin =float(sys.argv[4])
if len(sys.argv) >5:
    drMax =float(sys.argv[5])
if len(sys.argv) >6:
    ext=sys.argv[6]

print("Loading cfg file ",cfgFileName)
f=open(cfgFileName,'r')
cfgTxt=f.readlines()
f.close()

headers=getListOfStringsFromConfigs(cfgTxt,"#HEADER_BEG","#HEADER_END")
for header in headers:
    print("Loading cfg file ",header)
    if not os.path.exists(header):
        print("\t file  : ",header," do not exists")
        continue
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
maskSignalMgg=getBoolFromConfigs(cfgTxt,"maskSignalMgg",default=True)
doBjetCounting=getBoolFromConfigs(cfgTxt,"doBjetCounting",default=False)
minNBjetsFromMC=int(getValueFromConfigs(cfgTxt,"minNBjetsFromMC",default="0"))
maxNBjetsFromMC=int(getValueFromConfigs(cfgTxt,"maxNBjetsFromMC",default="1000"))
weightScale=float(getValueFromConfigs(cfgTxt,"WeightScale",default="1.0"))
resetWeight=float(getValueFromConfigs(cfgTxt,"resetWeight",default=-1e5))
doSR=getBoolFromConfigs(cfgTxt,"doSR",default=False)
doOverlapRemoval =getValueFromConfigs(cfgTxt,"doOverlapRemoval",default="1") ; doOverlapRemoval = int(doOverlapRemoval) > 0.5
isData =getValueFromConfigs(cfgTxt,"isData",default="1") ; isData = int(isData) > 0.5
etaMax =float(getValueFromConfigs(cfgTxt,"etaMax",default="2.5"))
pTMin =float(getValueFromConfigs(cfgTxt,"pTMin",default="25.0"))
overlapRemovalDRMax =float(getValueFromConfigs(cfgTxt,"overlapRemovalDRMax",default="0.4"))

doPtReWeighting=getBoolFromConfigs(cfgTxt,"doPtReWeighting",default=False)
pTReweitingFile=getValueFromConfigs(cfgTxt,"pTReweitingFile",default="")
pTReweitingHistName=getValueFromConfigs(cfgTxt,"pTReweitingHistName",default="")
pTReweitingHistCatagories=getValueFromConfigs(cfgTxt,"pTReweitingHistCatagories",default="")

pTReweitingHistName=pTReweitingHistName.split(',')
pTReweitingHistCatagories=pTReweitingHistCatagories.split(',')

treesToStore=getValueFromConfigs(cfgTxt,"TreesToStore",default="background")
treesToStore=treesToStore.split(',')

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
print("doSR     :       ", doSR)
print('etaMax   :       ', etaMax)
print('pTMin    :       ', pTMin)
print('drMax    :       ', drMax)

scaler=hhhUtil.mcScaler()
if doPtReWeighting:
    print(reweighter)
    scaler.setSFHistFromFile(pTReweitingFile,reweighter)

print("")


skipBad= 'allEvents' not in treesToStore


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
filemode="RECREATE"
foutName=foutName.replace('.root',ext+'.root')
fout = ROOT.TFile(foutName, filemode)
dir_ = fout.mkdir('trees')

branches=brList.allMCBranchesV2
branches.append('matchedJet_b1')
branches.append('matchedJet_b2')
branches.append('matchedJet_b3')
branches.append('matchedJet_b4')
for i in range(8):
    branches.append('label_'+str(i))
for i in range(8):
    branches.append('label_idx_'+str(i))

for i in range(8):
    branches.append('jet_'+str(i)+'_drWithDP')
    branches.append('jet_'+str(i)+'_drWithDP_leadG')
    branches.append('jet_'+str(i)+'_dr_WithDPSubleadG')
    branches.append('jet_'+str(i)+'_dEta_WithDP')
    branches.append('jet_'+str(i)+'_dEta_WithDPLeadG')
    branches.append('jet_'+str(i)+'_dEta_WithDPSubleadG')
    branches.append('jet_'+str(i)+'_dPhi_WithDP')
    branches.append('jet_'+str(i)+'_dPhi_WithDPLeadG')
    branches.append('jet_'+str(i)+'_dPhi_WithDPSubleadG')
    branches.append('jet_'+str(i)+'_mass_WithDP')
    branches.append('jet_'+str(i)+'_mass_WithDPLeadG')
    branches.append('jet_'+str(i)+'_mass_WithDPSubleadG')
    branches.append('jet_'+str(i)+'_massRatio_WithDP')
    branches.append('jet_'+str(i)+'_ptRatio_WithDP')
    branches.append('jet_'+str(i)+'_ptRatio_WithDPLeadG')
    branches.append('jet_'+str(i)+'_ptRatio_WithDPSubleadG')
 

ntuple={}
branches=np.unique(branches)
ntuple['allEvents']  = ROOT.TNtuple(outTreeName+'_allEvents', outTreeName, ':'.join(branches))
ntuple['background'] = ROOT.TNtuple(outTreeName+'_background', outTreeName, ':'.join(branches))

tofill = OrderedDict(zip(branches, [np.nan]*len(branches)))

kyList = [ ky for ky in tofill ]

print("len(branches) : " , len(branches))

m1m2 = ROOT.TH2F("m1jj_m2jj","H1bb , H2bb mass",300,0.0,300.0,300,0.0,300. )

sumEntries=ROOT.TH1F("sumEvts","sumEvts",1,0.0,1.0)
sumEntries.SetCanExtend(ROOT.TH1.kAllAxes)
sumWeights=ROOT.TH1F("sumWeighs","sumWeighs",1,0.0,1.0)
sumWeights.SetCanExtend(ROOT.TH1.kAllAxes)

th1Store={}

LVStore={}
LVStore['jetLV']=ROOT.TLorentzVector();
LVStore['g1LV']=ROOT.TLorentzVector()
LVStore['g2LV']=ROOT.TLorentzVector()
LVStore['HggLV'] =ROOT.TLorentzVector()

beg=datetime.datetime.now()

for fname in allFnames:
    print("Opening file : ",fname)
    simFile = ROOT.TFile(fname,'READ')
    eTree=simFile.Get(treeName)
    print(" NEntries = ", eTree.GetEntries())
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
    for i in range(maxEvents_):
        
        eTree.GetEntry(i)
        wei=eTree.weight
        sumEntries.Fill('total', 1)
        sumWeights.Fill('total', wei)

        if(i%250==0):
            now=datetime.datetime.now()
            timeSpendSec=np.round((now-beg).total_seconds() , 2)
            timeLeftSec =np.round(1.0*(maxEvents_-i)*timeSpendSec/( i +1e-3),2)
            print("   Doing i = ",i," / ",maxEvents_)
            print("      time left : ", str(datetime.timedelta(seconds= timeLeftSec)),
                    " [ time elapsed : ",datetime.timedelta(seconds= timeSpendSec), " s ]")
            print(" gg mass   : ",eTree.CMS_hgg_mass)
            
        for ky in tofill:
            if ky in allBranches:
                tofill[ky]=getattr(eTree,ky)
            else:
                tofill[ky]=0.0
        #continue
        tofill['matchedJet_b1']=-1
        tofill['matchedJet_b2']=-1
        tofill['matchedJet_b3']=-1
        tofill['matchedJet_b4']=-1
        
        ##   JET PRE-SELECTION
        jetMask=hhhSelector.getSelectedJetCollectionMaskEta(eTree,jetMask=[],etaMax=etaMax)
        if sum(jetMask) < 3 :
            sumEntries.Fill('nJetPreselectionEta',1)
            sumWeights.Fill('nJetPreselectionEta',wei)
            if skipBad: continue

        jetMask=hhhSelector.getSelectedJetCollectionMaskPt(eTree,jetMask=jetMask,pTMin=pTMin)
        if sum(jetMask) < 3 :
            sumEntries.Fill('nJetPreselectionPt',1)
            sumWeights.Fill('nJetPreselectionPt',wei)
            if skipBad: continue
        if doOverlapRemoval:
            jetMask=hhhSelector.getSelectedJetCollectionMaskOverLap(eTree,jetMask=jetMask,overlapRemovalDRMax=overlapRemovalDRMax)
            if sum(jetMask) < 3 :
                sumEntries.Fill('nJetPreselectionOR',1)
                sumWeights.Fill('nJetPreselectionOR',wei)
                if skipBad: continue
    
        nVld=0
        for i in range(8):
            tofill['label_'+str(i)]=0
            tofill['label_idx_'+str(i)]=0
            if not jetMask[i]:
                tofill['jet_'+str(i)+'_isValid']=0
            if tofill['jet_'+str(i)+'_isValid'] > 0.5:
                nVld+=1

        LVStore['g1LV'].SetPtEtaPhiM(eTree.leadingPhoton_pt,eTree.leadingPhoton_eta,eTree.leadingPhoton_phi,0.0)
        LVStore['g2LV'].SetPtEtaPhiM(eTree.subleadingPhoton_pt,eTree.subleadingPhoton_eta,eTree.subleadingPhoton_phi,0.0)
        LVStore['HggLV'].SetPtEtaPhiM( eTree.diphoton_pt , eTree.diphoton_eta , eTree.diphoton_phi , eTree.CMS_hgg_mass )

        for i in range(8):
            if  getattr(eTree,'jet_'+str(i)+'_isValid') >0.3:
                LVStore['jetLV'].SetPtEtaPhiM(  
                                                getattr(eTree,'jet_'+str(i)+'_pt'),
                                                getattr(eTree,'jet_'+str(i)+'_eta'),
                                                getattr(eTree,'jet_'+str(i)+'_phi'),
                                                getattr(eTree,'jet_'+str(i)+'_mass')
                                             )
                tofill['jet_'+str(i)+'_drWithDP'] = LVStore['jetLV'].DeltaR( LVStore['HggLV'] )
                tofill['jet_'+str(i)+'_drWithDP_leadG'] = LVStore['jetLV'].DeltaR( LVStore['g1LV'] )
                tofill['jet_'+str(i)+'_dr_WithDPSubleadG'] = LVStore['jetLV'].DeltaR( LVStore['g2LV'] )
                tofill['jet_'+str(i)+'_dEta_WithDP'] = abs( LVStore['jetLV'].Eta() - LVStore['HggLV'].Eta() ) 
                tofill['jet_'+str(i)+'_dEta_WithDPLeadG'] = abs( LVStore['jetLV'].Eta() - LVStore['g1LV'].Eta() )
                tofill['jet_'+str(i)+'_dEta_WithDPSubleadG'] = abs( LVStore['jetLV'].Eta() - LVStore['g2LV'].Eta() )
                
                tofill['jet_'+str(i)+'_dPhi_WithDP'] = LVStore['jetLV'].DeltaPhi( LVStore['HggLV'] )
                tofill['jet_'+str(i)+'_dPhi_WithDPLeadG'] = LVStore['jetLV'].DeltaPhi( LVStore['g1LV'] )
                tofill['jet_'+str(i)+'_dPhi_WithDPSubleadG'] = LVStore['jetLV'].DeltaPhi( LVStore['g2LV'] )
                
                tofill['jet_'+str(i)+'_mass_WithDP'] = ( LVStore['jetLV'] + LVStore['HggLV'] ).M()
                tofill['jet_'+str(i)+'_mass_WithDPLeadG'] = ( LVStore['jetLV'] + LVStore['g1LV'] ).M()
                tofill['jet_'+str(i)+'_mass_WithDPSubleadG'] = ( LVStore['jetLV'] + LVStore['g2LV'] ).M()
                
                tofill['jet_'+str(i)+'_massRatio_WithDP'] = ( LVStore['jetLV'] + LVStore['HggLV'] ).M() /  LVStore['HggLV'].M()
                tofill['jet_'+str(i)+'_ptRatio_WithDP'] = ( LVStore['jetLV'] + LVStore['HggLV'] ).Pt() / LVStore['HggLV'].Pt()
                tofill['jet_'+str(i)+'_ptRatio_WithDPLeadG'] = ( LVStore['jetLV'] + LVStore['g1LV'] ).Pt() / LVStore['g1LV'].Pt()
                tofill['jet_'+str(i)+'_ptRatio_WithDPSubleadG'] = ( LVStore['jetLV'] + LVStore['g2LV'] ).Pt() / LVStore['g2LV'].Pt()

        for ky in branches:
            if ky not in tofill:
                print(ky)
        for ky in tofill:
            if ky not in branches:
                print(ky)
        if 'allEvents' in treesToStore:
            arr=array('f', tofill.values())
            ntuple['allEvents'].Fill(arr)
        elif  'background' in treesToStore:
            ntuple['background'].Fill(array('f', tofill.values()))

    simFile.Close()           
    print("Closing file : ",fname)
dir_.cd()    

sumEntries.Write()
sumWeights.Write()

for ky in th1Store:
    th1Store[ky].Write()
for ky in ntuple:
    ntuple[ky].Write()

fout.Close()
print(" File written out  : ",foutName)

