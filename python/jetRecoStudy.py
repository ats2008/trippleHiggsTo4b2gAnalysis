from __future__ import print_function
import ROOT 
import warnings
from collections import OrderedDict
import numpy as np
import branches as brList
from array import array
import Util  as utl
import trippleHiggsUtils as hhhUtil
import trippleHiggsSelector as hhhSelector
import trippleHiggsGenAnlayzer as hhhGen
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

headers=utl.getListOfStringsFromConfigs(cfgTxt,"#HEADER_BEG","#HEADER_END")
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
 
allFnames =utl.getListOfStringsFromConfigs(cfgTxt,"#FNAMES_BEG","#FNAMES_END")
foutName  =utl.getValueFromConfigs(cfgTxt,"OutpuFileName","fggHists.root")
processID =utl.getValueFromConfigs(cfgTxt,"processID",default="DATA")
treeName  =utl.getValueFromConfigs(cfgTxt,"treeName",default="tagsDumper/trees/Data_13TeV_TrippleHTag_0")
outTreeName=utl.getValueFromConfigs(cfgTxt,"outTreeName",default="Data_13TeV_TrippleHTag_0")
weightScale=float(utl.getValueFromConfigs(cfgTxt,"WeightScale",default="1.0"))
resetWeight=float(utl.getValueFromConfigs(cfgTxt,"resetWeight",default=-1e5))
doWeight=utl.getBoolFromConfigs(cfgTxt,"doWeight",default=False)
doPtReWeighting=utl.getBoolFromConfigs(cfgTxt,"doPtReWeighting",default=False)
pTReweitingFile=utl.getValueFromConfigs(cfgTxt,"pTReweitingFile",default="")
pTReweitingHistName=utl.getValueFromConfigs(cfgTxt,"pTReweitingHistName",default="")
mlScoreTag=utl.getValueFromConfigs(cfgTxt,"mlScoreTag",default="")
methord=utl.getValueFromConfigs(cfgTxt,"methord",default="DHH")
doOverlapRemoval =utl.getValueFromConfigs(cfgTxt,"doOverlapRemoval",default="1") ; doOverlapRemoval = int(doOverlapRemoval) > 0.5
isData =utl.getValueFromConfigs(cfgTxt,"isData",default="1") ; isData = int(isData) > 0.5
etaMax =float(utl.getValueFromConfigs(cfgTxt,"etaMax",default="2.5"))
pTMin =float(utl.getValueFromConfigs(cfgTxt,"pTMin",default="25.0"))
overlapRemovalDRMax =float(utl.getValueFromConfigs(cfgTxt,"overlapRemovalDRMax",default="0.4"))

print("allFnames   :  ",              allFnames)
print("foutName   :  ",               foutName)
print("processID   :  ",              processID)
print("treeName   :  ",               treeName)
print("weightScale   :  ",            weightScale)
print("doPtReWeighting   :  ",        doPtReWeighting)
print("pTReweitingFile   :  ",        pTReweitingFile)
print("pTReweitingHistName   :  ",    pTReweitingHistName)
print('etaMax   :       ', etaMax)
print('pTMin    :       ', pTMin)
print('drMax    :       ', drMax)
print('methord    :       ', methord)
print('mlScoreTag    :       ', mlScoreTag)

scaler=hhhUtil.mcScaler()
if doPtReWeighting:
    scaler.setSFHistFromFile(pTReweitingFile,pTReweitingHistName)
print("")


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
filemode="RECREATE"
foutName=foutName.replace('.root',ext+'.root')
fout = ROOT.TFile(foutName, filemode)
dir_ = fout.mkdir('trees')

branches=[
     'CMS_hgg_mass', 'D_HH', 'H1H2JetAbsCosThetaMax', 'H1H2JetAbsCosThetaMin',
     'H1H2JetDrMax', 'H1H2JetDrMin', 'H1bbToH2bbAbsCosTheta',
     'H1bbToH2bbAbsCosTheta', 'HggTo4bAbsCosTheta',
     'LeadJetAbsCosThetaMax', 'LeadJetAbsCosThetaMin',
     'LeadJetAbsCosThetaMin','LeadJetDrMaxWithOtherJets',
     'LeadJetDrMaxWithOtherJets', 'LeadJetDrMinWithOtherJets', 'PhoJetMinDr','PhoJetMinDrOther','PhoJetMaxDr','PhoJetMaxDrOther',
     'absCosThetaH4bHgg', 'customLeadingPhotonIDMVA','pThgg_overMgg',
     'customSubLeadingPhotonIDMVA','h1bb_mass', 'diphotonPtOverDiphotonMass',
     'H1H2JetDrMax','h1bbCosThetaLeadJet','H2bbCosThetaLeadJet','PhoJetMinDr',
     'h1bb_mass', 'h1bb_phi', 'h1bb_pt','h2bb_eta' , 'h2bb_mass', 'h1bb_eta', 'h2bb_phi','h2bb_pt', 
     'HH4bCosThetaLeadJet','leadingPhotonSigOverE', 'pT_4b', 'pTh1leadJ_overMh1',
     'pTh1subleadJ_overMh1', 'pTh2leadJ_overMh2', 'pTh2subleadJ_overMh2',
     'pTh2subleadJ_overMh2','absCosThetaH4bHgg', 'pTleadG_overMgg',
     'pTsubleadG_overMgg', 'quadjet_0_deepJetScore', 'quadjet_1_deepJetScore',
     'quadjet_2_deepJetScore', 'quadjet_3_deepJetScore', 'r_HH', 'scalarPtSum4b',
     'scalarPtSum4b2g', 'scalarPtSum4b2g','HggTo4bAbsCosTheta', 'scalarPtSumHHH',
     'h1_dijetSigmaMOverM', 'h2_dijetSigmaMOverM', 'sigmaMOverM',
     'subleadingPhotonSigOverE', 'sumScore_3j', 'sumScore_4j', 'trihiggs_mass',
     'ttH_MET', 'weight','rho','trihiggs_pt',
     "CosThetaH1_hhhF",  "HH4bCosTheta_hhhF",  "HggCosTheta_hhhF",  "HH4bCosThetaLeadJet_hhhF",  "absCosThetaH4bHgg_hhhF",  
     'H1bbCosTheta_hhhF','H2bbCosTheta_hhhF'
]
branches+=['quadjet_0_mlScore', 'quadjet_0_mlScoreY1s', 'quadjet_0_mlScoreY2s' ,
    'quadjet_1_mlScore', 'quadjet_1_mlScoreY1s', 'quadjet_1_mlScoreY2s' ,
    'quadjet_2_mlScore', 'quadjet_2_mlScoreY1s', 'quadjet_2_mlScoreY2s' ,
    'quadjet_3_mlScore', 'quadjet_3_mlScoreY1s', 'quadjet_3_mlScoreY2s' ]

ntuple={}
branches=np.unique(branches)
#print(branches)

ntuple['eventTree']  = ROOT.TNtuple(outTreeName, outTreeName, ':'.join(branches))

tofill = OrderedDict(zip(branches, [np.nan]*len(branches)))
outputDataDict=OrderedDict(zip(branches, [np.nan]*len(branches))) 

kyList = [ ky for ky in tofill ]

print("len(branches) : " , len(branches))

sumEntries=ROOT.TH1F("sumEvts","sumEvts",1,0.0,1.0)
sumEntries.SetCanExtend(ROOT.TH1.kAllAxes)
sumWeights=ROOT.TH1F("sumWeighs","sumWeighs",1,0.0,1.0)
sumWeights.SetCanExtend(ROOT.TH1.kAllAxes)
isDataHist=ROOT.TH1F("isData","isData",1,0.0,1.0)
isDataHist.SetCanExtend(ROOT.TH1.kAllAxes)

th1Store={}


beg=datetime.datetime.now()
has_printed=True
for fname in allFnames:
    print("Opening file : ",fname)
    simFile = ROOT.TFile(fname,'READ')
    eTree=simFile.Get(treeName)
    if not eTree:
        treeName='tagsDumper/trees/EGamma_13TeV_TrippleHTag_0'
        print("Trying tree name : ",treeName)
        eTree=simFile.Get(treeName)
    if not eTree:
        treeName='tagsDumper/trees/ggHHH_125_13TeV'
        print("Trying tree name : ",treeName)
        eTree=simFile.Get(treeName)
    print(" NEntries = ", eTree.GetEntries())
    maxEvents_ = eTree.GetEntries()
    if(maxEvents >0  and (totalEvents+maxEvents_) > maxEvents):
        maxEvents_= (maxEvents - totalEvents)
    totalEvents+=maxEvents_
    
    allBranches=[]
    for ky in eTree.GetListOfBranches():
        allBranches.append(ky.GetName())
    isData_= isData or 'EGamma' in fname.split('/')[-1] 
    isData_= isData_ or ('DoubleEG' in fname.split('/')[-1])
    print(" isData : ",isData_)

    for evt in range(maxEvents_):
        eTree.GetEntry(evt)
        
        wei=eTree.weight
        if doWeight:
            if resetWeight > -1e4:
                wei=resetWeight
            if weightScale > -1e4:
                wei*=weightScale
            if doPtReWeighting:
                pT=eTree.diphoton_pt
                scaleFactor=scaler.getSFForX(pT)
                wei*=scaleFactor

        sumEntries.Fill('total', 1)
        sumWeights.Fill('total', wei)
        i=evt
        if(i%250==0):
            now=datetime.datetime.now()
            timeSpendSec=np.round((now-beg).total_seconds() , 2)
            timeLeftSec =np.round(1.0*(maxEvents_-i)*timeSpendSec/( i +1e-3),2)
            print("   Doing i = ",i," / ",maxEvents_)
            print("      time left : ", str(datetime.timedelta(seconds= timeLeftSec)),
                    " [ time elapsed : ",datetime.timedelta(seconds= timeSpendSec), " s ]")
            print(" gg mass   : ",eTree.CMS_hgg_mass)
         
        _tmp=False
        for ky in tofill:
            if ky in allBranches:
                tofill[ky]=getattr(eTree,ky)
            elif not has_printed:
                warnings.warn(ky+" : variable  not avalable in allBranches ! ",category=UserWarning,stacklevel=0)
                tofill[ky]=0.0
                _tmp=True
        if _tmp:
            has_printed=True
        
        ##   JET PRE-SELECTION
        print()
        print("Processing event : ",i)
        hhhGen.printGenInfo(eTree)

        text=""
        jetMask=[]
        jetMask=hhhSelector.getSelectedJetCollectionMaskEta(eTree,jetMask=[],etaMax=etaMax)
        if sum(jetMask) < 4 :
            sumEntries.Fill('nJetPreselectionEta',1)
            sumWeights.Fill('nJetPreselectionEta',wei)
            text+='Eta Sel. : 0 '
        jetMask=hhhSelector.getSelectedJetCollectionMaskPt(eTree,jetMask=jetMask,pTMin=pTMin)
        if sum(jetMask) < 4 :
            sumEntries.Fill('nJetPreselectionPt',1)
            sumWeights.Fill('nJetPreselectionPt',wei)
            text+='Pt Sel. : 0 '
        if doOverlapRemoval:
            jetMask=hhhSelector.getSelectedJetCollectionMaskOverLap(eTree,jetMask=jetMask,overlapRemovalDRMax=overlapRemovalDRMax)
            if sum(jetMask) < 4 :
                sumEntries.Fill('nJetPreselectionOR',1)
                sumWeights.Fill('nJetPreselectionOR',wei)
                text+='OLR Sel. : 0 '
        #if len(text)<1:
        #    text+="Preselected"
        
        rslt=hhhSelector.getBestGetMatchesGlobalFGG(eTree,jetMask=jetMask)
        idxs=rslt['fgg_idxs']
        isRecoed={}
        jetMatchesString=""
        for i in range(4):
            isRecoed['isRecoed_'+str(i)]=True
        for i in range(4):
            if idxs[i] < 0:
                isAllRecoed=False
                isRecoed['isRecoed_'+str(i)]=False
                jetMatchesString+=" b"+str(i+1)+"("+u'\u2717'+")" 
            else:
                jetMatchesString+=" b"+str(i+1)+"(j"+str(idxs[i])+")"
        h1Out= isRecoed['isRecoed_0'] and isRecoed['isRecoed_1']
        h2Out=  isRecoed['isRecoed_2'] and isRecoed['isRecoed_3']
        
        nout=0
        if not isRecoed['isRecoed_0']:            nout+=1;    
        if not isRecoed['isRecoed_1']:            nout+=1;    
        if not isRecoed['isRecoed_2']:            nout+=1;    
        if not isRecoed['isRecoed_3']:            nout+=1;

        if h1Out:
            jetMatchesString+="| H1("+u'\u2713'+")"
        else:
            jetMatchesString+="| H1("+u'\u2717'+")"
        if h2Out:
            jetMatchesString+=" H2("+u'\u2713'+")"
        else:
            jetMatchesString+=" H2("+u'\u2717'+")"
        
        x=0
        for i in range(8):
            if i in idxs:
                x+=1
        cat=' |'    
        if nout==0 and x==4 :
            cat+='allRecoed |'
        elif nout==1 and x==3 :
            cat+=' SingleJetMiss |'
        if nout<=1 and x>=3 :
            cat+=' Max1JetMiss |'
        if nout>1  :
            cat+=' MultiJetMiss |'
        idx_=[k for k in idxs]
        while -1 in idx_:
            idx_.remove(-1)
        if x!=len(np.unique(idx_)) :
            cat+=' MergedJetMiss |'

        text='Event '+str(evt)+" "+text +" , "+jetMatchesString+cat
        #if rslt['idxs'][1] >1:
        #    continue
        #if rslt['idxs'][1] <0:
        #    continue
        #print("b2 issue  idx ",rslt['idxs'],"  rslt : ",rslt['fgg_idxs']," | btag : ",rslt['bJetScoreRanks'])
        #
        fpicname='results/plots/genPlots/evtDepiction_'+str(evt)+'.pdf'
        print(fpicname)
        hhhUtil.visualizeEvents(eTree,text=text,outputFname=fpicname,jetMask=jetMask)
        

    simFile.Close()           
    print("Closing file : ",fname)
dir_.cd()    

sumEntries.Write()
sumWeights.Write()
isDataHist.Write()
for ky in th1Store:
    th1Store[ky].Write()
for ky in ntuple:
    ntuple[ky].Write()

fout.Close()

print(" File written out  : ",foutName)
