from __future__ import print_function
import ROOT 
import warnings
from collections import OrderedDict
import numpy as np
import trippleHiggsUtils as hhhUtil
import Util  as utl
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
#branches+=['quadjet_0_mlScore', 'quadjet_0_mlScoreY1s', 'quadjet_0_mlScoreY2s' ,
#    'quadjet_1_mlScore', 'quadjet_1_mlScoreY1s', 'quadjet_1_mlScoreY2s' ,
#    'quadjet_2_mlScore', 'quadjet_2_mlScoreY1s', 'quadjet_2_mlScoreY2s' ,
#    'quadjet_3_mlScore', 'quadjet_3_mlScoreY1s', 'quadjet_3_mlScoreY2s' ]

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
    print(" NEntries = ", eTree.GetEntries())
    if not eTree:
        eTree=simFile.Get('tagsDumper/trees/EGamma_13TeV_TrippleHTag_0')
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

    for i in range(maxEvents_):
        eTree.GetEntry(i)
        
        isMasked = isData_ and  (eTree.CMS_hgg_mass > 115.0) and (eTree.CMS_hgg_mass < 135.0)
        if isMasked:
            continue
        
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
        jetMask=hhhSelector.getSelectedJetCollectionMaskEta(eTree,etaMax=etaMax)
        if sum(jetMask) < 4 :
            sumEntries.Fill('nJetPreselectionEta',1)
            sumWeights.Fill('nJetPreselectionEta',wei)
            continue

        jetMask=hhhSelector.getSelectedJetCollectionMaskPt(eTree,jetMask=jetMask,pTMin=pTMin)
        if sum(jetMask) < 4 :
            sumEntries.Fill('nJetPreselectionPt',1)
            sumWeights.Fill('nJetPreselectionPt',wei)
            continue
        if doOverlapRemoval:
            jetMask=hhhSelector.getSelectedJetCollectionMaskOverLap(eTree,jetMask=jetMask,overlapRemovalDRMax=overlapRemovalDRMax)
            if sum(jetMask) < 4 :
                sumEntries.Fill('nJetPreselectionOR',1)
                sumWeights.Fill('nJetPreselectionOR',wei)
                continue
    
 
        verbose=False
        sumScore_3j=0.0
        sumScore_4j=0.0

        if methord=='DHH':
            allQuads=hhhSelector.getBJetParisFGG(eTree,mask=jetMask)
            quad=allQuads['bJetQuad']
            i=0
            for idx in quad['fgg_idxs']:
                val=max(getattr(eTree,'jet_'+str( idx )+'_deepCSVScore'),0.0)
                tofill['quadjet_'+str(i)+'_deepJetScore']   = val
                if i< 3:
                    sumScore_3j+=val
                sumScore_4j+=val
                i+=1

        elif methord=='mha':
            allQuads=hhhSelector.getBJetParis_wrapper(eTree, mask= jetMask , methord='mha',mlScoreTag='H3SIN61',threshold=-1e3,doMisclassificationCorrection = False)
            if not allQuads['isValid']:
                print("anity check failed after the quad finding")
                exit(1)
                continue
            if allQuads['jetsValid']:
                i=0
                for idx in allQuads['allJetsSelected']:
                    tofill['quadjet_'+str(i)+'_deepJetScore']   =getattr(eTree,'jet_'+str( idx )+'_deepCSVScore') 
                    tofill['quadjet_'+str(i)+'_mlScore']        =getattr(eTree,'jet_'+str( idx )+mlScoreTag+'_score')
                    tofill['quadjet_'+str(i)+'_mlScoreY1s']     =getattr(eTree,'jet_'+str( idx )+mlScoreTag+'_y0s')
                    tofill['quadjet_'+str(i)+'_mlScoreY2s']     =getattr(eTree,'jet_'+str( idx )+mlScoreTag+'_y1s')
                    if tofill['quadjet_'+str(i)+'_mlScore']   < -1:
                        verbose=True
                    if i < 3 :
                        sumScore_3j+=getattr(eTree,'jet_'+str( idx )+mlScoreTag+'_score')
                    if i < 4 :
                        sumScore_4j+=getattr(eTree,'jet_'+str( idx )+mlScoreTag+'_score')
                        
                    i+=1

        if verbose:
            sumEntries.Fill('badEvents' , 1)
            sumWeights.Fill('badEvents' , wei)
            continue
            #print()
            #print("Event  : ",eTree.event)
            #print(jetMask)
            #for i in range(8):
            #    print("Jet : ",i)
            #    print("\t pT ",getattr(eTree,'jet_'+str(i)+'_pt'))
            #    print("\t eta ",getattr(eTree,'jet_'+str(i)+'_eta'))
            #    print("\t phi ",getattr(eTree,'jet_'+str(i)+'_phi'))
            #    print("\t m ",getattr(eTree,'jet_'+str(i)+'_mass'))
            #    print("\t deepFlav ",getattr(eTree,'jet_'+str( i )+'_deepCSVScore') )
            #    print("\t score ",getattr(eTree,'jet_'+str(i)+mlScoreTag+'_score'))
            #i=0
            #for idx in allQuads['allJetsSelected']:
            #    print(tofill['quadjet_'+str(i)+'_deepJetScore'])
            #    print(tofill['quadjet_'+str(i)+'_mlScore']     )
            #    print(tofill['quadjet_'+str(i)+'_mlScoreY1s']  )
            #    print(tofill['quadjet_'+str(i)+'_mlScoreY2s']  )
            #    i+=1
            
        tofill["sumScore_4j"]=sumScore_4j
        tofill["sumScore_3j"]=sumScore_3j

        quad=allQuads['bJetQuad']
        LVStore = hhhUtil.getLVStoreFromTreeAndQuad(eTree,quad)
        j1CosTheta,k1CosTheta,ggCostheta,drMin,drOther=hhhUtil.getCosthetaVars(eTree,LVStore)
        
        tofill['r_HH'] = quad['r_HH']
        tofill['D_HH'] = quad['D_HH']
        tofill['weight'] =  wei

        varDict=hhhUtil.getOtherDerivedVariables(eTree,LVStore,quad)
        #for ky in tofill:
        #    if ky not in varDict:
        #       print("filling default barnch value  for branch",ky)
        for ky in varDict:
            tofill[ky]=varDict[ky]
        tofill['ttH_MET'] = eTree.ttH_MET


        for ky in outputDataDict:
            if ky not in tofill:
                print("not in tofill ",ky)

        for i in outputDataDict:
            outputDataDict[i]=tofill[i]
        arr=array('f', outputDataDict.values())
        ntuple['eventTree'].Fill(arr)

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
