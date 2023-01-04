from __future__ import print_function
from collections import OrderedDict
import ROOT 
import numpy as np
from array import array
from TMVA_Model import *
import datetime

import Util  as utl
import  trippleHiggsMLInterface as hhhMLI
import trippleHiggsUtils as hhhUtil
import trippleHiggsSelector as hhhSelector

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

allFnames= utl.getListOfStringsFromConfigs(cfgTxt,"#FNAMES_BEG","#FNAMES_END")
foutName= utl.getValueFromConfigs(cfgTxt,"OutpuFileName","fggHists.root")
processID= utl.getValueFromConfigs(cfgTxt,"processID",default= "DATA")
treeName= utl.getValueFromConfigs(cfgTxt,"treeName",default= "tagsDumper/trees/Data_13TeV_TrippleHTag_0")
outTreeName= utl.getValueFromConfigs(cfgTxt,"outTreeName",default="Data_13TeV_TrippleHTag_0")
eventIdxOffset= int(utl.getValueFromConfigs(cfgTxt,"EventIdxOffset",default= 0))
mlScoreTag= utl.getValueFromConfigs(cfgTxt,"mlScoreTag",default= "")
doOverlapRemoval = utl.getValueFromConfigs(cfgTxt,"doOverlapRemoval",default= "1") ; doOverlapRemoval =  int(doOverlapRemoval) > 0.5
isData =utl.getValueFromConfigs(cfgTxt,"isData",default="1") ; isData = int(isData) > 0.5
etaMax =float(utl.getValueFromConfigs(cfgTxt,"etaMax",default="2.5"))
pTMin =float(utl.getValueFromConfigs(cfgTxt,"pTMin",default="25.0"))
overlapRemovalDRMax =float(utl.getValueFromConfigs(cfgTxt,"overlapRemovalDRMax",default="0.4"))
methord=utl.getValueFromConfigs(cfgTxt,"methord",default="DHH")


MVAWeightFile= utl.getValueFromConfigs(cfgTxt,"MVAWeightFile",default= "")
MVATag= utl.getValueFromConfigs(cfgTxt,"mvaScoreTag","")
MVABranches= utl.getListOfStringsFromConfigs(cfgTxt,"#MVAVARLIST_BEG","#MVAVARLIST_END")
MVASpecBranches= utl.getListOfStringsFromConfigs(cfgTxt,"#SPECTATORLIST_BEG","#SPECTATORLIST_END")

print("allFnames   :  ",              allFnames)
print("foutName   :  ",               foutName)
print("processID   :  ",              processID)
print("treeName   :  ",               treeName)
print("MVAWeightFile    :  ", MVAWeightFile)
print("MVATag    :  ", MVATag)

nonResonantMVA=TMVAModel()
nonResonantMVA.setupTMVAModel("MVA_"+MVATag, MVAWeightFile , MVABranches , MVASpecBranches )


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
nMatchFound=0
fname=allFnames[0]
simFile = ROOT.TFile(fname,'READ')
eTree=simFile.Get(treeName)
allBranchesFromFile=[]
for ky in eTree.GetListOfBranches():
    allBranchesFromFile.append(ky.GetName())

branchesToFill=allBranchesFromFile

branches=[ 'pTh1jj_overMh1','pTh2jj_overMh2','pThgg_overMgg','pTleadG_overMgg','pTsubleadG_overMgg','pTh1leadJ_overMh1',
'pTh1subleadJ_overMh1','pTh2leadJ_overMh2','pTh2subleadJ_overMh2',"h1bbCosThetaLeadJet",
'h1bb_pt','h2bb_pt','h1bb_eta','h2bb_eta','h1bb_phi','h2bb_phi','h1bb_mass','h2bb_mass',
"HH4bCosThetaLeadJet","absCosThetaH4bHgg","CosThetaH1_hhhF","HH4bCosTheta_hhhF","HggCosTheta_hhhF",
"H1bbCosTheta_hhhF","H2bbCosTheta_hhhF","HH4bCosThetaLeadJet_hhhF","absCosThetaH4bHgg_hhhF",
"LeadJetAbsCosThetaMax","LeadJetAbsCosThetaMin","LeadJetDrMaxWithOtherJets","LeadJetDrMinWithOtherJets",
"H1H2JetAbsCosThetaMax","H1H2JetAbsCosThetaMin","H1H2JetDrMax","H1H2JetDrMin","PhoJetMinDr","PhoJetMinDrOther",
"PhoJetMaxDr","PhoJetMaxDrOther",'dije1CandidatePtOverHiggsM','dije2CandidatePtOverHiggsM',"H1bbCosThetaLeadJet",
'H2bbCosThetaLeadJet',"pT_4b","scalarPtSumHHH","scalarPtSum4b","scalarPtSum4b2g","HggTo4bAbsCosTheta",
"H1bbToH2bbAbsCosTheta",'h1leadJ_deepjetScore','h1subleadJ_deepjetScore','h2leadJ_deepjetScore','h2subleadJ_deepjetScore',
'h1_dijetSigmaMOverM','h2_dijetSigmaMOverM','trihiggs_mass','trihiggs_pt','ttH_MET']

for br in branches:
    branchesToFill.append(br)

branches+=['isBadEvent']
branches+=['MVA_'+MVATag]


filemode="RECREATE"
fout = ROOT.TFile(foutName, filemode)
dir_ = fout.mkdir('tagsDumper')
dir_ = dir_.mkdir('trees')

branches=np.unique(branchesToFill)
tofill = OrderedDict(zip(branches, [np.nan]*len(branches)))
outputDataDict=OrderedDict(zip(branches, [np.nan]*len(branches))) 

ntuple={}
ntuple['all'] = ROOT.TNtuple(outTreeName, outTreeName, ':'.join(branches))


kyList = [ ky for ky in tofill ]
print("len(branches) : " , len(branches))


m1m2 = ROOT.TH2F("m1jj_m2jj","H1bb , H2bb mass",300,0.0,300.0,300,0.0,300. )


sumEntries=ROOT.TH1F("sumEntries","sumWeights",1,0.0,1.0)
sumEntries.SetCanExtend(ROOT.TH1.kAllAxes)
sumWeights=ROOT.TH1F("sumWeights","sumWeights",1,0.0,1.0)
sumWeights.SetCanExtend(ROOT.TH1.kAllAxes)
th1Store={}

LVStore={}
LVStore['jetLV']=ROOT.TLorentzVector();
LVStore['g1LV']=ROOT.TLorentzVector()
LVStore['g2LV']=ROOT.TLorentzVector()
LVStore['HggLV'] =ROOT.TLorentzVector()

beg=datetime.datetime.now()
eventIdx=eventIdxOffset
nProcessed=0
for fname in allFnames:
    print("Opening file : ",fname)
    simFile = ROOT.TFile(fname,'READ')
    eTree=simFile.Get(treeName)
    print(" NEntries = ", eTree.GetEntries() )
    if not eTree:
        eTree=simFile.Get('tagsDumper/trees/EGamma_13TeV_TrippleHTag_0')
    maxEvents_ = eTree.GetEntries()
    if(maxEvents >0  and (totalEvents+maxEvents_) > maxEvents):
        maxEvents_= (maxEvents - totalEvents)
    totalEvents+=maxEvents_
    allBranches=[]
    for ky in eTree.GetListOfBranches():
        allBranches.append(ky.GetName())
    for evt in range(maxEvents_):
        eTree.GetEntry(evt)
        wei=eTree.weight
        sumWeights.Fill('total', wei)

        if(evt%250==0):
            now=datetime.datetime.now()
            i=evt
            timeSpendSec=np.round((now-beg).total_seconds() , 2)
            timeLeftSec =np.round(1.0*(maxEvents_-i)*timeSpendSec/( i +1e-3),2)
            print("   Doing i = ",i," / ",maxEvents_, "  Processed : ",nProcessed )
            print("      time left : ", str(datetime.timedelta(seconds= timeLeftSec)),
                    " [ time elapsed : ",datetime.timedelta(seconds= timeSpendSec), " s ]")
            print(" gg mass   : ",eTree.CMS_hgg_mass)
         
        for ky in tofill:
            if ky in allBranches:
                tofill[ky]=getattr(eTree,ky)
            else:
                tofill[ky]=0.0
        
        
        ##   JET PRE-SELECTION
        isBadEvents=False
        jetMask=hhhSelector.getSelectedJetCollectionMaskEta(eTree,etaMax=etaMax)
        if sum(jetMask) < 4 :
            sumEntries.Fill('nJetPreselectionEta',1)
            sumWeights.Fill('nJetPreselectionEta',wei)
            isBadEvents=True

        jetMask=hhhSelector.getSelectedJetCollectionMaskPt(eTree,jetMask=jetMask,pTMin=pTMin)
        if sum(jetMask) < 4 :
            sumEntries.Fill('nJetPreselectionPt',1)
            sumWeights.Fill('nJetPreselectionPt',wei)
            isBadEvents=True
        if doOverlapRemoval:
            jetMask=hhhSelector.getSelectedJetCollectionMaskOverLap(eTree,jetMask=jetMask,overlapRemovalDRMax=overlapRemovalDRMax)
            if sum(jetMask) < 4 :
                sumEntries.Fill('nJetPreselectionOR',1)
                sumWeights.Fill('nJetPreselectionOR',wei)
                isBadEvents=True

 
        verbose=False
        sumScore_3j=0.0
        sumScore_4j=0.0

        if methord=='DHH':
            allQuads=hhhSelector.getBJetParisFGG(eTree,mask=jetMask)
            quad=allQuads['bJetQuad']
            if allQuads['isValid']:
                i=0
                for idx in quad['fgg_idxs']:
                    val=max(getattr(eTree,'jet_'+str( idx )+'_deepCSVScore'),0.0)
                    tofill['quadjet_'+str(i)+'_deepJetScore']   = val
                    if i< 3:
                        sumScore_3j+=val
                    sumScore_4j+=val
                    i+=1
            else:
                isBadEvents=True

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
            isBadEvents=True
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
        if not isBadEvents:
            quad=allQuads['bJetQuad']
            LVStore = hhhUtil.getLVStoreFromTreeAndQuad(eTree,quad)
            varDict=hhhUtil.getOtherDerivedVariables(eTree,LVStore,quad)
            for ky in varDict:
                tofill[ky]=varDict[ky]
            tofill['r_HH'] = quad['r_HH']
            tofill['D_HH'] = quad['D_HH']
        tofill['weight'] =  wei

        if not isBadEvents:
            v1=0.5*( 1.0 + nonResonantMVA.predict(tofill) )
        else:
            v1=-0.5
        tofill["MVA_"+MVATag]=v1
        #print("Mva score : ",v1)
        tofill['isBadEvent']=isBadEvents
        for i in outputDataDict:
            outputDataDict[i]=tofill[i]
        
        ntuple['all'].Fill(array('f', outputDataDict.values()))
        if not isBadEvents : nProcessed+=1

    simFile.Close()           
    print("Closing file : ",fname)
print("Number of good events ",nProcessed)

dir_.cd()    
sumWeights.Write()
m1m2.Write()

for cat in ntuple:
    ntuple[cat].Write()
fout.Close()
print(" File written out  : ",foutName)

