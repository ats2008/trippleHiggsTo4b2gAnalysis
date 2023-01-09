from __future__ import print_function
from collections import OrderedDict
import ROOT 
import numpy as np
from trippleHiggsUtils import *
from Util  import *
from branches import *
from array import array
import trippleHiggsSelector as hhhSelector
from TMVA_Model import *

import  trippleHiggsMLInterface as hhhMLI

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

headers=getListOfStringsFromConfigs(cfgTxt,"#HEADER_BEG","#HEADER_END")
for header in headers:
    print("Loading cfg file ",header)
    f=open(header,'r')
    tmp=f.readlines()
    f.close()
    for l in tmp:
        cfgTxt.append(l)

def insideTheEllipse( x,y,x1,y1,x2,y2,a):
    return np.sqrt( (x1-x)*(x1-x)+(y1-y)*(y1-y) ) + np.sqrt( (x2-x)*(x2-x)+(y2-y)*(y2-y) ) < a           

allFnames=getListOfStringsFromConfigs(cfgTxt,"#FNAMES_BEG","#FNAMES_END")
foutName=getValueFromConfigs(cfgTxt,"OutpuFileName","fggHists.root")
mlScoreFileNames=getValueFromConfigs(cfgTxt,"mlScoreFileName","ats.plk").split(",")
mlScoreTags=getValueFromConfigs(cfgTxt,"mlScoreTag","").split(",")
processID=getValueFromConfigs(cfgTxt,"processID",default="DATA")
treeName=getValueFromConfigs(cfgTxt,"treeName",default="tagsDumper/trees/Data_13TeV_TrippleHTag_0")
outTreeName=getValueFromConfigs(cfgTxt,"outTreeName",default="Data_13TeV_TrippleHTag_0")
eventIdxOffset=int(getValueFromConfigs(cfgTxt,"EventIdxOffset",default=0))
doOverlapRemoval =getValueFromConfigs(cfgTxt,"doOverlapRemoval",default="1") ; doOverlapRemoval = int(doOverlapRemoval) > 0.5
isData =getValueFromConfigs(cfgTxt,"isData",default="1") ; isData = int(isData) > 0.5
etaMax =float(getValueFromConfigs(cfgTxt,"etaMax",default="2.5"))
pTMin =float(getValueFromConfigs(cfgTxt,"pTMin",default="25.0"))
overlapRemovalDRMax =float(getValueFromConfigs(cfgTxt,"overlapRemovalDRMax",default="0.4"))
skipBad=getValueFromConfigs(cfgTxt,"skipBad",default="0") ; skipBad = int(skipBad) > 0.5
mlScoreFileNameMap={ tagM : val for tagM,val in zip(mlScoreTags,mlScoreFileNames) }

etaMax=2.5
pTMin=25.0

print("allFnames   :  ",              allFnames)
print("foutName   :  ",               foutName)
print("processID   :  ",              processID)
print("treeName   :  ",               treeName)
print("eventIdxOffset   :  ", eventIdxOffset)
print("ML Score file    :  ", mlScoreFileNames)
print("ML Score tags    :  ", mlScoreTags)

maxEvents=-1
tmp_=getValueFromConfigs(cfgTxt,"MaxEvents")
if tmp_!='':
    maxEvents=int(tmp_)

for i in allFnames:
    print(" file : ",i)
if(maxEvtsSuperSeeder > 0):
    maxEvents=maxEvtsSuperSeeder
print("maxevents : ",maxEvents)


mlScoreInterface={}
for ky in mlScoreFileNameMap:
    if os.path.isfile(mlScoreFileNameMap[ky]):
        mlScoreInterface[ky]=hhhMLI.mlScoreManger(mlScoreFileNameMap[ky])
        print("Loaded ML Sores for ",mlScoreInterface[ky].GetEntries()," events ")
        print()
    #    mlScoreInterface.printAllEventsIn()
    else:
        print("ML Score file not found !! ",mlScoreFileNameMap[ky])
        exit(1)

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
branchesToFill.append('nVldJets')
branchesToFill.append('matchedJet_b1')
branchesToFill.append('matchedJet_b2')
branchesToFill.append('matchedJet_b3')
branchesToFill.append('matchedJet_b4')

branchesToFill.append('preselectedAt15GeV')
branchesToFill.append('preselectedAt20GeV')
branchesToFill.append('preselectedAtOLR')

for i in range(8):
    branchesToFill.append('label_'+str(i))
for i in range(8):
    branchesToFill.append('label_idx_'+str(i))

for i in range(8):
    branchesToFill.append('jet_'+str(i)+'_drWithDP')
    branchesToFill.append('jet_'+str(i)+'_drWithDP_leadG')
    branchesToFill.append('jet_'+str(i)+'_dr_WithDPSubleadG')
    branchesToFill.append('jet_'+str(i)+'_dEta_WithDP')
    branchesToFill.append('jet_'+str(i)+'_dEta_WithDPLeadG')
    branchesToFill.append('jet_'+str(i)+'_dEta_WithDPSubleadG')
    branchesToFill.append('jet_'+str(i)+'_dPhi_WithDP')
    branchesToFill.append('jet_'+str(i)+'_dPhi_WithDPLeadG')
    branchesToFill.append('jet_'+str(i)+'_dPhi_WithDPSubleadG')
    branchesToFill.append('jet_'+str(i)+'_mass_WithDP')
    branchesToFill.append('jet_'+str(i)+'_mass_WithDPLeadG')
    branchesToFill.append('jet_'+str(i)+'_mass_WithDPSubleadG')
    branchesToFill.append('jet_'+str(i)+'_massRatio_WithDP')
    branchesToFill.append('jet_'+str(i)+'_ptRatio_WithDP')
    branchesToFill.append('jet_'+str(i)+'_ptRatio_WithDPLeadG')
    branchesToFill.append('jet_'+str(i)+'_ptRatio_WithDPSubleadG')
    for mlScoreTag in mlScoreInterface:
        branchesToFill.append('jet_'+str(i)+mlScoreTag+'_y0')
        branchesToFill.append('jet_'+str(i)+mlScoreTag+'_y1')
        branchesToFill.append('jet_'+str(i)+mlScoreTag+'_y0s')
        branchesToFill.append('jet_'+str(i)+mlScoreTag+'_y1s')
        branchesToFill.append('jet_'+str(i)+mlScoreTag+'_score')

branches=np.unique(branchesToFill)

filemode="RECREATE"
fout = ROOT.TFile(foutName, filemode)
dir_ = fout.mkdir('tagsDumper')
dir_ = dir_.mkdir('trees')

tofill = OrderedDict(zip(branches, [np.nan]*len(branches)))

ntuple={}
outputDataDict=OrderedDict(zip(branches, [np.nan]*len(branches))) 
ntuple['all'] = ROOT.TNtuple(outTreeName, outTreeName, ':'.join(branches))


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

eventIdx=eventIdxOffset
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
        if(evt%100==0):
            print("   Doing i = ",evt," / ",maxEvents_,
                  " n Found : ", nMatchFound,
                  " n Miss  : ", nMiss)
            print(" gg mass , h1MassVsh2Mass  : ",eTree.CMS_hgg_mass ,"( ",eTree.M1jj , eTree.M2jj,")" )
            
        wei=eTree.weight
        for ky in tofill:
            if ky in allBranches:
                tofill[ky]=getattr(eTree,ky)
            else:
                tofill[ky]=-999e3
        
        tofill['matchedJet_b1']=-1
        tofill['matchedJet_b2']=-1
        tofill['matchedJet_b3']=-1
        tofill['matchedJet_b4']=-1
        nVldJets=0
        nVldJets2=0
        
        ##   JET PRE-SELECTION
        jetMask=hhhSelector.getSelectedJetCollectionMaskEta(eTree,etaMax=etaMax)

        if sum(jetMask) < 4 :
            sumEntries.Fill('nJetPreselectionEta',1)
            sumWeights.Fill('nJetPreselectionEta',wei)
            if skipBad: continue

        tofill['preselectedAt20GeV']=1.0
        jetMask=hhhSelector.getSelectedJetCollectionMaskPt(eTree,jetMask=jetMask,pTMin=20)
        if sum(jetMask) < 4 :
            sumEntries.Fill('nJetPreselectionPt',1)
            sumWeights.Fill('nJetPreselectionPt',wei)
            tofill['preselectedAt20GeV']=0.0
            if skipBad: continue

        tofill['preselectedAt15GeV']=1.0
        jetMask=hhhSelector.getSelectedJetCollectionMaskPt(eTree,jetMask=jetMask,pTMin=15)
        if sum(jetMask) < 4 :
            sumEntries.Fill('nJetPreselectionPt',1)
            sumWeights.Fill('nJetPreselectionPt',wei)
            tofill['preselectedAt15GeV']=0.0
            if skipBad: continue

        tofill['preselectedAtOLR']=1.0
        if doOverlapRemoval:
            jetMask=hhhSelector.getSelectedJetCollectionMaskOverLap(eTree,jetMask=jetMask,overlapRemovalDRMax=overlapRemovalDRMax)
            if sum(jetMask) < 4 :
                sumEntries.Fill('nJetPreselectionOR',1)
                sumWeights.Fill('nJetPreselectionOR',wei)
                tofill['preselectedAtOLR']=0.0
                if skipBad: continue

        for i in range(8):
            tofill['label_'+str(i)]=0
            tofill['label_idx_'+str(i)]=0
            if tofill['jet_'+str(i)+'_isValid'] > 0.5:
                nVldJets+=1
           

        tofill['nVldJets']=nVldJets
       
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

        
        eventIdx=eTree.event
        for mlScoreTag in mlScoreInterface:
            rslt=mlScoreInterface[mlScoreTag].getEntry(eventIdx)
            if rslt['valid']:
                for i in range(8):
                    tofill['jet_'+str(i)+mlScoreTag+'_y0']=rslt['y0'][i] 
                    tofill['jet_'+str(i)+mlScoreTag+'_y1']=rslt['y1'][i] 
                    tofill['jet_'+str(i)+mlScoreTag+'_y0s']=hhhMLI.sigmoid(rslt['y0'][i] )
                    tofill['jet_'+str(i)+mlScoreTag+'_y1s']=hhhMLI.sigmoid(rslt['y1'][i] )
                    tofill['jet_'+str(i)+mlScoreTag+'_score']=rslt['score'][i] 
                    #print(mlScoreTag," i ",i," ", rslt['score'][i] )
                    
                    for i in range(NJETS):
                        tofill['label_'+str(i)]=float(rslt['label'][i]  >0.1)
                        tofill['label_idx_'+str(i)]=rslt['label'][i]
                        if rslt['label'][i]>0 and rslt['label'][i] <5:
                            k=int(rslt['label'][i])
                            tofill['matchedJet_b'+str(k)]=i
        
                        ## SANITY CHECK 
                        #if  tofill['jet_'+str(i)+'_isValid'] != rslt['vldMask'][i]:
                        #    x=[]
                        #    for kk in range(NJETS):
                        #        x.append((tofill['jet_'+str(i)+'_isValid']))
                        #    print("FISHY FISY FISHY !! idx ",evt," eventIdx : ",eventIdx," vldMask from Json ",rslt['vldMask']," | from tree ",x)
                nMatchFound+=1
            else:
                nMiss+=1
                print("Event not found ! " , eventIdx," nvld Jets : ",nVldJets)

        for i in outputDataDict:
            outputDataDict[i]=tofill[i]
        
        cat='cat0'
        sumWeights.Fill('total_'+cat, 1)
        sumWeights.Fill('total_wei'+cat, wei)
        ntuple['all'].Fill(array('f', outputDataDict.values()))

    simFile.Close()           
    print("Closing file : ",fname)
print("Number of matches found ",nMatchFound)
print("Number of misses ",nMiss)

dir_.cd()    
sumWeights.Write()
m1m2.Write()

for cat in ntuple:
    ntuple[cat].Write()
fout.Close()
print(" File written out  : ",foutName)

