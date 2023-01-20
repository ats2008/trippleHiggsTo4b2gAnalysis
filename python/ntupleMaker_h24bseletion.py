from __future__ import print_function
from collections import OrderedDict
import ROOT 
import numpy as np
import trippleHiggsUtils as hhhUtil
from Util  import *
from branches import *
from array import array
import trippleHiggsSelector  as hhhSelector
from TMVA_Model import *

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
processID=getValueFromConfigs(cfgTxt,"processID",default="DATA")
treeName=getValueFromConfigs(cfgTxt,"treeName",default="tagsDumper/trees/Data_13TeV_TrippleHTag_0")
outTreeName=getValueFromConfigs(cfgTxt,"outTreeName",default="Data_13TeV_TrippleHTag_0")
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
resetWeight=float(getValueFromConfigs(cfgTxt,"resetWeight",default=-1e5))
doSR=getBoolFromConfigs(cfgTxt,"doSR",default=False)
etaMax =float(getValueFromConfigs(cfgTxt,"etaMax",default="2.5"))
pTMin =float(getValueFromConfigs(cfgTxt,"pTMin",default="25.0"))
doOverlapRemoval =getValueFromConfigs(cfgTxt,"doOverlapRemoval",default="1") ; doOverlapRemoval = int(doOverlapRemoval) > 0.5
overlapRemovalDRMax =float(getValueFromConfigs(cfgTxt,"overlapRemovalDRMax",default="0.4"))

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
print("doSR     :       ", doSR)
print('etaMax   :       ', etaMax)
print('pTMin    :       ', pTMin)
print('drMax    :       ', drMax)

scaler=hhhUtil.mcScaler()
if doPtReWeighting:
    scaler.setSFHistFromFile(pTReweitingFile,pTReweitingHistName)
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


nNoHtoGammaGamma=0
nNoHHto4B=0
totalEvents=0
branchesToFill=allMCBranches
branchesToFill.append("nonResonantMVA_v0")
branchesToFill.append("peakingMVA_v0")
branchesToFill.append("nonResonantMVA_v1")
branchesToFill.append("nonResonantMVA_v2")
branchesToFill.append("peakingMVA_v1")

branches=np.unique(branchesToFill)
filemode="RECREATE"
foutName=foutName.replace('.root',ext+'.root')
fout = ROOT.TFile(foutName, filemode)
dir_ = fout.mkdir('tagsDumper')
dir_ = dir_.mkdir('trees')

branches_skimmed=['M1jj','M2jj','CMS_hgg_mass','peakingMVA_v1','nonResonantMVA_v1','dZ','r_HH','weight','nonResonantMVA_v2']

ntuple={}
ntuple['all'] = ROOT.TNtuple(outTreeName+'_NOTAG', outTreeName, ':'.join(branches_skimmed))
tofill = OrderedDict(zip(branches, [np.nan]*len(branches)))

skimmedDataDict = OrderedDict(zip(branches_skimmed, [np.nan]*len(branches_skimmed)))
kyList = [ ky for ky in tofill ]

print("len(branches) : " , len(branches))
m1m2 = ROOT.TH2F("m1jj_m2jj","H1bb , H2bb mass",300,0.0,300.0,300,0.0,300. )
mvaAll=ROOT.TH1F("mvaPk_v1_all","",24,-0.1,1.1)
mvaPass=ROOT.TH1F("mvaPk_v1_pass","",24,-0.1,1.1)


edges=np.array([0.03099544, 0.63089365, 0.72614938, 0.77085077, 0.80215021,
       0.82709223, 0.84269528, 0.85683767, 0.86924267, 0.88029956,
       0.88964659, 0.89783251, 0.90468966, 0.91155702, 0.91749775,
       0.92227721, 0.92692335, 0.93142081, 0.93516294, 0.93902266,
       0.94250536, 0.94561276, 0.94860082, 0.95128779, 0.95418559,
       0.95662987, 0.95925509, 0.96143816, 0.96357881, 0.96554942,
       0.96749735, 0.96924451, 0.97082723, 0.97247566, 0.97380439,
       0.97518504, 0.97649124, 0.97777721, 0.97895595, 0.98027201,
       0.98148018, 0.98266215, 0.98377918, 0.98489203, 0.98616692,
       0.98737586, 0.98870801, 0.98988895, 0.99142853, 0.99329859])

def getCMVAScore(x,edges):
    if x > edges[-1]:
        return 1.0
    n=edges.searchsorted(x)
    return (n-1.0*( edges[n] - x)/(edges[n] - edges[n-1])  )/edges.shape[0]

sumEntries=ROOT.TH1F("sumEvts","sumEvts",1,0.0,1.0)
sumEntries.SetCanExtend(ROOT.TH1.kAllAxes)
sumWeights=ROOT.TH1F("sumWeights","sumWeights",1,0.0,1.0)
sumWeights.SetCanExtend(ROOT.TH1.kAllAxes)

th1Store={}
th1Store["allGenDR_H1bb"]= ROOT.TH1F("allGenDR_H1bb","",60,0.0,6.0)
th1Store["allGenDR_H2bb"]= ROOT.TH1F("allGenDR_H2bb","",60,0.0,6.0)
th1Store["allGenPt_H1bb"]= ROOT.TH1F("allGenPt_H1bb","",600,0.0,600.0)
th1Store["allGenPt_H2bb"]= ROOT.TH1F("allGenPt_H2bb","",600,0.0,600.0)
th1Store["allGenEta_H1bb"]= ROOT.TH1F("allGenEta_H1bb","",100,-5.0,5.0)
th1Store["allGenEta_H2bb"]= ROOT.TH1F("allGenEta_H2bb","",100,-5.0,5.0)


th1Store["noGenMatchMass_H1bb"]= ROOT.TH1F("noGenMatchMass_H1bb","",600,0.0,300.0)
th1Store["noGenMatchMass_H2bb"]= ROOT.TH1F("noGenMatchMass_H2bb","",600,0.0,300.0)
th1Store["recoMatchMass_H1bb"]= ROOT.TH1F("recoMatchMass_H1bb","",600,0.0,300.0)
th1Store["recoMatchMass_H2bb"]= ROOT.TH1F("recoMatchMass_H2bb","",600,0.0,300.0)
th1Store["quadRecoMass_H1bb"]= ROOT.TH1F("quadRecoMass_H1bb","",600,0.0,300.0)
th1Store["quadRecoMass_H2bb"]= ROOT.TH1F("quadRecoMass_H2bb","",600,0.0,300.0)
th1Store["quadV2RecoMass_H1bb"]= ROOT.TH1F("quadV2RecoMass_H1bb","",600,0.0,300.0)
th1Store["quadV2RecoMass_H2bb"]= ROOT.TH1F("quadV2RecoMass_H2bb","",600,0.0,300.0)
th1Store["noGenMatch_H1bb_b0HFlav"]= ROOT.TH1F("noGenMatch_H1bb_b0HFlav","",7,-0.5,6.5)
th1Store["noGenMatch_H1bb_b1HFlav"]= ROOT.TH1F("noGenMatch_H1bb_b1HFlav","",7,-0.5,6.5)
th1Store["noGenMatch_H2bb_b2HFlav"]= ROOT.TH1F("noGenMatch_H2bb_b2HFlav","",7,-0.5,6.5)
th1Store["noGenMatch_H2bb_b3HFlav"]= ROOT.TH1F("noGenMatch_H2bb_b3HFlav","",7,-0.5,6.5)
th1Store["recoMatch_H1bb_b0HFlav"]= ROOT.TH1F("recoMatch_H1bb_b0HFlav","",7,-0.5,6.5)
th1Store["recoMatch_H1bb_b1HFlav"]= ROOT.TH1F("recoMatch_H1bb_b1HFlav","",7,-0.5,6.5)
th1Store["recoMatch_H2bb_b2HFlav"]= ROOT.TH1F("recoMatch_H2bb_b2HFlav","",7,-0.5,6.5)
th1Store["recoMatch_H2bb_b3HFlav"]= ROOT.TH1F("recoMatch_H2bb_b3HFlav","",7,-0.5,6.5)
th1Store["quadReco_b0HFlav"]= ROOT.TH1F("quadReco_b0HFlav","",7,-0.5,6.5)
th1Store["quadReco_b1HFlav"]= ROOT.TH1F("quadReco_b1HFlav","",7,-0.5,6.5)
th1Store["quadReco_b2HFlav"]= ROOT.TH1F("quadReco_b2HFlav","",7,-0.5,6.5)
th1Store["quadReco_b3HFlav"]= ROOT.TH1F("quadReco_b3HFlav","",7,-0.5,6.5)
th1Store["quadV2Reco_b0HFlav"]= ROOT.TH1F("quadV2Reco_b0HFlav","",7,-0.5,6.5)
th1Store["quadV2Reco_b1HFlav"]= ROOT.TH1F("quadV2Reco_b1HFlav","",7,-0.5,6.5)
th1Store["quadV2Reco_b2HFlav"]= ROOT.TH1F("quadV2Reco_b2HFlav","",7,-0.5,6.5)
th1Store["quadV2Reco_b3HFlav"]= ROOT.TH1F("quadV2Reco_b3HFlav","",7,-0.5,6.5)

th1Store["recoMatchMass_H1H2"]= ROOT.TH2F("recoMatchMass_H1H2","",150,0.0,300.0,150,0.0,300.0)
th1Store["quadRecoMass_H1H2"] = ROOT.TH2F("quadRecoMass_H1H2" ,"",150,0.0,300.0,150,0.0,300.0)
th1Store["quadV2RecoMass_H1H2"] = ROOT.TH2F("quadV2RecoMass_H1H2" ,"",150,0.0,300.0,150,0.0,300.0)


th1Store['acceptance2d']=ROOT.TH1F("acceptance2d","acceptance2d",1,0.0,1.0)
th1Store['acceptance2d'].SetCanExtend(ROOT.TH1.kAllAxes)

for i in range(4):
    th1Store["genMatchDr_"+str(i)]= ROOT.TH1F("genMatchDr_"+str(i),"",120,0.0,6.0)
    th1Store["genMatchWithFlavourDr_"+str(i)]= ROOT.TH1F("genMatchWithFlavourDr_"+str(i),"",120,0.0,6.0)
    th1Store["recoMatchDr_"+str(i)]= ROOT.TH1F("recoMatchDr_"+str(i),"",120,0.0,6.0)
    th1Store["noGenMatchJetMass_"+str(i)]= ROOT.TH1F("noGenMatchJetMass_"+str(i),"",600,0.0,300.0)
    th1Store["recoMatchJetMass_"+str(i)]= ROOT.TH1F("recoMatchJetMass_"+str(i),"",600,0.0,300.0)
    th1Store["noGenMatchDeepJet_"+str(i)]= ROOT.TH1F("noGenMatchDeepJet_"+str(i),"",50,-2.5,2.5)
    th1Store["recoMatchDeepJet_"+str(i)]= ROOT.TH1F("recoMatchDeepJet_"+str(i),"",50,-2.5,2.5)
    th1Store["noGenMatchHFlavour_"+str(i)]= ROOT.TH1F("noGenMatchHFlavour_"+str(i),"",7,-0.5,6.5)
    th1Store["recoMatchHFlavour_"+str(i)]= ROOT.TH1F("recoMatchHFlavour_"+str(i),"",7,-0.5,6.5)
    th1Store["noGenMatchJetIdIsLoose_"+str(i)]= ROOT.TH1F("noGenMatchJetIdIsLoose_"+str(i),"",5,-2.5,2.5)
    th1Store["recoMatchJetIdIsLoose_"+str(i)]= ROOT.TH1F("recoMatchJetIdIsLoose_"+str(i),"",5,-2.5,2.5)
    th1Store["noGenMatchJetIdIsTight_"+str(i)]= ROOT.TH1F("noGenMatchJetIdIsTight_"+str(i),"",5,-2.5,2.5)
    th1Store["recoMatchJetIdIsTight_"+str(i)]= ROOT.TH1F("recoMatchJetIdIsTight_"+str(i),"",5,-2.5,2.5)
    th1Store["noGenMatchJetIdIsTight2017_"+str(i)]= ROOT.TH1F("noGenMatchJetIdIsTight2017_"+str(i),"",5,-2.5,2.5)
    th1Store["recoMatchJetIdIsTight2017_"+str(i)]= ROOT.TH1F("recoMatchJetIdIsTight2017_"+str(i),"",5,-2.5,2.5)
    th1Store["noGenMatchJetIdIsTight2018_"+str(i)]= ROOT.TH1F("noGenMatchJetIdIsTight2018_"+str(i),"",5,-2.5,2.5)
    th1Store["recoMatchJetIdIsTight2018_"+str(i)]= ROOT.TH1F("recoMatchJetIdIsTight2018_"+str(i),"",5,-2.5,2.5)
    th1Store["noGenMatchJetIdPUScore_"+str(i)]= ROOT.TH1F("noGenMatchJetIdPUScore_"+str(i),"",25,-1.0,1.0)
    th1Store["recoMatchJetIdIPUScore_"+str(i)]= ROOT.TH1F("recoMatchJetIdIPUScore_"+str(i),"",25,-1.0,1.0)

    th1Store["allGenPt_"+str(i)]= ROOT.TH1F("allGenPt_"+str(i),"",200,0.0,500.0)
    th1Store["recoGenPt_"+str(i)]= ROOT.TH1F("recoGenPt_"+str(i),"",200,0.0,500.0)
    th1Store["recoPt_"+str(i)]= ROOT.TH1F("recoPt_"+str(i),"",200,0.0,500.0)
    th1Store["quadRecoPt_"+str(i)]= ROOT.TH1F("quadRecoPt_"+str(i),"",200,0.0,500.0)
    th1Store["quadV2RecoPt_"+str(i)]= ROOT.TH1F("quadV2RecoPt_"+str(i),"",200,0.0,500.0)
    
    th1Store["allGenEta_"+str(i)]= ROOT.TH1F("allGenEta_"+str(i),"",100,-5.0,5.0)
    th1Store["recoGenEta_"+str(i)]= ROOT.TH1F("recoGenEta_"+str(i),"",100,-5.0,5.0)
    th1Store["recoEta_"+str(i)]= ROOT.TH1F("recoEta_"+str(i),"",100,-5.0,5.0)
    th1Store["quadRecoEta_"+str(i)]= ROOT.TH1F("quadRecoEta_"+str(i),"",100,-5.0,5.0)
    th1Store["quadV2RecoEta_"+str(i)]= ROOT.TH1F("quadV2RecoEta_"+str(i),"",100,-5.0,5.0)
    
    th1Store["fggIndex_matched_b"+str(i)]= ROOT.TH1F("fggIndex_matched_b"+str(i),"",11,-2.5,8.5)
    th1Store["btagOrderIndex_matched_b"+str(i)]= ROOT.TH1F("btagOrderIndex_matched_b"+str(i),"",11,-2.5,8.5)
    th1Store["pTOrderIndex_matched_b"+str(i)]  = ROOT.TH1F("pTorderIndex_matched_b"+str(i),"",11,-2.5,8.5)


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
    for evt in range(maxEvents_):
        i=evt
        #print("Event :    ",evt)
        if(i%250==0):
            print("   Doing i = ",i," / ",maxEvents_,
                  " n No HiggsCand : ", nNoHtoGammaGamma,
                  " n No 4BCands  : ", nNoHHto4B)
            print(" gg mass , h1MassVsh2Mass  : ",eTree.CMS_hgg_mass ,"( ",eTree.M1jj , eTree.M2jj,")" )
            
        eTree.GetEntry(evt)
        wei=eTree.weight
        sumEntries.Fill('total', 1)
        sumWeights.Fill('total', wei)
        ## Control and signal region
        allGammaDaus=hhhSelector.getHiggsDauP4s(eTree,22)
        allBDaus=hhhSelector.getHiggsDauP4s(eTree,5)
        
        
        #if doSR and isSR:
        #    continue
        #if( (allGammaDaus[0]+allGammaDaus[1]).Pt() < 100 ): continue;
        text=''
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
        if sum(jetMask) < 1:
            continue
        rslt=hhhSelector.getBestGetMatchesGlobalFGG(eTree,drMax=drMax,jetMask=jetMask)
        idxs=rslt['idxs']
        fgg_idxs=rslt['fgg_idxs']
        drMins=rslt['drMins']
        drMinsPostFlavourMatch=rslt['drMinsPostFlavourMatch']
        jetMass=rslt['jetMass']
        allJetMatchP4=rslt['p4s']
        deepJetScores=rslt['deepJetScores']
        hFlavour=rslt['hFlavour']
        jetIdLoose = rslt['isLoose']
        jetIdTight = rslt['isTight']
        jetIdTight2017 = rslt['isTight2017']
        jetIdTight2018 = rslt['isTight2018']
        hFlavourAll= rslt['hFlavourAll']
        btagRank= rslt['bJetScoreRanks']
        isAllRecoed=True
        isRecoed={}
        jetMatchesString=""
        for i in range(4):
            if idxs[i] < 0:
                jetMatchesString+=" b"+str(i+1)+"(j"+u'\u2717'+")" 
            else:
                jetMatchesString+=" b"+str(i+1)+"(j"+str(fgg_idxs[i])+")"
        
        #print("fgg_idxs : ",fgg_idxs)
        #print("drMinsPostFlavourMatch : ",drMinsPostFlavourMatch)
        #print("drMins : ",drMins)
        #print("hFlavour : ",hFlavour)
        #print("hFlavourAll : ",hFlavourAll)
        #text='Event '+str(evt)+" , "+jetMatchesString
        #fpicname='results/plots/genPlots4bSel/'+str(evt)+'.png'
        #print(fpicname)
        #hhhUtil.visualizeEvents(eTree,text=text,outputFname=fpicname,jetMask=jetMask)
        #print("Event : ",evt)
        #print("Idxs matched  : ",idxs,fgg_idxs)
        #print(drMinsPostFlavourMatch)
        for i in range(4):
            isRecoed['isRecoed_'+str(i)]=True

            th1Store["allGenPt_"+str(i)].Fill( allBDaus[i].Pt()  )
            th1Store["allGenEta_"+str(i)].Fill( allBDaus[i].Eta() )
            if len(drMins) >0:
                th1Store["genMatchDr_"+str(i)].Fill( drMins[i])
                if hFlavourAll[i]> 4.8:
                    th1Store["genMatchWithFlavourDr_"+str(i)].Fill( drMinsPostFlavourMatch[i])
            th1Store['fggIndex_matched_b'+str(i)].Fill(fgg_idxs[i])
            th1Store['pTOrderIndex_matched_b'+str(i)].Fill(idxs[i])
            th1Store['btagOrderIndex_matched_b'+str(i)].Fill(btagRank[i])
            #histStore['pTOrderIndex_matched_b'+str(i)].Fill(-2)
            if idxs[i] < 0:
         #       print("\t Filling no Match for i = ",i)
                isAllRecoed=False
                jetMatchesString+=" b"+str(i+1)+"(j"+u'\u2717'+")" 
                isRecoed['isRecoed_'+str(i)]=False
                th1Store["noGenMatchJetMass_"+str(i)].Fill( jetMass[i] )
                th1Store["noGenMatchDeepJet_"+str(i)].Fill( deepJetScores[i] )
                th1Store["noGenMatchHFlavour_"+str(i)].Fill( hFlavour[i] )
                th1Store["noGenMatchJetIdIsLoose_"+str(i)].Fill( jetIdLoose[i] )
                th1Store["noGenMatchJetIdIsTight_"+str(i)].Fill( jetIdTight[i] )
                th1Store["noGenMatchJetIdIsTight2017_"+str(i)].Fill( jetIdTight2017[i] )
                th1Store["noGenMatchJetIdIsTight2018_"+str(i)].Fill( jetIdTight2018[i] )
            else:
                jetMatchesString+=" b"+str(i+1)+"(j"+str(idxs[i])+")"
                sumEntries.Fill('isRecoed_'+str(i), 1)
                sumWeights.Fill('isRecoed_'+str(i), wei)
                th1Store["recoMatchDr_"+str(i)].Fill( drMins[i])
                th1Store["recoMatchJetMass_"+str(i)].Fill( jetMass[i] )
                th1Store["recoMatchDeepJet_"+str(i)].Fill( deepJetScores[i] )
                th1Store["recoMatchHFlavour_"+str(i)].Fill( hFlavour[i] )
                th1Store["recoMatchJetIdIsLoose_"+str(i)].Fill( jetIdLoose[i] )
                th1Store["recoMatchJetIdIsTight_"+str(i)].Fill( jetIdTight[i] )
                th1Store["recoMatchJetIdIsTight2017_"+str(i)].Fill( jetIdTight2017[i] )
                th1Store["recoMatchJetIdIsTight2018_"+str(i)].Fill( jetIdTight2018[i] )
                th1Store["recoGenPt_"+str(i)].Fill( allBDaus[i].Pt()  )
                th1Store["recoGenEta_"+str(i)].Fill( allBDaus[i].Eta() )
                th1Store["recoPt_"+str(i) ].Fill( getattr(eTree,'jet_'+str(fgg_idxs[i])+'_pt')  )
                th1Store["recoEta_"+str(i)].Fill( getattr(eTree,'jet_'+str(fgg_idxs[i])+'_eta') )

        th1Store["allGenPt_H1bb"].Fill((allBDaus[0]+allBDaus[1]).Pt() )        
        th1Store["allGenPt_H2bb"].Fill((allBDaus[2]+allBDaus[3]).Pt() )        
        th1Store["allGenEta_H1bb"].Fill((allBDaus[0]+allBDaus[1]).Eta() )        
        th1Store["allGenEta_H2bb"].Fill((allBDaus[2]+allBDaus[3]).Eta() )        
        th1Store["allGenDR_H1bb"].Fill(allBDaus[0].DeltaR(allBDaus[1]) )        
        th1Store["allGenDR_H2bb"].Fill(allBDaus[2].DeltaR(allBDaus[3]) )        
         
        nout=0
        if not isRecoed['isRecoed_0']:            nout+=1;    
        if not isRecoed['isRecoed_1']:            nout+=1;    
        if not isRecoed['isRecoed_2']:            nout+=1;    
        if not isRecoed['isRecoed_3']:            nout+=1;
        
        nleadOut=0
        if not isRecoed['isRecoed_0']:            nleadOut+=1;    
        if not isRecoed['isRecoed_2']:            nleadOut+=1;    
        
        nSubLeadOut=0
        if not isRecoed['isRecoed_1']:            nSubLeadOut+=1;    
        if not isRecoed['isRecoed_3']:            nSubLeadOut+=1;    
        
        h1Out= not ( isRecoed['isRecoed_0'] and isRecoed['isRecoed_1'])
        h2Out= not ( isRecoed['isRecoed_2'] and isRecoed['isRecoed_3'])
        if not h1Out:
            jetMatchesString+="| H1("+u'\u2713'+")"
        else:
            jetMatchesString+="| H1("+u'\u2717'+")"
        if not h2Out:
            jetMatchesString+=" H2("+u'\u2713'+")"
        else:
            jetMatchesString+=" H2("+u'\u2717'+")"
        
        
        sumEntries.Fill(str(nout)+'_JetsNotIsRecoed',1)
        sumEntries.Fill(str(nleadOut)+'_LeadNotIsRecoed',1)
        sumEntries.Fill(str(nSubLeadOut)+'_SubLeadNotIsRecoed',1)
        if h1Out:
            sumEntries.Fill('H1NotIsRecoed', 1)
            th1Store["noGenMatchMass_H1bb"].Fill( (allJetMatchP4[0] + allJetMatchP4[1]).M()  )
            th1Store["noGenMatch_H1bb_b0HFlav"].Fill( hFlavour[0]  )
            th1Store["noGenMatch_H1bb_b1HFlav"].Fill( hFlavour[1]  )
        else:
            th1Store["recoMatchMass_H1bb"].Fill( (allJetMatchP4[0] + allJetMatchP4[1]).M()  )
            th1Store["recoMatch_H1bb_b0HFlav"].Fill( hFlavour[0]  )
            th1Store["recoMatch_H1bb_b1HFlav"].Fill( hFlavour[1]  )
            
        if h2Out:
            sumEntries.Fill('H2NotIsRecoed', 1)
            th1Store["noGenMatchMass_H2bb"].Fill( (allJetMatchP4[2] + allJetMatchP4[3]).M()  )
            th1Store["noGenMatch_H2bb_b2HFlav"].Fill( hFlavour[2]  )
            th1Store["noGenMatch_H2bb_b3HFlav"].Fill( hFlavour[3]  )
        else:
            th1Store["recoMatchMass_H2bb"].Fill( (allJetMatchP4[2] + allJetMatchP4[3]).M()  )
            th1Store["recoMatch_H2bb_b2HFlav"].Fill( hFlavour[2]  )
            th1Store["recoMatch_H2bb_b3HFlav"].Fill( hFlavour[3]  )

        if isRecoed['isRecoed_0'] and isRecoed['isRecoed_1'] and isRecoed['isRecoed_2'] and isRecoed['isRecoed_3'] :
            sumEntries.Fill('allBJetsIsRecoed', 1)
            th1Store["recoMatchMass_H1H2"].Fill( (allJetMatchP4[0] + allJetMatchP4[1]).M() , (allJetMatchP4[2] + allJetMatchP4[3]).M() )
        else:
            pass


        allQuads=hhhSelector.getBJetParisFGG(eTree,mask=jetMask)
        #print("\nRecoed bjets : ",idxs)
        if allQuads['isValid']:
            dr = allQuads['bJetQuad']['r_HH']
            isSR = dr < 25.0
            if isSR:
                sumEntries.Fill('isROI', 1)
                sumWeights.Fill('isROI', wei)

            sumEntries.Fill('hasQuad', 1)
            #print("selection : ",allQuads['bJetQuad']['fgg_idxs'])
            th1Store["quadRecoMass_H1bb"].Fill(allQuads['bJetQuad']['m1'])
            th1Store["quadRecoMass_H2bb"].Fill(allQuads['bJetQuad']['m2'])
            th1Store["quadRecoMass_H1H2"].Fill(allQuads['bJetQuad']['m1'], allQuads['bJetQuad']['m2'])

            hasAllGenMatchedBJetsInQuad=True
            for i in range(4):
                ix=allQuads['bJetQuad']['fgg_idxs'][i]
                th1Store["quadRecoPt_"+str(i)].Fill( getattr(eTree,'jet_'+str(ix)+'_pt')  )
                th1Store["quadRecoEta_"+str(i)].Fill( getattr(eTree,'jet_'+str(ix)+'_eta')  )
                th1Store["quadReco_b"+str(i)+"HFlav"].Fill( getattr(eTree,'jet_'+str(ix)+'_flavour')  )

                if ix not in fgg_idxs:
                    hasAllGenMatchedBJetsInQuad=False
            if hasAllGenMatchedBJetsInQuad:
                sumEntries.Fill("hasAllGenMatchedBJetsInQuad",1)

            hasH1GenMatchedInQuad=False
            if allQuads['bJetQuad']['fgg_idxs'][0] in fgg_idxs[0:2] and allQuads['bJetQuad']['fgg_idxs'][1] in  fgg_idxs[0:2]:
                hasH1GenMatchedInQuad=True
            if allQuads['bJetQuad']['fgg_idxs'][0] in fgg_idxs[2:] and allQuads['bJetQuad']['fgg_idxs'][1] in  fgg_idxs[2:]:
                hasH1GenMatchedInQuad=True
            if hasH1GenMatchedInQuad:
                sumEntries.Fill("hasH1GenMatchedInQuad",1)
            
            hasH2GenMatchedInQuad=False
            if allQuads['bJetQuad']['fgg_idxs'][2] in fgg_idxs[0:2] and allQuads['bJetQuad']['fgg_idxs'][3] in  fgg_idxs[0:2]:
                hasH2GenMatchedInQuad=True
            if allQuads['bJetQuad']['fgg_idxs'][2] in fgg_idxs[2:] and allQuads['bJetQuad']['fgg_idxs'][3] in  fgg_idxs[2:]:
                hasH2GenMatchedInQuad=True
            if hasH2GenMatchedInQuad:
                sumEntries.Fill("hasH2GenMatchedInQuad",1)
            
            hasAllJetsMatchedProperly=True
            num=0
            for k in fgg_idxs:
                if k in  allQuads['bJetQuad']['fgg_idxs']:
                    num+=1
            sumEntries.Fill('hasMatchInQuad_'+str(num),1)

            for i in range(4):
                if fgg_idxs[i]!=allQuads['bJetQuad']['fgg_idxs'][i]:
                    hasAllJetsMatchedProperly=False
                    break
            if hasAllJetsMatchedProperly:
                sumEntries.Fill("hasAllJetsMatchedOneToOneProperly",1)

            if hasH2GenMatchedInQuad and hasH1GenMatchedInQuad:
                sumEntries.Fill("hasEachHMatchedProperly",1)
        else:
            sumEntries.Fill('noQuad', 1)
        

        ###   ML selection scheme
        #allQuads=hhhSelector.getBJetParisFGG_MHA(eTree,etaMax,pTMin,mlScoreTag='H3SIN61')
        ##print("\nRecoed bjets : ",idxs)
        #if allQuads['isValid']:
        #    dr = allQuads['bJetQuad']['r_HH']
        #    isSR = dr < 25.0
        #    if isSR:
        #        sumEntries.Fill('isROI_Q2', 1)
        #        sumWeights.Fill('isROI_Q2', wei)

        #    sumEntries.Fill('hasQuad', 1)
        #    #print("selection : ",allQuads['bJetQuad']['fgg_idxs'])
        #    th1Store["quadV2RecoMass_H1bb"].Fill(allQuads['bJetQuad']['m1'])
        #    th1Store["quadV2RecoMass_H2bb"].Fill(allQuads['bJetQuad']['m2'])
        #    th1Store["quadV2RecoMass_H1H2"].Fill(allQuads['bJetQuad']['m1'], allQuads['bJetQuad']['m2'])

        #    hasAllGenMatchedBJetsInQuad=True
        #    for i in range(4):
        #        ix=allQuads['bJetQuad']['fgg_idxs'][i]
        #        th1Store["quadV2RecoPt_"+str(i)].Fill( getattr(eTree,'jet_'+str(ix)+'_pt')  )
        #        th1Store["quadV2RecoEta_"+str(i)].Fill( getattr(eTree,'jet_'+str(ix)+'_eta')  )
        #        th1Store["quadV2Reco_b"+str(i)+"HFlav"].Fill( getattr(eTree,'jet_'+str(ix)+'_flavour')  )

        #        if ix not in idxs:
        #            hasAllGenMatchedBJetsInQuad=False
        #    if hasAllGenMatchedBJetsInQuad:
        #        sumEntries.Fill("hasAllGenMatchedBJetsInQuadV2",1)

        #    hasH1GenMatchedInQuad=False
        #    if allQuads['bJetQuad']['fgg_idxs'][0] in idxs[0:2] and allQuads['bJetQuad']['fgg_idxs'][1] in  idxs[0:2]:
        #        hasH1GenMatchedInQuad=True
        #    if allQuads['bJetQuad']['fgg_idxs'][0] in idxs[2:] and allQuads['bJetQuad']['fgg_idxs'][1] in  idxs[2:]:
        #        hasH1GenMatchedInQuad=True
        #    if hasH1GenMatchedInQuad:
        #        sumEntries.Fill("hasH1GenMatchedInQuadV2",1)
        #    
        #    hasH2GenMatchedInQuad=False
        #    if allQuads['bJetQuad']['fgg_idxs'][2] in idxs[0:2] and allQuads['bJetQuad']['fgg_idxs'][3] in  idxs[0:2]:
        #        hasH2GenMatchedInQuad=True
        #    if allQuads['bJetQuad']['fgg_idxs'][2] in idxs[2:] and allQuads['bJetQuad']['fgg_idxs'][3] in  idxs[2:]:
        #        hasH2GenMatchedInQuad=True
        #    if hasH2GenMatchedInQuad:
        #        sumEntries.Fill("hasH2GenMatchedInQuadV2",1)
        #    
        #    hasAllJetsMatchedProperly=True
        #    for i in range(4):
        #        if idxs[i]!=allQuads['bJetQuad']['fgg_idxs'][i]:
        #            hasAllJetsMatchedProperly=False
        #            break
        #    if hasAllJetsMatchedProperly:
        #        sumEntries.Fill("hasAllJetsMatchedOneToOneProperlyQuadV2",1)

        #    if hasH2GenMatchedInQuad and hasH1GenMatchedInQuad:
        #        sumEntries.Fill("hasEachHMatchedProperlyQuadV2",1)
        #else:
        #    sumEntries.Fill('noQuadV2', 1)
    
        if not isAllRecoed:
            continue
        sumEntries.Fill('isAllRecoed', 1)
        sumWeights.Fill('isAllRecoed', wei)
        
        isMerged=False
        for i in range(4):
            for j in range(4):
                if i <= j :
                    continue
                if idxs[i]==idxs[j]:
                    sumEntries.Fill('isMerged_'+str(i)+'_'+str(j), 1)
                    sumWeights.Fill('isMerged_'+str(i)+'_'+str(j), wei)
                    isMerged=True
        if isMerged:
            sumEntries.Fill('isMerged', 1)
            sumWeights.Fill('isMerged', wei)
            continue        
        for i in range(4):
            puJetIdMVA=getattr(eTree,'jet_'+str(idxs[i])+'_puJetIdMVA')
            pt=getattr(eTree,'jet_'+str(idxs[i])+'_pt')
            eta=getattr(eTree,'jet_'+str(idxs[i])+'_eta')

            if hhhUtil.puJetIDLoose(pt,puJetIdMVA):
                sumEntries.Fill('isPUJetIDLoosePass_'+str(i), 1)
                sumWeights.Fill('isPUJetIDLoosePass_'+str(i), wei)
        
            if hhhUtil.puJetIDTight(pt,puJetIdMVA):
                sumEntries.Fill('isPUJetIDTightPass_'+str(i), 1)
                sumWeights.Fill('isPUJetIDTightPass_'+str(i), wei)
        
            if hhhUtil.puJetIDMedium(pt,puJetIdMVA):
                sumEntries.Fill('isPUJetIDMediumPass_'+str(i), 1)
                sumWeights.Fill('isPUJetIDMediumPass_'+str(i), wei)
        
        isLeadsLoose=True
        for i in [0,2]:
            isLeadsLoose = isLeadsLoose and (getattr(eTree,'jet_'+str(idxs[i])+'_isLoose') > 0.25 )
        if isLeadsLoose:
            sumEntries.Fill('isLeadsLoose_'+str(i), 1)
            sumWeights.Fill('isLeadsLoose_'+str(i), wei)
                
        isLeadsTight=True
        for i in [0,2]:
            isLeadsTight = isLeadsTight and (getattr(eTree,'jet_'+str(idxs[i])+'_isTight') > 0.25 )
        if isLeadsTight:
            sumEntries.Fill('isLeadsTight_'+str(i), 1)
            sumWeights.Fill('isLeadsTight_'+str(i), wei)
        
        isSubLeadsLoose=True
        for i in [1,3]:
            isSubLeadsLoose = isSubLeadsLoose and (getattr(eTree,'jet_'+str(idxs[i])+'_isLoose') > 0.25 )
        if isSubLeadsLoose:
            sumEntries.Fill('isSubLeadsLoose_'+str(i), 1)
            sumWeights.Fill('isSubLeadsLoose_'+str(i), wei)
    
        isSubLeadsTight=True
        for i in [1,3]:
            isSubLeadsTight = isSubLeadsTight and (getattr(eTree,'jet_'+str(idxs[i])+'_isTight') > 0.25 )
        if isSubLeadsTight:
            sumEntries.Fill('isSubLeadsTight_'+str(i), 1)
            sumWeights.Fill('isSubLeadsTight_'+str(i), wei)
    
        if isLeadsTight and isSubLeadsTight:
            sumEntries.Fill('isAllTight_'+str(i), 1)
            sumWeights.Fill('isAllTight_'+str(i), wei)
    
        if isLeadsTight and isSubLeadsLoose:
            sumEntries.Fill('isLeadTightSLeadLoose_'+str(i), 1)
            sumWeights.Fill('isLeadTightSLeadLoose_'+str(i), wei)
           
    simFile.Close()           
    print("Closing file : ",fname)
dir_.cd()    

sumEntries.Write()
sumWeights.Write()

for i in range(4):
    th1Store["efficiencyPt_"+str(i)] =ROOT.TEfficiency( th1Store["recoGenPt_"+str(i)] , th1Store["allGenPt_"+str(i)] ) 
    th1Store["efficiencyPt_"+str(i)].SetName("efficiencyPt_"+str(i))
    
    th1Store["th1f_efficiencyPt_"+str(i)] = th1Store["recoGenPt_"+str(i)].Clone()
    for k in range(th1Store["th1f_efficiencyPt_"+str(i)].GetNbinsX()):
        th1Store["th1f_efficiencyPt_"+str(i)].SetBinContent( k , th1Store["efficiencyPt_"+str(i)].GetEfficiency(k) )
        err=0.5* (th1Store["efficiencyPt_"+str(i)].GetEfficiencyErrorLow(k)  +  th1Store["efficiencyPt_"+str(i)].GetEfficiencyErrorUp(k) )
        th1Store["th1f_efficiencyPt_"+str(i)].SetBinError( k , err)
    th1Store["th1f_efficiencyPt_"+str(i)].SetName("th1_efficiencyPt_"+str(i))

    th1Store["efficiencyEta_"+str(i)] = ROOT.TEfficiency(  th1Store["recoGenEta_"+str(i)] ,  th1Store["allGenEta_"+str(i)] )
    th1Store["efficiencyEta_"+str(i)].SetName("efficiencyEta_"+str(i))
    
    th1Store["th1f_efficiencyEta_"+str(i)] = th1Store["recoGenEta_"+str(i)].Clone()
    for k in range(th1Store["th1f_efficiencyEta_"+str(i)].GetNbinsX()):
        th1Store["th1f_efficiencyEta_"+str(i)].SetBinContent( k , th1Store["efficiencyEta_"+str(i)].GetEfficiency(k) )
        err=0.5* (th1Store["efficiencyEta_"+str(i)].GetEfficiencyErrorLow(k)  +  th1Store["efficiencyEta_"+str(i)].GetEfficiencyErrorUp(k) )
        th1Store["th1f_efficiencyEta_"+str(i)].SetBinError( k , err)
    th1Store["th1f_efficiencyEta_"+str(i)].SetName("th1_efficiencyEta_"+str(i))

    th1Store["cumulativeGenMatchDr_"+str(i)] = th1Store["genMatchDr_"+str(i)].GetCumulative()
    th1Store["cumulativeGenMatchDr_"+str(i)].SetName( "cumulativeGenMatchDr_"+str(i))
    th1Store["cumulativeGenMatchDr_"+str(i)].Scale( 1.0/th1Store["genMatchDr_"+str(i)].Integral()  )

    th1Store["cumulativeGenMatchWithFlavourDr_"+str(i)] = th1Store["genMatchWithFlavourDr_"+str(i)].GetCumulative()
    th1Store["cumulativeGenMatchWithFlavourDr_"+str(i)].SetName( "cumulativeGenMatchWithFlavourDr_"+str(i))
    th1Store["cumulativeGenMatchWithFlavourDr_"+str(i)].Scale( 1.0/th1Store["genMatchDr_"+str(i)].Integral()  )

for ky in th1Store:
    th1Store[ky].Write()
fout.Close()
print(" File written out  : ",foutName)

