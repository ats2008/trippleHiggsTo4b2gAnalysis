from __future__ import print_function
from collections import OrderedDict
import ROOT 
import numpy as np
from trippleHiggsUtils import *
from Util  import *
from branches import *
from array import array
from trippleHiggsSelector import *
from TMVA_Model import *

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
if len(sys.argv) >2:
    maxEvtsSuperSeeder=int(sys.argv[2])
if len(sys.argv) >3:
    etaMax=float(sys.argv[3])
if len(sys.argv) >4:
    ext=sys.argv[4]


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

MVAWeightFile=getValueFromConfigs(cfgTxt,"MVAWeightFile",default="")
MVABranches=getListOfStringsFromConfigs(cfgTxt,"#MVAVARLIST_BEG","#MVAVARLIST_END")
MVASpecBranches=getListOfStringsFromConfigs(cfgTxt,"#SPECTATORLIST_BEG","#SPECTATORLIST_END")

nonResonantMVA=TMVAModel()
#print(MVAWeightFile, MVABranches , MVASpecBranches)
nonResonantMVA.setupTMVAModel("aMVA", MVAWeightFile , MVABranches , MVASpecBranches )

tthMVAWeightFile=getValueFromConfigs(cfgTxt,"tthMVAWeightFile",default="")
tthMVABranches=getListOfStringsFromConfigs(cfgTxt,"#tthMVAVARLIST_BEG","#tthMVAVARLIST_END")
tthMVASpecBranches=getListOfStringsFromConfigs(cfgTxt,"#tthSPECTATORLIST_BEG","#tthSPECTATORLIST_END")

tthMVA=TMVAModel()
#print(tthMVAWeightFile, tthMVABranches , tthMVASpecBranches)
tthMVA.setupTMVAModel("aMVA", tthMVAWeightFile , tthMVABranches , tthMVASpecBranches )

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
sumWeights=ROOT.TH1F("sumWeighs","sumWeighs",1,0.0,1.0)
sumWeights.SetCanExtend(ROOT.TH1.kAllAxes)

th1Store={}
th1Store["allGenDR_H1bb"]= ROOT.TH1F("allGenDR_H1bb","",60,0.0,6.0)
th1Store["allGenDR_H2bb"]= ROOT.TH1F("allGenDR_H2bb","",60,0.0,6.0)
th1Store["allGenPt_H1bb"]= ROOT.TH1F("allGenPt_H1bb","",600,0.0,600.0)
th1Store["allGenPt_H2bb"]= ROOT.TH1F("allGenPt_H2bb","",600,0.0,600.0)
th1Store["allGenEta_H1bb"]= ROOT.TH1F("allGenEta_H1bb","",100,-5.0,5.0)
th1Store["allGenEta_H2bb"]= ROOT.TH1F("allGenEta_H2bb","",100,-5.0,5.0)

th1Store['acceptance2d']=ROOT.TH1F("acceptance2d","acceptance2d",1,0.0,1.0)
th1Store['acceptance2d'].SetCanExtend(ROOT.TH1.kAllAxes)

for i in range(4):
    th1Store["allGenPt_"+str(i)]= ROOT.TH1F("allGenPt_"+str(i),"",200,0.0,500.0)
    th1Store["recoGenPt_"+str(i)]= ROOT.TH1F("recoGenPt_"+str(i),"",200,0.0,500.0)
    th1Store["recoPt_"+str(i)]= ROOT.TH1F("recoPt_"+str(i),"",200,0.0,500.0)
for i in range(4):
    th1Store["allGenEta_"+str(i)]= ROOT.TH1F("allGenEta_"+str(i),"",100,-5.0,5.0)
    th1Store["recoGenEta_"+str(i)]= ROOT.TH1F("recoGenEta_"+str(i),"",100,-5.0,5.0)
    th1Store["recoEta_"+str(i)]= ROOT.TH1F("recoEta_"+str(i),"",100,-5.0,5.0)

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
        wei=1.0 #eTree.weight
        sumEntries.Fill('total', 1)
        sumWeights.Fill('total', wei)
        if(i%500==0):
            print("   Doing i = ",i," / ",maxEvents_,
                  " n No HiggsCand : ", nNoHtoGammaGamma,
                  " n No 4BCands  : ", nNoHHto4B)
            
        
        
        ## Control and signal region
        dr = 1.0 # np.sqrt((eTree.M1jj-125.0)**2 + (eTree.M2jj - 125.0)**2)
        allGammaDaus=getHiggsDauP4s(eTree,22)
        allBDaus=getHiggsDauP4s(eTree,5)
        isSR = dr < 25.0
        
        
        if doSR and isSR:
            continue
        #if( (allGammaDaus[0]+allGammaDaus[1]).Pt() < 100 ): continue;

        sumEntries.Fill('isROI', 1)
        sumWeights.Fill('isROI', wei)
    
        isMasked =  (dr < 25.0) and maskSignalMgg and processID=="DATA" and ( eTree.CMS_hgg_mass > 115.0) and (eTree.CMS_hgg_mass < 135.0)
        

        idxs=getBestGetMatchesGlobal(eTree)
        isAllRecoed=True
        
        #print(idxs)       
        isRecoed={}
        for i in range(4):
            isRecoed['isRecoed_'+str(i)]=True
        for i in range(4):
            th1Store["allGenPt_"+str(i)].Fill( allBDaus[i].Pt()  )
            th1Store["allGenEta_"+str(i)].Fill( allBDaus[i].Eta() )
            if idxs[i] < 0:
                isAllRecoed=False
                isRecoed['isRecoed_'+str(i)]=False
            else:
                sumEntries.Fill('isRecoed_'+str(i), 1)
                sumWeights.Fill('isRecoed_'+str(i), wei)
                th1Store["recoGenPt_"+str(i)].Fill( allBDaus[i].Pt()  )
                th1Store["recoGenEta_"+str(i)].Fill( allBDaus[i].Eta() )
                th1Store["recoPt_"+str(i) ].Fill( eTree.jets_pt[i])  )
                th1Store["recoEta_"+str(i)].Fill( eTree.jets_eta[i]) )
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
        
        sumEntries.Fill(str(nout)+'_JetsNotIsRecoed',1)
        sumEntries.Fill(str(nleadOut)+'_LeadNotIsRecoed',1)
        sumEntries.Fill(str(nSubLeadOut)+'_SubLeadNotIsRecoed',1)
        if h1Out:
            sumEntries.Fill('H1NotIsRecoed', 1)
        if h2Out:
            sumEntries.Fill('H2NotIsRecoed', 1)
            
        if isRecoed['isRecoed_0'] and isRecoed['isRecoed_1'] and isRecoed['isRecoed_2'] and isRecoed['isRecoed_3'] :
            sumEntries.Fill('allBJetsIsRecoed', 1)
        else:
            pass

    
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
        #for i in range(4):
        #    puJetIdMVA=getattr(eTree,'jet_'+str(idxs[i])+'_puJetIdMVA')
        #    pt=getattr(eTree,'jet_'+str(idxs[i])+'_pt')
        #    eta=getattr(eTree,'jet_'+str(idxs[i])+'_eta')

        #    if puJetIDLoose(pt,puJetIdMVA):
        #        sumEntries.Fill('isPUJetIDLoosePass_'+str(i), 1)
        #        sumWeights.Fill('isPUJetIDLoosePass_'+str(i), wei)
        #
        #    if puJetIDTight(pt,puJetIdMVA):
        #        sumEntries.Fill('isPUJetIDTightPass_'+str(i), 1)
        #        sumWeights.Fill('isPUJetIDTightPass_'+str(i), wei)
        #
        #    if puJetIDMedium(pt,puJetIdMVA):
        #        sumEntries.Fill('isPUJetIDMediumPass_'+str(i), 1)
        #        sumWeights.Fill('isPUJetIDMediumPass_'+str(i), wei)
        #
        #isLeadsLoose=True
        # for i in [0,2]:
        #     isLeadsLoose = isLeadsLoose and (getattr(eTree,'jet_'+str(idxs[i])+'_isLoose') > 0.25 )
        # if isLeadsLoose:
        #     sumEntries.Fill('isLeadsLoose_'+str(i), 1)
        #     sumWeights.Fill('isLeadsLoose_'+str(i), wei)
        #         
        # isLeadsTight=True
        # for i in [0,2]:
        #     isLeadsTight = isLeadsTight and (getattr(eTree,'jet_'+str(idxs[i])+'_isTight') > 0.25 )
        # if isLeadsTight:
        #     sumEntries.Fill('isLeadsTight_'+str(i), 1)
        #     sumWeights.Fill('isLeadsTight_'+str(i), wei)
        # 
        # isSubLeadsLoose=True
        # for i in [1,3]:
        #     isSubLeadsLoose = isSubLeadsLoose and (getattr(eTree,'jet_'+str(idxs[i])+'_isLoose') > 0.25 )
        # if isSubLeadsLoose:
        #     sumEntries.Fill('isSubLeadsLoose_'+str(i), 1)
        #     sumWeights.Fill('isSubLeadsLoose_'+str(i), wei)
    
        # isSubLeadsTight=True
        # for i in [1,3]:
        #     isSubLeadsTight = isSubLeadsTight and (getattr(eTree,'jet_'+str(idxs[i])+'_isTight') > 0.25 )
        # if isSubLeadsTight:
        #     sumEntries.Fill('isSubLeadsTight_'+str(i), 1)
        #     sumWeights.Fill('isSubLeadsTight_'+str(i), wei)
    
        # if isLeadsTight and isSubLeadsTight:
        #     sumEntries.Fill('isAllTight_'+str(i), 1)
        #     sumWeights.Fill('isAllTight_'+str(i), wei)
    
        # if isLeadsTight and isSubLeadsLoose:
        #     sumEntries.Fill('isLeadTightSLeadLoose_'+str(i), 1)
        #     sumWeights.Fill('isLeadTightSLeadLoose_'+str(i), wei)
    simFile.Close()           
    print("Closing file : ",fname)
dir_.cd()    

sumEntries.Write()
sumWeights.Write()

for i in range(4):
    th1Store["efficiencyPt_"+str(i)] =  th1Store["recoGenPt_"+str(i)].Clone()
    th1Store["efficiencyPt_"+str(i)].SetName("efficiencyPt_"+str(i))
    th1Store["efficiencyPt_"+str(i)].Divide( th1Store["allGenPt_"+str(i)]  )
    th1Store["efficiencyEta_"+str(i)] =  th1Store["recoGenEta_"+str(i)].Clone()
    th1Store["efficiencyEta_"+str(i)].SetName("efficiencyEta_"+str(i))
    th1Store["efficiencyEta_"+str(i)].Divide( th1Store["allGenEta_"+str(i)]  )


for ky in th1Store:
    th1Store[ky].Write()
fout.Close()
print(" File written out  : ",foutName)

