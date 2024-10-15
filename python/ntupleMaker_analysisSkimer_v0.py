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

import os,sys

etaMax=2.5
pTMin=25.0


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
pTReweitingHistCatagories=getValueFromConfigs(cfgTxt,"pTReweitingHistCatagories",default="")

pTReweitingHistName=pTReweitingHistName.split(',')
pTReweitingHistCatagories=pTReweitingHistCatagories.split(',')

reweighter={}
for i,j in zip(pTReweitingHistCatagories,pTReweitingHistName):
    reweighter[i]=j
    print("set rew ",i,j)
print('reweighter : ',reweighter)
resetWeight=float(getValueFromConfigs(cfgTxt,"resetWeight",default=-1e5))


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
print("pTReweitingHistCatagories   :  ",    pTReweitingHistCatagories)
#print("pTReweitingValFile   :  ",     pTReweitingValFile)
#print("pTReweitingHistValName   :  ", pTReweitingHistValName)


scaler=mcScaler()
scalerVal=mcScaler()
if doPtReWeighting:
    scaler.setSFHistFromFile(pTReweitingFile,reweighter)

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
fout = ROOT.TFile(foutName, filemode)
dir_ = fout.mkdir('tagsDumper')
dir_ = dir_.mkdir('trees')

branches_skimmed=['M1jj','M2jj','CMS_hgg_mass','peakingMVA_v1','nonResonantMVA_v1','dZ','r_HH','weight','nonResonantMVA_v2']

ntuple={}
ntuple['all'] = ROOT.TNtuple(outTreeName+'_NOTAG', outTreeName, ':'.join(branches_skimmed))

histStore={}
tagList=['']
if processID=='DATA':
    tagList.append('_fitRange')
for tag in tagList:
    if filemode=='update':
        for cat in ['cat0','cat1']:
            ntuple[cat+tag] = fout.Get(outTreeName+'_'+cat)
    else:
        for cat in ['cat0','cat1']:
            ntuple[cat+tag] = ROOT.TNtuple(outTreeName+tag+'_'+cat, outTreeName+tag+'_'+cat, ':'.join(branches_skimmed))
            histStore[cat+tag]={}
            histStore[cat+tag]['m1jj'] = ROOT.TH1F(cat+'_'+tag+"_m1jj","H1bb mass",300,0.0,300.0 )
            histStore[cat+tag]['m2jj'] = ROOT.TH1F(cat+'_'+tag+"_m2jj","H1bb mass",300,0.0,300.0 )
            histStore[cat+tag]['m1m2'] = ROOT.TH2F(cat+'_'+tag+"_m1jj_m2jj","H1bb , H2bb mass",300,0.0,300.0,300,0.0,300. )
            histStore[cat+tag]['mvaNR_v0']=ROOT.TH1F(cat+'_'+tag+"mvaNR_v0","",44,-1.1,1.1)
            histStore[cat+tag]['mvaNR_v1']=ROOT.TH1F(cat+'_'+tag+"mvaNR_v1","",44,-1.1,1.1)

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

sumWeights=ROOT.TH1F("sumEvts","sumEvts",1,0.0,1.0)
sumWeights.SetCanExtend(ROOT.TH1.kAllAxes)

cutFlow=ROOT.TH1F("cutFlow","cutFlow",1,0.0,1.0)
cutFlow.SetCanExtend(ROOT.TH1.kAllAxes)

for fname in allFnames:
    print("Opening file : ",fname)
    simFile = ROOT.TFile(fname,'READ')
    eTree=simFile.Get(treeName)
    print(" NEntries = ", eTree.GetEntries() )
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
            
        dr = np.sqrt((eTree.M1jj-125.0)**2 + (eTree.M2jj - 125.0)**2)
        isMasked =  (dr < 25.0) and maskSignalMgg and processID=="DATA" and ( eTree.CMS_hgg_mass > 115.0) and (eTree.CMS_hgg_mass < 135.0)
        
        
        #if not insideTheEllipse( eTree.M1jj  , eTree.M2jj , 70.0  , 70.0 , 150.0 , 150.0 , 141.0):
        #if dr > 25.0:
        #    continue
        
        allQuads=hhhSelector.getBJetParisFGG(eTree,etaMax,pTMin)
        if not allQuads['isValid']:
            cutFlow.Fill(allQuads['fail'],1)
            continue

        LVStore = getLVStoreFromTreeAndQuad(eTree,allQuads['bJetQuad'])
        j1CosTheta,k1CosTheta,ggCostheta,drMin,drOther=getCosthetaVars(eTree,LVStore)
        
        tag=''
        if isMasked:
            tag='_fitRange'

        wei=eTree.weight
        if resetWeight > -1e-4:
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

            #scaleFactor=scalerVal.getSFForX(pT)
            #weiVal*=scaleFactor
        
        for ky in tofill:
            if ky in allBranches:
                tofill[ky]=getattr(eTree,ky)
            else:
                tofill[ky]=-1.111e3
        
        tofill['weight'] =  wei
        tofill['r_HH'] = allQuads['bJetQuad']['r_HH']
       
        mvaAll.Fill(tofill["peakingMVA_v1"])

        sumWeights.Fill('total', 1)
        sumWeights.Fill('total_wei', wei)
        m1m2.Fill(eTree.M1jj,eTree.M2jj)
        mvaPass.Fill(tofill["peakingMVA_v1"])

        for i in skimmedDataDict:
            skimmedDataDict[i]=tofill[i]
        
        cat='cat0'
        if tofill['r_HH'] < 40.0:
            cat = 'cat1'
        else:
            cutFlow.Fill('r_HH',1)
        sumWeights.Fill('total_'+cat, 1)
        sumWeights.Fill('total_wei'+cat, wei)
        

        #if tag!='':
        #ntuple[cat+tag].Fill(array('f', skimmedDataDict.values()))
        histStore[cat+tag]['m1jj'].Fill( quad['m1'] ,wei)
        histStore[cat+tag]['m2jj'].Fill( quad['m2'] ,wei)
        histStore[cat+tag]['m1m2'].Fill( quad['m1'] , quad['m2'] ,wei)
        histStore[cat+tag]['m1m2_evt'].Fill( quad['m1'] , quad['m2'])
        
        tofill['M1jj']=quad['m1'] 
        tofill['M2jj']=quad['m2'] 
        histStore[cat+tag]['mvaNR_v0'].Fill( tofill["nonResonantMVA_v0"] ,wei )
        histStore[cat+tag]['mvaNR_v1'].Fill( tofill["nonResonantMVA_v1"] ,wei )
        
        ntuple[cat+tag].Fill(array('f', skimmedDataDict.values()))
        ntuple['all'].Fill(array('f', skimmedDataDict.values()))

    simFile.Close()           
    print("Closing file : ",fname)
dir_.cd()    
sumWeights.Write()
cutFlow.Write()
m1m2.Write()
mvaAll.Write()
mvaPass.Write()
for cat in ntuple:
    ntuple[cat].Write()
for cat in histStore:
    for ky in histStore[cat]:
        histStore[cat][ky].SetName(cat+'_'+ky)
        histStore[cat][ky].Write()

fout.Close()
print(" File written out  : ",foutName)

