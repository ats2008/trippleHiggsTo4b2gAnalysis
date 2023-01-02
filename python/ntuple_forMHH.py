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
print("doSR     :       ", doSR)
print('etaMax   :       ', etaMax)
print('pTMin    :       ', pTMin)
print('drMax    :       ', drMax)

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

sumEntries=ROOT.TH1F("sumEvts","sumEvts",1,0.0,1.0)
sumEntries.SetCanExtend(ROOT.TH1.kAllAxes)
sumWeights=ROOT.TH1F("sumWeighs","sumWeighs",1,0.0,1.0)
sumWeights.SetCanExtend(ROOT.TH1.kAllAxes)

th1Store={}

for fname in allFnames:
    print("Opening file : ",fname)
    simFile = ROOT.TFile(fname,'READ')
    eTree=simFile.Get(treeName)
    print(" NEntries = ", eTree.GetEntries())
    if not eTree:
        eTree=simFile.GehFlavourt('tagsDumper/trees/Data_13TeV_TrippleHTag_0')
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
            print("   Doing i = ",i," / ",maxEvents_,
                  " n No HiggsCand : ", nNoHtoGammaGamma,
                  " n No 4BCands  : ", nNoHHto4B)
            print(" gg mass , h1MassVsh2Mass  : ",eTree.CMS_hgg_mass ,"( ",eTree.M1jj , eTree.M2jj,")" )
            
        

        rslt=getBestGetMatchesGlobalFGG(eTree,drMax,etaMax,pTMin)
        idxs=rslt['idxs']
        drMins=rslt['drMins']
        jetMass=rslt['jetMass']
        allJetMatchP4=rslt['p4s']
        deepJetScores=rslt['deepJetScores']
        hFlavour=rslt['hFlavour']
        isAllRecoed=True
        
        isRecoed={}
        for i in range(4):
            isRecoed['isRecoed_'+str(i)]=True
        for i in range(4):
            th1Store["allGenPt_"+str(i)].Fill( allBDaus[i].Pt()  )
            th1Store["allGenEta_"+str(i)].Fill( allBDaus[i].Eta() )
            th1Store["genMatchDr_"+str(i)].Fill( drMins[i])
            if idxs[i] < 0:
                isAllRecoed=False
                isRecoed['isRecoed_'+str(i)]=False
                th1Store["noGenMatchJetMass_"+str(i)].Fill( jetMass[i] )
                th1Store["noGenMatchDeepJet_"+str(i)].Fill( deepJetScores[i] )
                th1Store["noGenMatchHFlavour_"+str(i)].Fill( hFlavour[i] )
            else:
                sumEntries.Fill('isRecoed_'+str(i), 1)
                sumWeights.Fill('isRecoed_'+str(i), wei)
                th1Store["recoMatchDr_"+str(i)].Fill( drMins[i])
                th1Store["recoMatchJetMass_"+str(i)].Fill( jetMass[i] )
                th1Store["recoMatchDeepJet_"+str(i)].Fill( deepJetScores[i] )
                th1Store["recoMatchHFlavour_"+str(i)].Fill( hFlavour[i] )
                th1Store["recoGenPt_"+str(i)].Fill( allBDaus[i].Pt()  )
                th1Store["recoGenEta_"+str(i)].Fill( allBDaus[i].Eta() )
                th1Store["recoPt_"+str(i) ].Fill( getattr(eTree,'jet_'+str(idxs[i])+'_pt')  )
                th1Store["recoEta_"+str(i)].Fill( getattr(eTree,'jet_'+str(idxs[i])+'_eta') )
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
        
        allQuads=getBJetParisFGG(eTree,etaMax,pTMin)
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

                if ix not in idxs:
                    hasAllGenMatchedBJetsInQuad=False
            if hasAllGenMatchedBJetsInQuad:
                sumEntries.Fill("hasAllGenMatchedBJetsInQuad",1)

            hasH1GenMatchedInQuad=False
            if allQuads['bJetQuad']['fgg_idxs'][0] in idxs[0:2] and allQuads['bJetQuad']['fgg_idxs'][1] in  idxs[0:2]:
                hasH1GenMatchedInQuad=True
            if allQuads['bJetQuad']['fgg_idxs'][0] in idxs[2:] and allQuads['bJetQuad']['fgg_idxs'][1] in  idxs[2:]:
                hasH1GenMatchedInQuad=True
            if hasH1GenMatchedInQuad:
                sumEntries.Fill("hasH1GenMatchedInQuad",1)
            
            hasH2GenMatchedInQuad=False
            if allQuads['bJetQuad']['fgg_idxs'][2] in idxs[0:2] and allQuads['bJetQuad']['fgg_idxs'][3] in  idxs[0:2]:
                hasH2GenMatchedInQuad=True
            if allQuads['bJetQuad']['fgg_idxs'][2] in idxs[2:] and allQuads['bJetQuad']['fgg_idxs'][3] in  idxs[2:]:
                hasH2GenMatchedInQuad=True
            if hasH2GenMatchedInQuad:
                sumEntries.Fill("hasH2GenMatchedInQuad",1)
            
            hasAllJetsMatchedProperly=True
            for i in range(4):
                if idxs[i]!=allQuads['bJetQuad']['fgg_idxs'][i]:
                    hasAllJetsMatchedProperly=False
                    break
            if hasAllJetsMatchedProperly:
                sumEntries.Fill("hasAllJetsMatchedOneToOneProperly",1)

            if hasH2GenMatchedInQuad and hasH1GenMatchedInQuad:
                sumEntries.Fill("hasEachHMatchedProperly",1)

                    
        else:
            sumEntries.Fill('noQuad', 1)
 
    
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

            if puJetIDLoose(pt,puJetIdMVA):
                sumEntries.Fill('isPUJetIDLoosePass_'+str(i), 1)
                sumWeights.Fill('isPUJetIDLoosePass_'+str(i), wei)
        
            if puJetIDTight(pt,puJetIdMVA):
                sumEntries.Fill('isPUJetIDTightPass_'+str(i), 1)
                sumWeights.Fill('isPUJetIDTightPass_'+str(i), wei)
        
            if puJetIDMedium(pt,puJetIdMVA):
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
    
    th1Store["efficiencyEta_"+str(i)] = ROOT.TEfficiency(  th1Store["recoGenEta_"+str(i)] ,  th1Store["allGenEta_"+str(i)] )
    th1Store["efficiencyEta_"+str(i)].SetName("efficiencyEta_"+str(i))

    th1Store["cumulativeGenMatchDr_"+str(i)] = th1Store["genMatchDr_"+str(i)].GetCumulative()
    th1Store["cumulativeGenMatchDr_"+str(i)].SetName( "cumulativeGenMatchDr_"+str(i))
    th1Store["cumulativeGenMatchDr_"+str(i)].Scale( 1.0/th1Store["genMatchDr_"+str(i)].Integral()  )

for ky in th1Store:
    th1Store[ky].Write()
fout.Close()
print(" File written out  : ",foutName)

