from __future__ import print_function
from collections import OrderedDict
import ROOT 
import numpy as np
from trippleHiggsUtils import *
from Util  import *
from branches import *
from array import array
from trippleHiggsSelector import *

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

for ky in ['j1','j2','k1','k2','g1','g1','g2','g2','j1','j2','k1','k2','Hgg','Hgg','H1bb','H2bb','HHH','H1','H2','H3']:
    for kyEx in ['px','py','pz','E','mass']:
        branchesTTH.append(ky+kyEx)
branches=np.unique(branchesTTH)


isMC = True
isMC = False
filemode="RECREATE"
fout = ROOT.TFile(foutName, filemode)
if filemode=='update':
    ntuple = fout.Get('tree')
else:
    ntuple = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
tofill = OrderedDict(zip(branches, [np.nan]*len(branches)))

kyList = [ ky for ky in tofill ]

print("len(branches) : " , len(branches))

sumWeights=ROOT.TH1F("sumEvts","sumEvts",1,0.0,1.0)
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
        if processID=="MC":
            #wei*=weightScale
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
            pT=eTree.diphoton_pt
            scaleFactor=scaler.getSFForX(pT)
            wei*=scaleFactor
            scaleFactor=scalerVal.getSFForX(pT)
            weiVal*=scaleFactor
        
        isBTaggedAllLoose = (eTree.h1LeadingJet_DeepFlavour > 0.0490 ) and ( eTree.h1SubleadingJet_DeepFlavour > 0.0490 ) and ( eTree.h2LeadingJet_DeepFlavour > 0.0490 ) and ( eTree.h2SubleadingJet_DeepFlavour > 0.0490 )
        #if not isBTaggedAllLoose:
        #    continue

        ## Control and signal region
        dr = np.sqrt((eTree.M1jj-125.0)**2 + (eTree.M2jj - 125.0)**2)

        if dr>50.0:
            continue
        
        #isMasked =  (dr < 25.0) and maskSignalMgg and processID=="DATA" and ( eTree.CMS_hgg_mass > 115.0) and (eTree.CMS_hgg_mass < 135.0)
        #if isMasked:
        #    continue

        for ky in tofill:
            if ky in allBranches:
                tofill[ky]=getattr(eTree,ky)
            else:
                tofill[ky]=-1.111e3

        #if eTree.nonResonantMVA_v1 <0.5:
        #    continue

        LVStore = getLVStore(eTree)
        j1CosTheta,k1CosTheta,ggCostheta,drMin,drOther=getCosthetaVars(eTree,LVStore)
        for ky in ['j1','j2','k1','k2','g1','g1','g2','g2','j1','j2','k1','k2','Hgg','Hgg','H1bb','H2bb','HHH','H1','H2','H3']:
            tofill[ky+'px']=LVStore[ky+'LV'].Px()
            tofill[ky+'py']=LVStore[ky+'LV'].Py()
            tofill[ky+'pz']=LVStore[ky+'LV'].Pz()
            tofill[ky+'E']=LVStore[ky+'LV'].E()
            tofill[ky+'mass']=LVStore[ky+'LV'].M()
        
        sumWeights.Fill("sumWeights",wei)

        tofill['weight'] =  wei
        tofill['r_HH'] =  dr
        tofill['hhh_pT'] =  LVStore['HHHLV'].Pt() 
        tofill['hhhCosThetaH1'] = eTree.absCosThetaStar_CS 
        tofill["hh4CosThetaLeadJet"] = eTree.absCosTheta_bb
        tofill["h1bbCosThetaLeadJet"]= abs(j1CosTheta)
        tofill["h2bbCosThetaLeadJet"]= abs(k1CosTheta)
        tofill["h2bbCosThetaLeadJet"]= abs(k1CosTheta)

        tofill["PhoJetMinDr"]= drMin
        tofill["PhoJetOtherDr"]= drOther
        tofill['pTleadG_overMgg'] =  eTree.leadingPhoton_pt / eTree.CMS_hgg_mass 
        tofill['pTh1leadJ_overMh1'] =  eTree.h1LeadingJet_pt / eTree.M1jj 
        tofill['pTh2leadJ_overMh2'] =  eTree.h2LeadingJet_pt / eTree.M2jj 

        tofill['pTsubleadG_overMgg'] =  eTree.subleadingPhoton_pt / eTree.CMS_hgg_mass 
        tofill['pTh1subleadJ_overMh1'] =  eTree.h1SubleadingJet_pt / eTree.M1jj 
        tofill['pTh2subleadJ_overMh2'] =  eTree.h2SubleadingJet_pt / eTree.M2jj 


        tofill["absCosThetaH4bHgg"] =  np.cos( (LVStore['H1bbLV'] + LVStore['H2bbLV']).Angle(LVStore['HggLV'].Vect()))
        vals=[
            LVStore['j1LV'].Angle( LVStore['j2LV'].Vect()),
            LVStore['j1LV'].Angle( LVStore['k1LV'].Vect()),
            LVStore['j1LV'].Angle( LVStore['k2LV'].Vect()),
        ]
        vals=abs(np.cos(vals))
        tofill["LeadJetAbsCosThetaMax"]  = max(vals)
        tofill["LeadJetAbsCosThetaMin"]  = min(vals)
       
        vals=[
            LVStore['j1LV'].DeltaR( LVStore['j2LV']),
            LVStore['j1LV'].DeltaR( LVStore['k1LV']),
            LVStore['j1LV'].DeltaR( LVStore['k2LV']),
        ]
        tofill["LeadJetDrMax"]  = max(vals)
        tofill["LeadJetDrMin"]  = min(vals)
         
        vals=[
            LVStore['j1LV'].Angle( LVStore['k1LV'].Vect()),
            LVStore['j1LV'].Angle( LVStore['k2LV'].Vect()),
            LVStore['j2LV'].Angle( LVStore['k1LV'].Vect()),
            LVStore['j2LV'].Angle( LVStore['k2LV'].Vect()),
        ]
        vals=abs(np.cos(vals))
        tofill["H1H2JetAbsCosThetaMax"]  = max(vals)
        tofill["H1H2JetAbsCosThetaMin"]  = min(vals)         
        
        vals=[
            LVStore['j1LV'].DeltaR( LVStore['k1LV']),
            LVStore['j1LV'].DeltaR( LVStore['k2LV']),
            LVStore['j2LV'].DeltaR( LVStore['k1LV']),
            LVStore['j2LV'].DeltaR( LVStore['k2LV']),
        ]
        tofill["H1H2JetDrMax"]  = max(vals)
        tofill["H1H2JetDrMin"]  = min(vals)
        tofill["pT_4b"]  = ( LVStore['H1bbLV'] + LVStore['H2bbLV']).Pt()
        
        tofill["scalarPtSumHHH"]  = LVStore["H1LV"].Pt()+LVStore["H2LV"].Pt()+LVStore["H3LV"].Pt()
        tofill["scalarPtSum4b"]   = LVStore["j1LV"].Pt() + LVStore["j2LV"].Pt() + LVStore["k1LV"].Pt() + LVStore["k2LV"].Pt()
        tofill["scalarPtSum4b2g"] = tofill["scalarPtSum4b"]  + LVStore["g1LV"].Pt() + LVStore["g2LV"].Pt()

        tofill["H1To4bAbsCosTheta"]  = np.cos( (LVStore['H1bbLV']+LVStore['H2bbLV']).Angle( LVStore['H1LV'].Vect())  )
        
        #print(tofill["H1To4bAbsCosTheta"])
        #LVStore['H1bbLV'].Print()
        #LVStore['H2bbLV'].Print() 
        #LVStore['H1LV'].Vect().Print()
        #print((LVStore['H1bbLV']+LVStore['H2bbLV']).Angle( LVStore['H1LV'].Vect()) ,  np.cos( (LVStore['H1bbLV']+LVStore['H2bbLV']).Angle( LVStore['H1LV'].Vect())  ))
        
        tofill["H1bbToH2bbAbsCosTheta"]  = np.cos(  LVStore['H1bbLV'].Angle( LVStore['H2bbLV'].Vect())  )         
        ntuple.Fill(array('f', tofill.values()))
        
        #print(len(tofill) , " / ",len(branches))
        #print("tofill not in key")
        #for ky in tofill:
        #    if ky not in kyList:
        #        print(ky)

        #print("key not in tofill")
        #for ky in kyList:
        #    if ky not in tofill:
        #        print(ky)
        #print("branch not in tofill")
        #for ky in branches:
        #    if ky not in tofill:
        #        print(ky)


    simFile.Close()           
    print("Closing file : ",fname)
fout.cd()    
sumWeights.Write()
ntuple.Write()
fout.Close()
print(" File written out  : ",foutName)

