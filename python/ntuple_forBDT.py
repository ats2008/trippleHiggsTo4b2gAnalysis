from __future__ import print_function
import ROOT 
from collections import OrderedDict
import numpy as np
from trippleHiggsUtils import *
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
 
allFnames=getListOfStringsFromConfigs(cfgTxt,"#FNAMES_BEG","#FNAMES_END")
foutName=getValueFromConfigs(cfgTxt,"OutpuFileName","fggHists.root")
processID=getValueFromConfigs(cfgTxt,"processID",default="DATA")
treeName=getValueFromConfigs(cfgTxt,"treeName",default="tagsDumper/trees/Data_13TeV_TrippleHTag_0")
outTreeName=getValueFromConfigs(cfgTxt,"outTreeName",default="Data_13TeV_TrippleHTag_0")
weightScale=float(getValueFromConfigs(cfgTxt,"WeightScale",default="1.0"))
resetWeight=float(getValueFromConfigs(cfgTxt,"resetWeight",default=-1e5))
doPtReWeighting=getBoolFromConfigs(cfgTxt,"doPtReWeighting",default=False)
pTReweitingFile=getValueFromConfigs(cfgTxt,"pTReweitingFile",default="")
pTReweitingHistName=getValueFromConfigs(cfgTxt,"pTReweitingHistName",default="")
mlScoreTag=getValueFromConfigs(cfgTxt,"mlScoreTag",default="")

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
dir_ = fout.mkdir('trees')

branches=[
    'quadjet_0_deepJetScore', 'quadjet_0_mlScore', 'quadjet_0_mlScoreY1s', 'quadjet_0_mlScoreY2s' ,
    'quadjet_1_deepJetScore', 'quadjet_1_mlScore', 'quadjet_1_mlScoreY1s', 'quadjet_1_mlScoreY2s' ,
    'quadjet_2_deepJetScore', 'quadjet_2_mlScore', 'quadjet_2_mlScoreY1s', 'quadjet_2_mlScoreY2s' ,
    'quadjet_3_deepJetScore', 'quadjet_3_mlScore', 'quadjet_3_mlScoreY1s', 'quadjet_3_mlScoreY2s' ,
    'r_HH', 'weight', 'hhh_pT', 'hhhCosThetaH1', 'hh4CosThetaLeadJet', 'sumScore_4j' ,
    'h1bbCosThetaLeadJet', 'h2bbCosThetaLeadJet', 'h2bbCosThetaLeadJet', 'sumScore_3j',
    'PhoJetMinDr', 'PhoJetOtherDr', 'h1bb_mass', 'h2bb_mass', 'h1bb_pt', 'r_HH','D_HH',
    'h2bb_pt', 'h1bb_eta', 'h2bb_eta', 'h1bb_phi', 'h2bb_phi', 'pTleadG_overMgg', 
    'pTh1leadJ_overMh1', 'pTh2leadJ_overMh2', 'pTsubleadG_overMgg', 'pTh1subleadJ_overMh1', 
    'pTh2subleadJ_overMh2', 'absCosThetaH4bHgg', 'LeadJetAbsCosThetaMax', 'LeadJetAbsCosThetaMin', 
    'LeadJetDrMaxWithOtherJets', 'LeadJetDrMinWithOtherJets', 'H1H2JetAbsCosThetaMax', 
    'H1H2JetAbsCosThetaMin', 'H1H2JetDrMax', 'H1H2JetDrMin', 'pT_4b', 'scalarPtSumHHH', 
    'scalarPtSum4b', 'scalarPtSum4b2g', 'HggTo4bAbsCosTheta', 'H1bbToH2bbAbsCosTheta', 'ttH_MET', 
    'diphotonCandidatePtOverdiHiggsM', 'dije1CandidatePtOverdiHiggsM', 'dije2CandidatePtOverdiHiggsM', 
    'leadingPhotonSigOverE', 'subleadingPhotonSigOverE', 'sigmaMOverM', 
    'sigmaM1OverMJets', 'sigmaM2OverMJets', 'CMS_hgg_mass', 'weight',
    'customLeadingPhotonIDMVA','customSubLeadingPhotonIDMVA'
]

ntuple={}
branches=np.unique(branches)
ntuple['eventTree']  = ROOT.TNtuple(outTreeName, outTreeName, ':'.join(branches))

tofill = OrderedDict(zip(branches, [np.nan]*len(branches)))
outputDataDict=OrderedDict(zip(branches, [np.nan]*len(branches))) 

kyList = [ ky for ky in tofill ]

print("len(branches) : " , len(branches))

sumEntries=ROOT.TH1F("sumEvts","sumEvts",1,0.0,1.0)
sumEntries.SetCanExtend(ROOT.TH1.kAllAxes)
sumWeights=ROOT.TH1F("sumWeighs","sumWeighs",1,0.0,1.0)
sumWeights.SetCanExtend(ROOT.TH1.kAllAxes)

th1Store={}


beg=datetime.datetime.now()

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
        
        wei=eTree.weight
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
         

          
        for ky in tofill:
            if ky in allBranches:
                tofill[ky]=getattr(eTree,ky)
            else:
                tofill[ky]=0.0
        
        
        ##   JET PRE-SELECTION
        nVld=0
        jetMask=[]
        for i in range(8):
            jetMask.append(True)
            tofill['jet_'+str(i)+'_isValid']=getattr(eTree,'jet_'+str(i)+'_isValid')
            if abs(getattr(eTree,'jet_'+str(i)+'_eta') ) > etaMax:
                tofill['jet_'+str(i)+'_isValid']=0
            if abs(getattr(eTree,'jet_'+str(i)+'_pt') )  < pTMin:
                tofill['jet_'+str(i)+'_isValid']=0
            #if deltaR(getattr(eTree,'jet_'+str(i)+'_eta') , getattr(eTree,'jet_'+str(i)+'_phi') ,eTree.leadingPhoton_eta,eTree.leadingPhoton_phi ) < 0.4 :
            #    tofill['jet_'+str(i)+'_isValid']=0
            #if deltaR(getattr(eTree,'jet_'+str(i)+'_eta') , getattr(eTree,'jet_'+str(i)+'_phi') ,eTree.subleadingPhoton_eta,eTree.subleadingPhoton_phi ) < 0.4 :
            #    tofill['jet_'+str(i)+'_isValid']=0
            if tofill['jet_'+str(i)+'_isValid'] > 0.5:
                nVld+=1
            else:
                jetMask[-1]=False

        jetMask=np.array(jetMask)      
        if nVld < 4:
            continue

        allQuads=hhhSelector.getBJetParis_wrapper(eTree, mask= jetMask , methord='mha',mlScoreTag='H3SIN61',threshold=-1e3,doMisclassificationCorrection = False)
 #   def getBJetParisFGG_MHA( eTree, mask=None, mlScoreTag='', threshold=-1e3, doMisclassificationCorrection=True,returnOnlyJets=False):
        if not allQuads['isValid']:
            print("anity check failed after the quad finding")
            exit(1)
            continue
        sumScore_3j=0
        sumScore_4j=0
        verbose=False
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
        LVStore = getLVStoreFromTreeAndQuad(eTree,quad)
        j1CosTheta,k1CosTheta,ggCostheta,drMin,drOther=getCosthetaVars(eTree,LVStore)
        
        tofill['r_HH'] = quad['r_HH']
        tofill['D_HH'] = quad['D_HH']
        tofill['weight'] =  wei
        tofill['hhh_pT'] =  LVStore['HHHLV'].Pt() 
        tofill['hhhCosThetaH1'] = eTree.absCosThetaStar_CS 
        tofill["hh4CosThetaLeadJet"] = eTree.absCosTheta_bb
        
        tofill["h1bbCosThetaLeadJet"]= abs(j1CosTheta)
        tofill["h2bbCosThetaLeadJet"]= abs(k1CosTheta)
        tofill["h2bbCosThetaLeadJet"]= abs(k1CosTheta)
        
        tofill["PhoJetMinDr"]= drMin
        tofill["PhoJetOtherDr"]= drOther
        
        tofill["h1bb_mass"] = LVStore['H1bbLV'].M()
        tofill["h2bb_mass"] = LVStore['H2bbLV'].M()

        tofill["h1bb_pt"] = LVStore['H1bbLV'].Pt()
        tofill["h2bb_pt"] = LVStore['H2bbLV'].Pt()

        tofill["h1bb_eta"] = LVStore['H1bbLV'].Eta()
        tofill["h2bb_eta"] = LVStore['H2bbLV'].Eta()

        tofill["h1bb_phi"] = LVStore['H1bbLV'].Phi()
        tofill["h2bb_phi"] = LVStore['H2bbLV'].Phi()

        tofill['pTleadG_overMgg'] =  eTree.leadingPhoton_pt  / eTree.CMS_hgg_mass 
        tofill['pTh1leadJ_overMh1'] =  eTree.h1LeadingJet_pt / LVStore['H1bbLV'].M()
        tofill['pTh2leadJ_overMh2'] =  eTree.h2LeadingJet_pt / LVStore['H2bbLV'].M()

        tofill['pTsubleadG_overMgg']   =  eTree.subleadingPhoton_pt / eTree.CMS_hgg_mass 
        tofill['pTh1subleadJ_overMh1'] =  eTree.h1SubleadingJet_pt / LVStore['H1bbLV'].M()
        tofill['pTh2subleadJ_overMh2'] =  eTree.h2SubleadingJet_pt / LVStore['H1bbLV'].M() 

        tofill["absCosThetaH4bHgg"]    =  np.cos( (LVStore['H1bbLV'] + LVStore['H2bbLV']).Angle(LVStore['HggLV'].Vect()))
        
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
        tofill["LeadJetDrMaxWithOtherJets"]  = max(vals)
        tofill["LeadJetDrMinWithOtherJets"]  = min(vals)
         
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

        tofill['dije1CandidatePtOverdiHiggsM'] =   LVStore['H1bbLV'].Pt() / 125.0
        tofill['dije2CandidatePtOverdiHiggsM'] =   LVStore['H2bbLV'].Pt() / 125.0
        tofill["H1H2JetDrMax"]  = max(vals)
        tofill["H1H2JetDrMin"]  = min(vals)
        tofill["pT_4b"]  = ( LVStore['H1bbLV'] + LVStore['H2bbLV'] ).Pt()
        tofill["scalarPtSumHHH"]  = LVStore["H1LV"].Pt()+LVStore["H2LV"].Pt()+LVStore["H3LV"].Pt()
        tofill["scalarPtSum4b"]   = LVStore["j1LV"].Pt() + LVStore["j2LV"].Pt() + LVStore["k1LV"].Pt() + LVStore["k2LV"].Pt()
        tofill["scalarPtSum4b2g"] = tofill["scalarPtSum4b"]  + LVStore["g1LV"].Pt() + LVStore["g2LV"].Pt()
        tofill["HggTo4bAbsCosTheta"]  = np.cos( (LVStore['H1bbLV']+LVStore['H2bbLV']).Angle( LVStore['H1LV'].Vect())  )
        tofill["H1bbToH2bbAbsCosTheta"]  = np.cos(  LVStore['H1bbLV'].Angle( LVStore['H2bbLV'].Vect())  )         
        j1_res = getattr(eTree,'jet_'+str(quad['idxs'][0])+'_bJetRegRes') ; 
        j2_res = getattr(eTree,'jet_'+str(quad['idxs'][1])+'_bJetRegRes') ; 
        k1_res = getattr(eTree,'jet_'+str(quad['idxs'][2])+'_bJetRegRes') ; 
        k2_res = getattr(eTree,'jet_'+str(quad['idxs'][3])+'_bJetRegRes') ; 
        tofill['sigmaM1OverMJets'] = getSigmaMOverM( LVStore['j1LV'],j1_res , LVStore['j2LV'] , j2_res )
        tofill['sigmaM2OverMJets'] = getSigmaMOverM( LVStore['k1LV'],k1_res , LVStore['k2LV'] , k2_res )

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

for ky in th1Store:
    th1Store[ky].Write()
for ky in ntuple:
    ntuple[ky].Write()

fout.Close()

print(" File written out  : ",foutName)
