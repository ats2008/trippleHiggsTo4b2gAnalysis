from __future__ import print_function
import sys,hashlib,json
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
from scaleFactorUtil import mcScaler
import classifierUtil as clsUtl
import scaleFactorUtil
from TMVA_Model import *

from correctionlib import _core


import datetime,pickle
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
allFnames = [ fls.replace('"','') for fls in allFnames ]
foutName  =utl.getValueFromConfigs(cfgTxt,"OutpuFileName","fggHists.root")
processID =utl.getValueFromConfigs(cfgTxt,"processID",default="DATA")
treeName  =utl.getValueFromConfigs(cfgTxt,"treeName",default="tagsDumper/trees/Data_13TeV_TrippleHTag_0")
outTreeName=utl.getValueFromConfigs(cfgTxt,"outTreeName",default="Data_13TeV_TrippleHTag_0")

weightScale=float(utl.getValueFromConfigs(cfgTxt,"WeightScale",default="1.0"))
resetWeight=float(utl.getValueFromConfigs(cfgTxt,"resetWeight",default=-1e5))
doWeight=utl.getBoolFromConfigs(cfgTxt,"doWeight",default=False)
doSR=utl.getBoolFromConfigs(cfgTxt,"doSR",default=False)

doBDTReWeighting=utl.getBoolFromConfigs(cfgTxt,"doBDTReWeighting",default=False)
bdtReweitingFile=utl.getValueFromConfigs(cfgTxt,"bdtReweitingFile",default="")

doPtReWeighting=utl.getBoolFromConfigs(cfgTxt,"doPtReWeighting",default=False)
binnedReweitingFile=utl.getValueFromConfigs(cfgTxt,"pTReweitingFile",default="")

mlScoreTag=utl.getValueFromConfigs(cfgTxt,"mlScoreTag",default="")
methord=utl.getValueFromConfigs(cfgTxt,"methord",default="DHH")

doOverlapRemoval =utl.getValueFromConfigs(cfgTxt,"doOverlapRemoval",default="1") ; doOverlapRemoval = int(doOverlapRemoval) > 0.5
doBadBJetRemoval =utl.getValueFromConfigs(cfgTxt,"doBadBJetRemoval",default="1") ; doBadBJetRemoval = int(doBadBJetRemoval) > 0.5
isData           =utl.getValueFromConfigs(cfgTxt,"isData",default="1") ; isData = int(isData) > 0.5
etaMax           =float(utl.getValueFromConfigs(cfgTxt,"etaMax",default="2.5"))
pTMin            =float(utl.getValueFromConfigs(cfgTxt,"pTMin",default="25.0"))
overlapRemovalDRMax =float(utl.getValueFromConfigs(cfgTxt,"overlapRemovalDRMax",default="0.4"))
doBjetCounting=True

bkgSupressionBDTFile = utl.getValueFromConfigs(cfgTxt,"bkgSupressionBDTFile",default=None)
bkgSupressionTags = utl.getValueFromConfigs(cfgTxt,"bkgSupressionBDTTags",default=None)

bkgSupressionBDTFileOdds       = utl.getListOfStringsFromConfigs(cfgTxt,"#bkgSupressionBDTFileOdds_BEG" ,"#bkgSupressionBDTFileOdds_END" )
bkgSupressionBDTFileEvens      = utl.getListOfStringsFromConfigs(cfgTxt,"#bkgSupressionBDTFileEvens_BEG","#bkgSupressionBDTFileEvens_END")
bkgSupressionBDTTagsExclusives = utl.getListOfStringsFromConfigs(cfgTxt,"#bkgSupressionBDTTagsExclusives_BEG","#bkgSupressionBDTTagsExclusives_END")

exportedModelsFileNames       = utl.getListOfStringsFromConfigs(cfgTxt,"#ExportedModels_BEG" ,"#ExportedModels_END" )
exportedModelTags             = utl.getListOfStringsFromConfigs(cfgTxt,"#ExportedModelTags_BEG" ,"#ExportedModelTags_END" )
exportedModelClassCounts      = utl.getListOfStringsFromConfigs(cfgTxt,"#ExportedModelClassCount_BEG" ,"#ExportedModelClassCount_END" ) 

doBTagScales =utl.getValueFromConfigs(cfgTxt,"doBTagScales",default="0") ; doBTagScales = int(doBTagScales) > 0.5

def getScore(model,features,allVarDict,classIdx=[]):
    x=np.array([allVarDict[var] for var in features ])
    scores = model.predict_proba([x])[0]
    
    if classIdx==[]:
        return scores
    else:
        return [ scores[i] for i in classIdx]


print("allFnames   :  ",              allFnames)
print("foutName   :  ",               foutName)
print("processID   :  ",              processID)
print("treeName   :  ",               treeName)
print("weightScale   :  ",            weightScale)
print("doBDTReWeighting   :  ",        doPtReWeighting)
print("bdtReweitingFile   :  ",        bdtReweitingFile)
print("doPtReWeighting   :  ",        doPtReWeighting)
print("pTReweitingFile   :  ",        binnedReweitingFile)
print('etaMax   :       ', etaMax)
print('pTMin    :       ', pTMin)
print('drMax    :       ', drMax)
print('methord    :       ', methord)
print('mlScoreTag    :       ', mlScoreTag)
print('unblind    :       ', doSR)
print('bkgSupressionBDTFile : ',bkgSupressionBDTFile)
print('bkgSupressionBDTTagsExclusives : ',bkgSupressionBDTTagsExclusives)
print('bkgSupressionBDTFileOdds  : ',bkgSupressionBDTFileOdds)
print('bkgSupressionBDTFileEvens : ',bkgSupressionBDTFileEvens)

scaler=mcScaler()
print(scaleFactorUtil.__file__)
bdt_scaler= scaleFactorUtil.bdtScaler()

unblind = doSR

fname=allFnames[0]
dataTag,dataId   = hhhUtil.getDataTag(fname)
dataTag,sampleID = hhhUtil.getDataTag(fname)
#if doPtReWeighting:
binnedReweitingF=None
binnedScaleFactorHist=None
if doPtReWeighting and hhhSelector.hasToRePtWeight( dataTag ):
    print("Loading the binned Reweighter model  from : ",binnedReweitingFile)
    binnedReweitingF=ROOT.TFile(binnedReweitingFile)
    binnedScaleFactorHist=binnedReweitingF.Get('histScaleFactor')
    binnedScaleFactorHist.Print()
    #f.Close()
print("")


#if doBDTReWeighting:
if doBDTReWeighting and hhhSelector.hasToBDTReWeight( dataTag ):
    result={}
    with open(bdtReweitingFile,'rb') as f:
        print("Loading the BDT Reweighter model  from : ",bdtReweitingFile)
        results = pickle.load(f)
        bdt_scaler_=results['model']
        bdt_scaler=bdt_scaler_
        bdt_scaler.lumi=float( bdt_scaler.lumi)  ## patch
        bdt_scaler.evalTestEvent()
        #bdt_scaler.loadFromScaler(bdt_scaler_)
        #bdt_scaler.var_list = [
        #        'leadingPhoton_pt','leadingPhoton_eta',
        #        'subleadingPhoton_pt','subleadingPhoton_eta',
        #        'pT_h1leadJ' ,'eta_h1leadJ',
        #        'pT_h1subleadJ' ,'eta_h1subleadJ',
        #        'pT_h2leadJ' ,'eta_h2leadJ',
        #        'pT_h2subleadJ' ,'eta_h2subleadJ'
        #    ] 
modelDict={}
supressionModel=None
if bkgSupressionBDTFile:
    supressionModel={}
    fileNmaes=bkgSupressionBDTFile.split(',')
    tagNmaes=bkgSupressionTags.split(',')
    if len(tagNmaes)!=len(fileNmaes):
        print("BDT TAGS and File Names are not one-to-one")
        exit(1)
    for fName,modelTag in zip(fileNmaes,tagNmaes):  
        print(f"Loading BDT {modelTag} from {fName}")
        with open(fName, 'rb') as f:
            modelDict = pickle.load(f)
            supressionModel[modelTag]=clsUtl.xgb_bdt()
            supressionModel[modelTag].setModel(modelDict['model'],modelDict['variables'])

modelDictExclusive={}
supressionModelExclusives=None

print()
print()
print(bkgSupressionBDTTagsExclusives)
print()
if bkgSupressionBDTTagsExclusives:
    supressionModelExclusives={}
    fileNmaesOdd=bkgSupressionBDTFileOdds
    fileNmaesEvens=bkgSupressionBDTFileEvens
    tagNmaes=bkgSupressionBDTTagsExclusives
    
    if len(tagNmaes)!=len(fileNmaesOdd):
        print("BDT TAGS and File Names are not one-to-one")
        exit(1)
    if len(fileNmaesOdd)!=len(fileNmaesEvens):
        print("BDT TAGS and File Names are not one-to-one")
        exit(1)
    for fNameO,fNameE,modelTag in zip(fileNmaesOdd,fileNmaesEvens,tagNmaes):  
        print(f"Loading BDT {modelTag} from {fNameO} , {fNameE}")
        supressionModelExclusives[modelTag]={}
        with open(fNameO, 'rb') as f:
            modelDict = pickle.load(f)
            supressionModelExclusives[modelTag]['odd']=clsUtl.xgb_bdt()
            supressionModelExclusives[modelTag]['odd'].setModel(modelDict['model'],modelDict['variables'])
        with open(fNameE, 'rb') as f:
            modelDict = pickle.load(f)
            supressionModelExclusives[modelTag]['even']=clsUtl.xgb_bdt()
            supressionModelExclusives[modelTag]['even'].setModel(modelDict['model'],modelDict['variables'])



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
     'run','event',
     'CMS_hgg_mass','hh_mass','hh_pt','hh_eta','hh_phi', 'D_HH', 'H1H2JetAbsCosThetaMax', 'H1H2JetAbsCosThetaMin','nJet',
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
     'pT_h1leadJ','pT_h1subleadJ', 'pT_h2leadJ', 'pT_h2subleadJ',
     'eta_h1leadJ','eta_h1subleadJ', 'eta_h2leadJ', 'eta_h2subleadJ',
     'phi_h1leadJ','phi_h1subleadJ', 'phi_h2leadJ', 'phi_h2subleadJ',
     'pTsubleadG_overMgg', 'quadjet_0_deepJetScore', 'quadjet_1_deepJetScore',
     'quadjet_2_deepJetScore', 'quadjet_3_deepJetScore', 'r_HH', 'scalarPtSum4b',
     'scalarPtSum4b2g', 'scalarPtSum4b2g','HggTo4bAbsCosTheta', 'scalarPtSumHHH',
     'h1_dijetSigmaMOverM', 'h2_dijetSigmaMOverM', 'sigmaMOverM',
     'subleadingPhotonSigOverE', 'sumScore_3j', 'sumScore_4j', 'sumScore_4j_Pnet','sumScore_3j_Pnet', 'trihiggs_mass',
     'ttH_MET', 'weight','rho','trihiggs_pt',
     "ttH_sumET", "ttH_MET", "ttH_phiMET", "ttH_dPhi1", "ttH_dPhi2", "ttH_PhoJetMinDr", 
     "ttH_njets", "ttH_Xtt0", "ttH_Xtt1", "ttH_pte1", "ttH_pte2", "ttH_ptmu1", "ttH_ptmu2", 
     "ttH_ptdipho", "ttH_etae1", "ttH_etae2", "ttH_etamu1", "ttH_etamu2", "ttH_etadipho", "ttH_phie1", "ttH_phie2", "ttH_phimu1", 
     "ttH_phimu2", "ttH_phidipho", "ttH_fabs_CosThetaStar_CS", "ttH_fabs_CosTheta_bb", "ttH_ptjet1", "ttH_ptjet2", "ttH_etajet1", 
     "ttH_etajet2", "ttH_phijet1", "ttH_phijet2", "ttHScore",
     "CosThetaH1_hhhF",  "HH4bCosTheta_hhhF",  "HggCosTheta_hhhF",  "HH4bCosThetaLeadJet_hhhF",  "absCosThetaH4bHgg_hhhF",  
     'H1bbCosTheta_hhhF','H2bbCosTheta_hhhF','dZ',
     'leadingPhoton_pt','leadingPhoton_eta','leadingPhoton_phi',
     'subleadingPhoton_pt','subleadingPhoton_eta','subleadingPhoton_phi',
     'diphoton_pt','diphoton_eta','diphoton_phi','weight_v0','weight_bdt','weight_binned',
     'weight_v0', 'weight','scale_factor','lumi','lumiID','year','yearID','dataTag','sample','scale_factor_bdt','scale_factor_binned',
     'corrMET','corrMET','corrMETPhi','diPhoMVAValueHHH'
] 
for i in range(4):
    branches.append('quadjet_'+str(i)+'_puJetIdMVA')
    branches.append('quadjet_'+str(i)+'_mass')
    branches.append('quadjet_'+str(i)+'_flavour')
    branches.append('quadjet_'+str(i)+'_bJetRegCorr')
    branches.append('quadjet_'+str(i)+'_bJetRegRes')
    branches.append('quadjet_'+str(i)+'_flavour')
    branches.append('quadjet_'+str(i)+'_isLoose')
    branches.append('quadjet_'+str(i)+'_isTight' )
    branches.append('quadjet_'+str(i)+'_isTight2017')
    branches.append('quadjet_'+str(i)+'_isTight2018' )
    branches.append('quadjet_'+str(i)+'_particleNetAK4_B')
    branches.append('quadjet_'+str(i)+'_isPUTight' )
    branches.append('quadjet_'+str(i)+'_isPUMedium' )
    branches.append('quadjet_'+str(i)+'_isPULoose' )
    branches.append('quadjet_'+str(i)+'_isBtagTight' )
    branches.append('quadjet_'+str(i)+'_isBtagMedium' )
    branches.append('quadjet_'+str(i)+'_isBtagLoose' )
branches+=['sumPNetScore_3j' ,'sumPNetScore_4j']



branches+=['quadjet_0_mlScore', 'quadjet_0_mlScoreY1s', 'quadjet_0_mlScoreY2s' ,
    'quadjet_1_mlScore', 'quadjet_1_mlScoreY1s', 'quadjet_1_mlScoreY2s' ,
    'quadjet_2_mlScore', 'quadjet_2_mlScoreY1s', 'quadjet_2_mlScoreY2s' ,
    'quadjet_3_mlScore', 'quadjet_3_mlScoreY1s', 'quadjet_3_mlScoreY2s' ]
branches+=['genRecoCategory']    
if supressionModel:
    branches+=list(supressionModel.keys())
if supressionModelExclusives:
    for ky in supressionModelExclusives:
        branches+= [ ky +'_even' ]
        branches+= [ ky +'_odd' ]
        branches+= [ ky  ]

pickedModels={}
for mName,fname,nClass in zip(exportedModelTags,exportedModelsFileNames,exportedModelClassCounts):
    with open(fname,'rb') as f:
        print( "Adding model ",mName," from ",fname , " with ",nClass)
        modelDict=pickle.load(f)
        pickedModels[mName] ={
                'model'    : modelDict['classifier'],
                'features' : modelDict['features']
        }
        for i in range(int(nClass)):
            branches.append(f"{mName}_{i}")

btvFiles={}
doBTagScales=True
if doBTagScales:
    keyMap={
              '2016'  :  '/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/data/jsonpog-integration/POG/BTV/2016preVFP_UL/btagging.json',
              '2016PreVFP'  :  '/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/data/jsonpog-integration/POG/BTV/2016preVFP_UL/btagging.json',
              '2016PostVFP' :  '/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/data/jsonpog-integration/POG/BTV/2016postVFP_UL/btagging.json',
              '2017'        :  '/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/data/jsonpog-integration/POG/BTV/2017_UL/btagging.json',
              '2018'        :  '/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/data/jsonpog-integration/POG/BTV/2018_UL/btagging.json'
          }
    for ky in keyMap:
        sfName=keyMap[ky]
        print("Opening BTV Shape scale factor file : ",sfName)
        btvFiles[ky] = _core.CorrectionSet.from_file(sfName)
        print(f" > Added the {sfName} for {ky}")
    #for v in ["up","down"]:
    #    for syst in ["lf", "hf", "hfstats1", "hfstats2", "lfstats1","lfstats2"]:
    #        syskey=f"{v}_{syst}"
    #        ucert=btvjson["deepJet_shape"].evaluate(syskey, flav,eta,pt,score)
    for fname in allFnames:
        print("Opening file : ",fname)
        if fname[-3:]=='roo':
            fname+='t'
            print("Ading bugfix")
        year,yearID     =hhhUtil.getYear(fname)
    
    JEC_SPLIT=[]
    if '18' in year:
        JEC_SPLIT=['jesAbsolute_2018','jesBBEC1_2018','jesEC2_2018','jesHF_2018','jesRelativeSample_2018']
    if '17' in year:
        JEC_SPLIT=['jesAbsolute_2017','jesBBEC1_2017','jesEC2_2017','jesHF_2017','jesRelativeSample_2017']
    if '16' in year:
        JEC_SPLIT=['jesAbsolute_2016','jesBBEC1_2016','jesEC2_2016','jesHF_2016','jesRelativeSample_2016']

    branches+=['btag_scale']
    for v in ["up","down"]:
        for syst in ["jes","lf", "hf", "hfstats1", "hfstats2", "lfstats1","lfstats2","cferr1","cferr2"]+JEC_SPLIT:
            branches+=[f'btag_scale_{v}_{syst}']

btagGlobalScalesFile="/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/results/analysisSystNtuples/v32p1_07AprSunday05PM04m37s/v32p1/btag_globalScal.json"
print("Loading the btagGlobalScalesFile ",btagGlobalScalesFile)
with open(btagGlobalScalesFile)  as f:
    btagNormalizationMap=json.load(f)

ntuple={}
branches=np.unique(branches)
#print(branches)

ntuple['eventTree']  = ROOT.TNtuple(outTreeName, outTreeName, ':'.join(branches))

tofill = OrderedDict(zip(branches, [np.nan]*len(branches)))
outputDataDict=OrderedDict(zip(branches, [np.nan]*len(branches))) 

kyList = [ ky for ky in tofill ]

print("len(branches) : " , len(branches))
print("len(tofill) : " , len(tofill) )

sumEntries=ROOT.TH1F("sumEvts","sumEvts",1,0.0,1.0)
sumEntries.SetCanExtend(ROOT.TH1.kAllAxes)
sumWeights=ROOT.TH1F("sumWeighs","sumWeighs",1,0.0,1.0)
sumWeights.SetCanExtend(ROOT.TH1.kAllAxes)
isDataHist=ROOT.TH1F("isData","isData",1,0.0,1.0)
isDataHist.SetCanExtend(ROOT.TH1.kAllAxes)

th1Store={}


beg=datetime.datetime.now()
has_printed=True
dataTag=""
countTotal=0
for fname in allFnames:
    print("Opening file : ",fname)
    if fname[-3:]=='roo':
        fname+='t'
        print("Ading bugfix")
    simFile = ROOT.TFile(fname,'READ')
    dataTag,dataId=hhhUtil.getDataTag(fname)
    dataTag,sampleID=hhhUtil.getDataTag(fname)
    year,yearID     =hhhUtil.getYear(fname)
    proc=hhhUtil.getProc(fname)
    if year:
        print(year)
        print(utl.lumiMap)
        print(utl.lumiMap[year])
        lumi=float(utl.lumiMap[year])
    if isData:
        lumi=1

    print(" data tag : "  , dataTag )
    print(" proc     : "  , proc    )
    print(" year     : "  , year    )
    print(" lumi     : "  , lumi    )
    btvjson=btvFiles[year]
    global_btagNormalization=1.0
    if proc!='data':
        global_btagNormalization=btagNormalizationMap[proc][year]["globalScale"]
    print("Gloabl btag Normalization loaded as  : ",global_btagNormalization)
    eTree=simFile.Get(treeName)
    if not eTree:
        treeName='tagsDumper/bkg_13TeV_TrippleHTag_0'
        print("Trying tree name : ",treeName)
        eTree=simFile.Get(treeName)
        if eTree:
            print("    got it ~")
    if not eTree:
        treeName='tagsDumper/trees/bkg_13TeV_TrippleHTag_0'
        print("Trying tree name : ",treeName)
        eTree=simFile.Get(treeName)
    if not eTree:
        treeName='tagsDumper/trees/ggHHH_125_13TeV'
        print("Trying tree name : ",treeName)
        eTree=simFile.Get(treeName)
    if not eTree:
        treeName='tagsDumper/trees/Data_13TeV_TrippleHTag_0'
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
    print(f'{len(tofill)=}')
    for i in range(maxEvents_):
        eTree.GetEntry(i)
        countTotal+=1
        if(i%250==0):
            now=datetime.datetime.now()
            timeSpendSec=np.round((now-beg).total_seconds() , 2)
            timeLeftSec =np.round(1.0*(maxEvents_-i)*timeSpendSec/( i +1e-3),2)
            print("   Doing i = ",i," / ",maxEvents_," in file  [  ", countTotal ," total ] ")
            print("      time left : ", str(datetime.timedelta(seconds= timeLeftSec)),
                    " [ time elapsed : ",datetime.timedelta(seconds= timeSpendSec), " s ]")
            print(" gg mass   : ",eTree.CMS_hgg_mass)        
        
        wei_base=eTree.weight
        wei=eTree.weight*lumi
        wei_binned = wei 
        wei_bdt    = wei
        if doWeight:
            if resetWeight > -1e4:
                wei=resetWeight
            if weightScale > -1e4:
                wei*=weightScale

        scaleFactor_binned = 1.0
        if doPtReWeighting and hhhSelector.hasToRePtWeight( dataTag ):
            pT=eTree.diphoton_pt
            scaleFactor_binned=scaleFactorUtil.getBinnedWeight( binnedScaleFactorHist,[pT])[0]
            wei_binned*=scaleFactor_binned
        #    scaleFactor_binned=scaler.getSFForX(pT,'inclusive')
        #    print("scaleFactor_binned : ",scaleFactor_binned,"[  ",pT," ]")
        sumEntries.Fill('total', 1)
        sumWeights.Fill('total', wei)

        isMasked = isData_ and  (eTree.CMS_hgg_mass > 115.0) and (eTree.CMS_hgg_mass < 135.0)
        if isMasked and not unblind:
            continue
        jetMask=hhhSelector.getCorrectJetCollection(eTree,jetMask=[]) 
        ##   JET PRE-SELECTION
        jetMask=hhhSelector.getSelectedJetCollectionMaskEta(eTree,jetMask=jetMask,etaMax=etaMax)
        #for ii in range(8):
        #    print(f"{ii}. pt {getattr(eTree,'jet_'+str( ii )+'_pt'):.3f} , eta  {getattr(eTree,'jet_'+str( ii )+'_eta'):.3f}")
        #print(jetMask)
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
        if doBadBJetRemoval:
            jetMask=hhhSelector.getSelectedJetCollectionMaskBadBJet(eTree,jetMask=jetMask)
            if sum(jetMask) < 4 :
                sumEntries.Fill('nJetPreselectionBadBJet',1)
                sumWeights.Fill('nJetPreselectionBadBJet',wei)
                continue
        
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
        tofill['lumiID']=tofill['lumi']
        verbose=False
        sumScore_3j=0.0
        sumScore_4j=0.0

        if methord=='DHH':
            allQuads=hhhSelector.getBJetParisFGG(eTree,year=year,mask=jetMask)
            if allQuads['isValid']:
                quad=allQuads['bJetQuad']
                i=0
                for idx in quad['fgg_idxs']:
                    val=max(getattr(eTree,'jet_'+str( idx )+'_deepJetScore'),0.0)
                    tofill['quadjet_'+str(i)+'_deepJetScore']   = val
                    if i< 3:
                        sumScore_3j+=val
                    sumScore_4j+=val
                    i+=1
            else:
                continue
        elif methord=='mha':
            allQuads=hhhSelector.getBJetParis_wrapper(eTree, mask= jetMask , methord='mha',mlScoreTag=mlScoreTag,threshold=-1e3,doMisclassificationCorrection = False)
            if not allQuads['isValid']:
                print("anity check failed after the quad finding")
                exit(1)
                continue
            if allQuads['jetsValid']:
                i=0
                for idx in allQuads['allJetsSelected']:
                    tofill['quadjet_'+str(i)+'_deepJetScore']   =getattr(eTree,'jet_'+str( idx )+'_deepJetScore') 
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
        i=0
        sumScore_3j_Pnet=0.0
        sumScore_4j_Pnet=0.0
        for idx in quad['fgg_idxs']:
           tofill['quadjet_'+str(i)+'_puJetIdMVA']  = getattr(eTree,'jet_'+str( idx )+'_puJetIdMVA')
           tofill['quadjet_'+str(i)+'_mass']        = getattr(eTree,'jet_'+str( idx )+'_mass')
           tofill['quadjet_'+str(i)+'_pt']     = getattr(eTree,'jet_'+str( idx )+'_pt')
           tofill['quadjet_'+str(i)+'_eta']     = getattr(eTree,'jet_'+str( idx )+'_eta')
           tofill['quadjet_'+str(i)+'_flavour']     = getattr(eTree,'jet_'+str( idx )+'_flavour')
           tofill['quadjet_'+str(i)+'_bJetRegCorr'] = getattr(eTree,'jet_'+str( idx )+'_bJetRegCorr')
           tofill['quadjet_'+str(i)+'_bJetRegRes']  = getattr(eTree,'jet_'+str( idx )+'_bJetRegRes')
           tofill['quadjet_'+str(i)+'_flavour']     = getattr(eTree,'jet_'+str( idx )+'_flavour')
           tofill['quadjet_'+str(i)+'_isLoose']     = getattr(eTree,'jet_'+str( idx )+'_isLoose')
           tofill['quadjet_'+str(i)+'_isTight']     = getattr(eTree,'jet_'+str( idx )+'_isTight')
           tofill['quadjet_'+str(i)+'_isTight2017'] = getattr(eTree,'jet_'+str( idx )+'_isTight2017')
           tofill['quadjet_'+str(i)+'_isTight2018'] = getattr(eTree,'jet_'+str( idx )+'_isTight2018')
           tofill['quadjet_'+str(i)+'_particleNetAK4_B'] = getattr(eTree,'jet_'+str( idx )+'_particleNetAK4_B')
           
           tofill['quadjet_'+str(i)+'_isPULoose' ]     = hhhUtil.puJetID(getattr(eTree,'jet_'+str( idx )+'_pt'),getattr(eTree,'jet_'+str( idx )+'_puJetIdMVA'),year,'loose' )
           tofill['quadjet_'+str(i)+'_isPUMedium']     = hhhUtil.puJetID(getattr(eTree,'jet_'+str( idx )+'_pt'),getattr(eTree,'jet_'+str( idx )+'_puJetIdMVA'),year,'medium')
           tofill['quadjet_'+str(i)+'_isPUTight' ]     = hhhUtil.puJetID(getattr(eTree,'jet_'+str( idx )+'_pt'),getattr(eTree,'jet_'+str( idx )+'_puJetIdMVA'),year,'tight' )

           tofill['quadjet_'+str(i)+'_isBtagLoose' ]     = hhhUtil.btagID(getattr(eTree,'jet_'+str( idx )+'_deepJetScore'),'loose' )
           tofill['quadjet_'+str(i)+'_isBtagMedium']     = hhhUtil.btagID(getattr(eTree,'jet_'+str( idx )+'_deepJetScore'),'medium')
           tofill['quadjet_'+str(i)+'_isBtagTight' ]     = hhhUtil.btagID(getattr(eTree,'jet_'+str( idx )+'_deepJetScore'),'tight' )
           
           val=max(getattr(eTree,'jet_'+str( idx )+'_particleNetAK4_B'),0.0) ;  
           
           if i< 3:
               sumScore_3j_Pnet+=val
           sumScore_4j_Pnet+=val
           i+=1



        genCat=0
        if (processID=='sig') and ( hasattr(eTree,'gen_H1_dau1_pdgId') ) and False :
            cat=1
            rslt=hhhSelector.getBestGetMatchesGlobalFGG(eTree,jetMask=jetMask)
            idxs=rslt['fgg_idxs']
            isRecoed={}
            for i in range(4):
                isRecoed['isRecoed_'+str(i)]=True
            for i in range(4):
                if idxs[i] < 0:
                    isAllRecoed=False
                    isRecoed['isRecoed_'+str(i)]=False
            
            h1Out=  isRecoed['isRecoed_0'] and isRecoed['isRecoed_1']
            h2Out=  isRecoed['isRecoed_2'] and isRecoed['isRecoed_3']
            
            nout=0
            if not isRecoed['isRecoed_0']:            nout+=1;    
            if not isRecoed['isRecoed_1']:            nout+=1;    
            if not isRecoed['isRecoed_2']:            nout+=1;    
            if not isRecoed['isRecoed_3']:            nout+=1;

            x=0
            for i in range(8):
                if i in idxs:
                    x+=1
            cat=1
            if nout==0 and x==4 :
                cat=2 #'allRecoed'
            elif nout==1 and x==3 :
                cat=3 #' SingleJetMiss |'
            if nout==2 and x==2  :
                cat=4 #' 2 JetMiss |'
            if nout==3 and x==1  :
                cat=5 #' 3 JetMiss |'
            genCat=cat
        
        tofill['genRecoCategory']=genCat*1.0
        tofill["sumScore_4j"]=sumScore_4j
        tofill["sumScore_3j"]=sumScore_3j
        tofill["sumScore_4j_Pnet"]=sumScore_4j_Pnet
        tofill["sumScore_3j_Pnet"]=sumScore_3j_Pnet
        tofill['minPhoIDMVA'] =  min(eTree.customSubLeadingPhotonIDMVA,eTree.customLeadingPhotonIDMVA)
        tofill['maxPhoIDMVA'] =  max(eTree.customSubLeadingPhotonIDMVA,eTree.customLeadingPhotonIDMVA)
        tofill['scalarPtSum4b2gAvg'] =  eTree.leadingPhoton_pt + eTree.subleadingPhoton_pt
        for idx in quad['fgg_idxs']:
            tofill['scalarPtSum4b2gAvg'] += getattr(eTree,'jet_'+str( idx )+'_pt')
        tofill['scalarPtSum4b2gAvg'] /= 6.0

        quad=allQuads['bJetQuad']

        if doBjetCounting:
            nb=hhhUtil.getNBsFromQuad(eTree,quad)
            if hhhSelector.vetoOverCountings(nb,dataTag):
                sumEntries.Fill('nFailedFromOverCounting',1)
                sumWeights.Fill('nFailedFromOverCounting',wei)
                continue


        LVStore = hhhUtil.getLVStoreFromTreeAndQuad(eTree,quad)
        j1CosTheta,k1CosTheta,ggCostheta,drMin,drOther=hhhUtil.getCosthetaVars(eTree,LVStore)
        
        tofill['r_HH'] = quad['r_HH']
        tofill['D_HH'] = quad['D_HH']

        #hhhUtil.printEventInfo(eTree,jetMask=jetMask)
        #hhhUtil.printEventInfoFromLVStore(LVStore,eid=eTree.event)
        #print("\n")
        #print("fgg indexs selected : ",quad['fgg_idxs'])
        #for ky in tofill:
        #    if ky not in varDict:
        #       print("filling default barnch value  for branch",ky)
        varDict=hhhUtil.getOtherDerivedVariables(eTree,LVStore,quad)
        for ky in varDict:
            tofill[ky]=varDict[ky]
        tofill['ttH_MET'] = eTree.ttH_MET


        for ky in outputDataDict:
            if ky not in tofill:
                print("not in tofill ",ky)
        
        scaleFactor_bdt=1.0
        if doBDTReWeighting and hhhSelector.hasToBDTReWeight( dataTag ):
            scaleFactor_bdt=bdt_scaler.getSFForX(tofill)[0]
            wei_bdt=scaleFactor_bdt #*lumi/bdt_scaler.lumi
        
        if supressionModel:
            for modelTag in supressionModel: 
                tofill[modelTag]=supressionModel[modelTag].predict_probaForX(tofill)[1]
        
        if supressionModelExclusives:
            evtIdNpy = np.array([ eTree.event  ])[0]
            hashId = int(hashlib.sha256( evtIdNpy ).hexdigest()[:-1],16)%2
            for modelTag in supressionModelExclusives: 
                tofill[modelTag+'_odd' ]=supressionModelExclusives[modelTag]['odd'].predict_probaForX(tofill)[1]
                tofill[modelTag+'_even']=supressionModelExclusives[modelTag]['even'].predict_probaForX(tofill)[1]
                #print(" event : ",eTree.event," ", tofill[modelTag+'_odd']  , tofill[modelTag+'_even']  )
                if hashId==0:
                    #print("\t even : ",tofill[modelTag+'_even'])
                    tofill[modelTag] = tofill[modelTag+'_even'] 
                else:
                    #print("\t odd : ",tofill[modelTag+'_odd'] )
                    tofill[modelTag] = tofill[modelTag+'_odd'] 
        for mName in pickedModels:
            #tofillX={ky:1.0+ij  for ij,ky in enumerate(pickedModels[mName]['features']) }
            #print(mName,"  :  ",[ tofill[ky]  for ij,ky in enumerate(pickedModels[mName]['features']) ])
            scores=getScore( pickedModels[mName]['model'] , pickedModels[mName]['features'], tofill)
            #print(mName)
            #print(tofillX)
            #print(scores)
            #print(mName , scores)
            for ii in range(len(scores)):
                name=f'{mName}_{ii}'
                tofill[name]=scores[ii]
        
        jet_sf=1.0*global_btagNormalization
        tofill['btag_scale']=global_btagNormalization
        #print("Global normalization set to : ",global_btagNormalization)
        for v in ["up","down"]:
            for syst in ["jes","lf", "hf", "hfstats1", "hfstats2", "lfstats1","lfstats2","cferr1","cferr2"]:
                tofill[f'btag_scale_{v}_{syst}']=global_btagNormalization
        
        if proc!='data':
            JEC_SPLIT=[]
            if '18' in year:
                JEC_SPLIT=['jesAbsolute_2018','jesBBEC1_2018','jesEC2_2018','jesHF_2018','jesRelativeSample_2018']
            if '17' in year:
                JEC_SPLIT=['jesAbsolute_2017','jesBBEC1_2017','jesEC2_2017','jesHF_2017','jesRelativeSample_2017']
            if '16' in year:
                JEC_SPLIT=['jesAbsolute_2016','jesBBEC1_2016','jesEC2_2016','jesHF_2016','jesRelativeSample_2016']
            for v in ["up","down"]:
                for syst in JEC_SPLIT:
                    tofill[f'btag_scale_{v}_{syst}']=1.0
            for i in range(4):
                flav,eta,pt,score=int(tofill['quadjet_'+str(i)+'_flavour']), abs(tofill['quadjet_'+str(i)+'_eta']), tofill['quadjet_'+str(i)+'_pt'], tofill['quadjet_'+str(i)+'_deepJetScore']
                if eta> (2.5 -1e-7):
                    eta=2.5-1e-8
                #if flav < 5:
                #    flav=0

                jet_sf_ = btvjson["deepJet_shape"].evaluate("central", flav,eta,pt,score)    
                jet_sf*=jet_sf_
                tofill['btag_scale']*=jet_sf_
                #print("Jet SF : ",jet_sf," ( ",jet_sf_," ) ")
                #print("   > ucerts : ",)
                for v in ["up","down"]:
                    for syst in [ "lf", "jes", "hf", "hfstats1", "hfstats2", "lfstats1", "lfstats2", "cferr1", "cferr2"]+JEC_SPLIT:
                        syskey=f"{v}_{syst}"
                        if ( abs(flav)==4) and ( syst not in ["cferr1", "cferr2"]  ) :
                            #print( syskey,flav  )
                            continue
                        if ( abs(flav)!=4) and ( syst in ["cferr1", "cferr2"]  )  :
                            continue
                        ucert=btvjson["deepJet_shape"].evaluate(syskey, flav,eta,pt,score)
                        tofill[f'btag_scale_{v}_{syst}']*=ucert
                tofill['quadjet_'+str(i)+'_flavour']
                        #print("     > ",syskey," -> ",ucert)
        
        #print("\njet SF for shape correction:")
        #print("SF: {jet_sf}".format(jet_sf=jet_sf))
        #for ky in tofill:
            #if 'quadjet_1' in ky:
        #    if True:
        #         print(ky)
        #exit(1)

        wei=wei_bdt*jet_sf
        scaleFactor=scaleFactor_bdt


        tofill['nJet'] =  sum(jetMask)
        tofill['weight'] =  wei
        tofill['weight_v0'] =  wei_base
        tofill['weight_binned'] =  wei_binned
        tofill['weight_bdt'] =  wei_bdt
        tofill['scale_factor'] =  scaleFactor
        tofill['scale_factor_bdt'] =  scaleFactor_bdt
        tofill['scale_factor_binned'] =  scaleFactor_binned
        tofill['lumi'] =  lumi
        tofill['yearID'] =  yearID
        tofill['sample'] =  sampleID
        #print(f"weight : {tofill['weight_v0'] = } {tofill['weight'] = } {tofill['weight_binned'] = } {tofill['weight_bdt'] =}")
        #print()
        #print()
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
