from __future__ import print_function
import ROOT ,copy
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools as itrTools
from Util import *
import warnings
import trippleHiggsSelector as hhhSelector
N_JET_MAX=8
def sigmoid(x):
    return 1 / (1 + np.exp(-x))

class mcScaler:
    def __init__(self):
        self.dataHist=None
        self.mcHist=None
        self.binEdges=[]
        self.scaleFactors=[]
        self.def_scaleFactor=0
        self.underflow_x={}
        self.overflow_x={}
        self.underflow_scaleFactor={}
        self.overflow_scaleFactor={}
        self.scaleFactorHist={}
        self.scaleFactorHistSet={}
        self.File=None
    
    def setSFHist(self,hist,cat='def'):
        self.scaleFactorHist[cat]=hist.Clone()
        self.underflow_x[cat]= self.scaleFactorHist[cat].GetBinLowEdge(1)
        self.overflow_x[cat]= self.scaleFactorHist[cat].GetBinLowEdge(self.scaleFactorHist[cat].GetNbinsX())
        self.underflow_scaleFactor[cat] = self.scaleFactorHist[cat].GetBinContent(1)
        self.overflow_scaleFactor[cat] = self.scaleFactorHist[cat].GetBinContent(self.scaleFactorHist[cat].GetNbinsX()-1)
        print("  Setting scale factor for category " , cat  )
        print("     Scale factors defined for x = [ ",self.underflow_x[cat],self.overflow_x[cat]," ]")
        print("     Under flow Scale factor     =  ",self.underflow_scaleFactor[cat])
        print("     Over flow Scale factor      = ",self.overflow_scaleFactor[cat])

        self.scaleFactorHistSet[cat]=True
    
    def setSFHistFromFile(self,fname,histNames={}):
        self.File=ROOT.TFile(fname,"READ")
        print(histNames)
        for cat in histNames:
            print("CAT : ",cat," || Setting SF hist ",histNames[cat]," from ",fname)
            hist=self.File.Get(histNames[cat])
            self.setSFHist(hist,cat)
    
    def getSFForX(self,x,cat):
        if cat not in self.scaleFactorHistSet:
            print("\tWARNING ! : scaleFactor Hist not set , returning -1.0  [ cat : ",cat," available cats : ",self.scaleFactorHist.keys(),"]")
            return -1.0
        if x < self.underflow_x[cat]:
            return self.underflow_scaleFactor[cat]
        if x > self.overflow_x[cat]:
            return self.overflow_x[cat]
        bid=self.scaleFactorHist[cat].FindBin(x)
        #print("bid = ",bid," for x ",x, " Bin edges --> ",
        #         self.scaleFactorHist.GetBinLowEdge(bid), 
        #         self.scaleFactorHist.GetBinLowEdge(bid)+self.scaleFactorHist.GetBinWidth(bid))
        return self.scaleFactorHist[cat].GetBinContent(bid)

def btagID( deepJetScore,TYPE='loose' ):

    if TYPE=='tight' and  deepJetScore  > 0.7346:
        return True;
    if TYPE=='medium' and  deepJetScore  > 0.3098:
        return True;
    if TYPE=='loose' and  deepJetScore  > 0.0594:
        return True;

    return False    


def puJetID(pT,mvaScore,TYPE='loose'):
    
    if TYPE=='loose':
        return puJetIDLoose(pT,mvaScore)

    if TYPE=='tight':
        return puJetIDTight(pT,mvaScore)

    if TYPE=='medium':
        return puJetIDMedium(pT,mvaScore)

def	puJetIDLoose(pT,mvaScore):
	if pT > 50.0:
		return True
	if pT < 10.0:
		return False
	elif pT < 20.0:
		return mvaScore > -0.95
	elif pT < 30.0:
		return mvaScore > -0.90
	elif pT < 40.0:
		return mvaScore > -0.71
	elif pT < 50.0:
		return mvaScore > -0.40
	return True

def	puJetIDMedium(pT,mvaScore):
	if pT > 50.0:
		return True
	if pT < 10.0:
		return False
	elif pT < 20.0:
		return mvaScore > 0.20
	elif pT < 30.0:
		return mvaScore > 0.62
	elif pT < 40.0:
		return mvaScore > 0.86
	elif pT < 50.0:
		return mvaScore > 0.93
	return True

def	puJetIDTight(pT,mvaScore):
	if pT > 50.0:
		return True
	if pT < 10.0:
		return False
	elif pT < 20.0:
		return mvaScore > 0.71
	elif pT < 30.0:
		return mvaScore > 0.87
	elif pT < 40.0:
		return mvaScore > 0.94
	elif pT < 50.0:
		return mvaScore > 0.97
	return True
          
def getPtDependentScaleFactor(name="scaleFactor",dataHist=None,mcHists=None,binEdges=None,def_scaleFactor=0.0):
    mcHistSum=mcHists[0].Clone()
    mcHistSum.Reset()
    scaleFactorHist=ROOT.TH1F("","",1,0,100.0)
    if binEdges!=None:
        nBinEs=len(binEdges)
        bEdges=np.asarray(binEdges)
        scaleFactorHist=ROOT.TH1F(name,"",nBinEs-1,bEdges)
    else:
        scaleFactorHist=mcHists[0].Clone()
        scaleFactorHist.Reset()
    scaleFactorHist.SetName(name+"_scaleFactor")
    dataHTmp=scaleFactorHist.Clone()
    mcHTmp=scaleFactorHist.Clone()
    for h in mcHists:
            mcHistSum.Add(h)
    if dataHist.GetNbinsX() != mcHistSum.GetNbinsX():
        print("bin counts dont match !!")
        return None
    print("Number of bins in Data : ",dataHist.GetNbinsX())
    print("Number of bins in MC : ",mcHistSum.GetNbinsX())
    print("Width of 1st bin in Data : ",dataHist.GetBinWidth(1))
    print("Width of 1st bin in MC : ",mcHistSum.GetBinWidth(1))
    print("Number MC histograms : ",len(mcHists))
    scaler=mcScaler()        
    for i in range(1,dataHist.GetNbinsX()):
        x=dataHist.GetBinCenter(i)
        dValue=dataHist.GetBinContent(i)
        mcValue=mcHistSum.GetBinContent(i)
        dataHTmp.Fill(x,dValue)
        mcHTmp.Fill(x,mcValue)
        
    for i in range(1,dataHTmp.GetNbinsX()):
        dValue=dataHTmp.GetBinContent(i)
        mcValue=mcHTmp.GetBinContent(i)
        if mcValue==0:
            mcValue=1e-4
            print("\tWARNING : MC value found to be 0 for [",
                  dataHTmp.GetBinLowEdge(i),",",dataHTmp.GetBinLowEdge(i+1),
                  "]", "data value = ",dValue)
        scl=dValue/mcValue
        scaleFactorHist.SetBinContent(i,scl)
        scaler.binEdges.append(dataHTmp.GetBinLowEdge(i))  
        scaler.scaleFactors.append(scl)
    scaler.binEdges.append(dataHist.GetBinLowEdge(dataHist.GetNbinsX()+1))
    scaler.mcHist   = mcHTmp
    scaler.mcHist.SetName(name+"_mcHist")
    scaler.dataHist = dataHTmp
    scaler.dataHist.SetName(name+"_dataHist")
    scaler.setSFHist(scaleFactorHist)
    scaler.def_scaleFactor = 0.0
    return scaler

def getHiggsDauP4s(eTree,pdgId):
    return hhhSelector.getHiggsDauP4s(eTree,pdgId)
def addFlasggVars(histStore): 
    
    histStore["flashggVars"]={}

    histStore["flashggVars"]['HHHCosThetaH1']= ROOT.TH1F("HHHCosThetaH1","",60,-0.1,1.1)
    histStore["flashggVars"]['HHHCosThetaHgg']= ROOT.TH1F("HHHCosThetaHgg","",60,-0.1,1.1)
    histStore["flashggVars"]['HggCosThetaLeadGamma']= ROOT.TH1F("HggCosThetaLeadGamma","",60,-0.1,1.1)
    histStore["flashggVars"]['pTgg_overMgg'] = ROOT.TH1F("pTgg_overMgg"  ,"",100,0.0,1.0)
    histStore["flashggVars"]['pTleadG_overMgg'] = ROOT.TH1F("pTleadG_overMgg"  ,"",100,0.0,10.0)
    histStore["flashggVars"]['pTsubleadG_overMgg'] = ROOT.TH1F("pTsubleadG_overMgg"  ,"",100,0.0,10.0)
    histStore["flashggVars"]['rho'] = ROOT.TH1F("rho","",50,0.0,50.0)
    histStore["flashggVars"]['phoJetMinDr'] = ROOT.TH1F("phoJetMinDr"  ,"",60,0.0,6.0)
    histStore["flashggVars"]['otherphoJetMinDr'] = ROOT.TH1F("otherphoJetMinDr"  ,"",60,0.0,6.0)
    histStore["flashggVars"]['h1leadG_SigOverE'] = ROOT.TH1F("h1leadG_SigOverE"  ,"",50,0.0,0.1)
    histStore["flashggVars"]['h1subleadG_SigOverE'] = ROOT.TH1F("h1subleadG_SigOverE"  ,"",50,0.0,0.1)
    histStore["flashggVars"]['hgg_SigMOverM'] = ROOT.TH1F("hgg_SigMOverM"  ,"",20,0.0,0.1)
    histStore["flashggVars"]['hgg_SigMOverMDecorrelated'] = ROOT.TH1F("hgg_SigMOverMDecorrelated"  ,"",20,0.0,0.1)
    histStore["flashggVars"]['customLeadingPhotonIDMVA'] = ROOT.TH1F("customLeadingPhotonIDMVA"  ,"",48,-0.1,1.1)
    histStore["flashggVars"]['customSubLeadingPhotonIDMVA'] = ROOT.TH1F("customSubLeadingPhotonIDMVA"  ,"",48,-0.1,1.1)
    
    """
        d    PhoJetMinDr                         
        d    PhoJetOtherDr                       
        d    absCosTheta_gg                      
        d    customLeadingPhotonIDMVA            
        d    customSubLeadingPhotonIDMVA         
        d    h1LeadingJet_DeepFlavour            
        d    h1SubleadingJet_DeepFlavour         
        d    h1bbCosThetaLeadJet                 
        d    h2LeadingJet_DeepFlavour            
        d    h2SubleadingJet_DeepFlavour         
        d    h2bbCosThetaLeadJet                 
        d    hhhCosThetaH1                       
        d    leadingPhotonSigOverE               
        d    pTh1leadJ_overMh1                   
        d    pTh2leadJ_overMh2                   
        d    pTleadG_overMgg                     
        d    rho                                 
        d    sigmaM1Jets                         
        d    sigmaM2Jets                         
        d    sigmaMOverM                         
        d    subleadingPhotonSigOverE            
    """
def addOtherDerivedVariables(histStore):
    histStore["miscVars"]={}
    histStore["miscVars"]['HH4bCosThetaLeadJet']= ROOT.TH1F("HH4bCosThetaLeadJet","",60,-0.1,1.1)
    histStore["miscVars"]['H1bbCosThetaLeadJet'] = ROOT.TH1F("H1bbCosThetaLeadJet"  ,"",60,-0.1,1.1)
    histStore["miscVars"]['H2bbCosThetaLeadJet'] = ROOT.TH1F("H2bbCosThetaLeadJet"  ,"",60,-0.1,1.1)
    histStore["miscVars"]['pTh1jj_overMh1'] = ROOT.TH1F("pTh1jj_overMh1"  ,"",100,0.0,50.0)
    histStore["miscVars"]['pTh2jj_overMh2'] = ROOT.TH1F("pTh2jj_overMh2"  ,"",100,0.0,50.0)
    
    histStore["miscVars"]['pTh1leadJ_overMh1'] = ROOT.TH1F("pTh1leadJ_overMh1"  ,"",100,0.0,10.0)
    histStore["miscVars"]['pTh2leadJ_overMh2'] = ROOT.TH1F("pTh2leadJ_overMh2"  ,"",100,0.0,10.0)

    histStore["miscVars"]['pTh1subleadJ_overMh1'] = ROOT.TH1F("pTh1subleadJ_overMh1"  ,"",100,0.0,10.0)
    histStore["miscVars"]['pTh2subleadJ_overMh2'] = ROOT.TH1F("pTh2subleadJ_overMh2"  ,"",100,0.0,10.0)

    histStore["miscVars"]['h1leadJ_deepjetScore'] = ROOT.TH1F("h1leadJ_deepjetScore"  ,"",12,-0.1,1.1)
    histStore["miscVars"]['h1subleadJ_deepjetScore'] = ROOT.TH1F("h1subleadJ_deepjetScore"  ,"",12,-0.1,1.1)
    histStore["miscVars"]['h2leadJ_deepjetScore'] = ROOT.TH1F("h2leadJ_deepjetScore"  ,"",12,-0.1,1.1)
    histStore["miscVars"]['h2subleadJ_deepjetScore'] = ROOT.TH1F("h2subleadJ_deepjetScore"  ,"",12,-0.1,1.1)

    histStore["miscVars"]['h1_dijetSigmaMOverM'] = ROOT.TH1F("h1_dijetSigmaMOverM"  ,"",50,0.0,0.5)
    histStore["miscVars"]['h2_dijetSigmaMOverM'] = ROOT.TH1F("h2_dijetSigmaMOverM"  ,"",50,0.0,0.5)
    
    histStore["miscVars"]['triHiggs_mass']=ROOT.TH1F("triHiggs_mass","",1500,0.0,1500)
    histStore["miscVars"]['dijet1CandidatePtOverdiHiggsM'] = ROOT.TH1F("dijet1CandidatePtOverdiHiggsM"  ,"",100,0.0,50.0)
    histStore["miscVars"]['dijet2CandidatePtOverdiHiggsM'] = ROOT.TH1F("dijet2CandidatePtOverdiHiggsM"  ,"",100,0.0,50.0)
    histStore["miscVars"]["absCosThetaH4bHgg"    ]            = ROOT.TH1F("absCosThetaH4bHgg"  ,"",60,-0.1,1.1)
    histStore["miscVars"]["H1H2JetAbsCosThetaMax"]            = ROOT.TH1F("H1H2JetAbsCosThetaMax"  ,"",60,-0.1,1.1)
    histStore["miscVars"]["H1H2JetAbsCosThetaMin"]            = ROOT.TH1F("H1H2JetAbsCosThetaMin"  ,"",60,-0.1,1.1)
    histStore["miscVars"]["LeadJetAbsCosThetaMax"]            = ROOT.TH1F("LeadJetAbsCosThetaMax"  ,"",60,-0.1,1.1)
    histStore["miscVars"]["LeadJetAbsCosThetaMin"]            = ROOT.TH1F("LeadJetAbsCosThetaMin"  ,"",60,-0.1,1.1)
    histStore["miscVars"]["HggTo4bAbsCosTheta"]               = ROOT.TH1F("HggTo4bAbsCosTheta"  ,"",60,-0.1,1.1)
    histStore["miscVars"]["H1bbToH2bbAbsCosTheta"]            = ROOT.TH1F("H1bbToH2bbAbsCosTheta"  ,"",60,-0.1,1.1)        
    histStore["miscVars"]["LeadJetDrMaxWithOtherJets"]        = ROOT.TH1F("LeadJetDrMaxWithOtherJets","",50,0.0,5.0)
    histStore["miscVars"]["LeadJetDrMinWithOtherJets"]        = ROOT.TH1F("LeadJetDrMinWithOtherJets","",50,0.0,5.0)
    histStore["miscVars"]["H1H2JetDrMax"] = ROOT.TH1F("H1H2JetDrMax","",50,0.0,5.0) 
    histStore["miscVars"]["H1H2JetDrMin"] = ROOT.TH1F("H1H2JetDrMin","",50,0.0,5.0) 

    histStore["miscVars"]["pT_4b"]  =  ROOT.TH1F("pT_4b","",1500,0.0,1500) 
    histStore["miscVars"]["scalarPtSumHHH"]  = ROOT.TH1F("scalarPtSumHHH","",1500,0.0,1500)
    histStore["miscVars"]["scalarPtSum4b"]   = ROOT.TH1F("scalarPtSum4b","",1500,0.0,1500)
    histStore["miscVars"]["scalarPtSum4b2g"] = ROOT.TH1F("scalarPtSum4b2g","",1500,0.0,1500)
    histStore["miscVars"]["ttH_MET"] = ROOT.TH1F("ttH_MET","",1500,0.0,1500)

def getOtherDerivedVariables(eTree,LVStore,quad):
    varDict={}
    varDict['pTh1jj_overMh1']  = LVStore['H1bbLV'].Pt() / LVStore['H1bbLV'].M()
    varDict['pTh2jj_overMh2']  = LVStore['H2bbLV'].Pt() / LVStore['H2bbLV'].M()
    varDict['pThgg_overMgg']   = LVStore['HggLV'].Pt() / LVStore['HggLV'].M()

    varDict['pTleadG_overMgg']    = LVStore['g1LV'].Pt() / LVStore['HggLV'].M()
    varDict['pTsubleadG_overMgg'] = LVStore['g2LV'].Pt() / LVStore['HggLV'].M()
    varDict['pTh1leadJ_overMh1'] = LVStore['j1LV'].Pt() / LVStore['H1bbLV'].M()
    varDict['pTh1subleadJ_overMh1'] = LVStore['j2LV'].Pt()  / LVStore['H1bbLV'].M()
    varDict['pTh2leadJ_overMh2'] = LVStore['k1LV'].Pt() /LVStore['H2bbLV'].M()
    varDict['pTh2subleadJ_overMh2'] = LVStore['k2LV'].Pt()  / LVStore['H2bbLV'].M() 
    
    varDict["h1bbCosThetaLeadJet"]    =  abs(np.cos( LVStore['H1LV'].Angle(LVStore['j1LV'].Vect())))
    varDict['h1bb_pt'] = LVStore['H1bbLV'].Pt() 
    varDict['h2bb_pt'] = LVStore['H2bbLV'].Pt() 
    varDict['h1bb_eta'] = LVStore['H1bbLV'].Eta() 
    varDict['h2bb_eta'] = LVStore['H2bbLV'].Eta() 
    varDict['h1bb_phi'] = LVStore['H1bbLV'].Phi() 
    varDict['h2bb_phi'] = LVStore['H2bbLV'].Phi() 
    varDict['h1bb_mass'] = LVStore['H1bbLV'].M() 
    varDict['h2bb_mass'] = LVStore['H2bbLV'].M() 
    

    varDict["HH4bCosThetaLeadJet"]    =  abs(np.cos( LVStore['HH4bLV'].Angle(LVStore['j1LV'].Vect())))
    varDict["absCosThetaH4bHgg"]    =  abs(np.cos( (LVStore['H1bbLV'] + LVStore['H2bbLV']).Angle(LVStore['HggLV'].Vect())))
        
    boost=-1.0*LVStore['HHHLV'].BoostVector()
    LVStore_hhhFrame=getBoostedLVs(LVStore,boost)
    varDict["CosThetaH1_hhhF"]    =  abs(LVStore_hhhFrame['H1LV'].CosTheta())
    varDict["HH4bCosTheta_hhhF"]    =  abs(LVStore_hhhFrame['HH4bLV'].CosTheta())
    varDict["HggCosTheta_hhhF"]    =  abs(LVStore_hhhFrame['HggLV'].CosTheta())
    varDict["H1bbCosTheta_hhhF"]    =  abs(LVStore_hhhFrame['H1bbLV'].CosTheta())
    varDict["H2bbCosTheta_hhhF"]    =  abs(LVStore_hhhFrame['H2bbLV'].CosTheta())
    varDict["HH4bCosThetaLeadJet_hhhF"]    =  abs(np.cos( LVStore_hhhFrame['HH4bLV'].Angle(LVStore_hhhFrame['j1LV'].Vect())))
    varDict["absCosThetaH4bHgg_hhhF"]    =  abs(np.cos( (LVStore_hhhFrame['H1bbLV'] + LVStore['H2bbLV']).Angle(LVStore_hhhFrame['HggLV'].Vect())))
    
    vals=[
        LVStore['j1LV'].Angle( LVStore['j2LV'].Vect()),
        LVStore['j1LV'].Angle( LVStore['k1LV'].Vect()),
        LVStore['j1LV'].Angle( LVStore['k2LV'].Vect()),
    ]
    vals=abs(np.cos(vals))
    varDict["LeadJetAbsCosThetaMax"]  = max(vals)
    varDict["LeadJetAbsCosThetaMin"]  = min(vals)
   
    vals=[
        LVStore['j1LV'].DeltaR( LVStore['j2LV']),
        LVStore['j1LV'].DeltaR( LVStore['k1LV']),
        LVStore['j1LV'].DeltaR( LVStore['k2LV']),
    ]
    varDict["LeadJetDrMaxWithOtherJets"]  = max(vals)
    varDict["LeadJetDrMinWithOtherJets"]  = min(vals)
     
    vals=[
        LVStore['j1LV'].Angle( LVStore['k1LV'].Vect()),
        LVStore['j1LV'].Angle( LVStore['k2LV'].Vect()),
        LVStore['j2LV'].Angle( LVStore['k1LV'].Vect()),
        LVStore['j2LV'].Angle( LVStore['k2LV'].Vect()),
    ]

    vals=abs(np.cos(vals))
    varDict["H1H2JetAbsCosThetaMax"]  = max(vals)
    varDict["H1H2JetAbsCosThetaMin"]  = min(vals)         
    
    vals=[
        LVStore['j1LV'].DeltaR( LVStore['k1LV']),
        LVStore['j1LV'].DeltaR( LVStore['k2LV']),
        LVStore['j2LV'].DeltaR( LVStore['k1LV']),
        LVStore['j2LV'].DeltaR( LVStore['k2LV']),
    ]
    varDict["H1H2JetDrMax"]  = max(vals)
    varDict["H1H2JetDrMin"]  = min(vals)
    
    vals1=[
        LVStore['j1LV'].DeltaR( LVStore['g1LV']),
        LVStore['j2LV'].DeltaR( LVStore['g1LV']),
        LVStore['k1LV'].DeltaR( LVStore['g1LV']),
        LVStore['k2LV'].DeltaR( LVStore['g1LV']),
    ]
    vals2=[
        LVStore['j1LV'].DeltaR( LVStore['g2LV']),
        LVStore['j2LV'].DeltaR( LVStore['g2LV']),
        LVStore['k1LV'].DeltaR( LVStore['g2LV']),
        LVStore['k2LV'].DeltaR( LVStore['g2LV']),
    ]
    
    v1=min(min(vals2),min(vals1))
    v2=max(min(vals2),min(vals1))
    varDict["PhoJetMinDr"]   =  v1
    varDict["PhoJetMinDrOther"]  = v2
    
    v1=max(max(vals2),max(vals1))
    v2=min(max(vals2),max(vals1))
    varDict["PhoJetMaxDr"]   = v1
    varDict["PhoJetMaxDrOther"]  = v2

    varDict['dijet1CandidatePtOverdiHiggsM'] =   LVStore['H1bbLV'].Pt() / LVStore['H1bbLV'].M()
    varDict['dijet2CandidatePtOverdiHiggsM'] =   LVStore['H2bbLV'].Pt() / LVStore['H2bbLV'].M()
    varDict["H1bbCosThetaLeadJet"] = abs(np.cos(  LVStore['H1bbLV'].Angle( LVStore['j1LV'].Vect())  )  )
    varDict['H2bbCosThetaLeadJet'] = abs(np.cos(  LVStore['H2bbLV'].Angle( LVStore['k1LV'].Vect())  )  )

    varDict["pT_4b"]  = ( LVStore['H1bbLV'] + LVStore['H2bbLV'] ).Pt()
    varDict["scalarPtSumHHH"]  = LVStore["H1LV"].Pt()+LVStore["H2LV"].Pt()+LVStore["H3LV"].Pt()
    varDict["scalarPtSum4b"]   = LVStore["j1LV"].Pt() + LVStore["j2LV"].Pt() + LVStore["k1LV"].Pt() + LVStore["k2LV"].Pt()
    varDict["scalarPtSum4b2g"] = varDict["scalarPtSum4b"]  + LVStore["g1LV"].Pt() + LVStore["g2LV"].Pt()
    varDict["HggTo4bAbsCosTheta"]    = abs(np.cos( (LVStore['H1bbLV']+LVStore['H2bbLV']).Angle( LVStore['H1LV'].Vect())  ) )
    varDict["H1bbToH2bbAbsCosTheta"] = abs(np.cos(  LVStore['H1bbLV'].Angle( LVStore['H2bbLV'].Vect())  )  )
    
    if eTree:
        varDict['h1leadJ_deepjetScore']    =  getattr(eTree,'jet_'+str(quad['fgg_idxs'][0])+'_deepJetScore') 
        varDict['h1subleadJ_deepjetScore'] =  getattr(eTree,'jet_'+str(quad['fgg_idxs'][1])+'_deepJetScore') 
        varDict['h2leadJ_deepjetScore']    =  getattr(eTree,'jet_'+str(quad['fgg_idxs'][2])+'_deepJetScore') 
        varDict['h2subleadJ_deepjetScore'] =  getattr(eTree,'jet_'+str(quad['fgg_idxs'][3])+'_deepJetScore') 
    else:    
        varDict['h1leadJ_deepjetScore']    =  -0.2
        varDict['h1subleadJ_deepjetScore'] =  -0.2
        varDict['h2leadJ_deepjetScore']    =  -0.2
        varDict['h2subleadJ_deepjetScore'] =  -0.2

    if quad:
        j1_res = getattr(eTree,'jet_'+str(quad['fgg_idxs'][0])+'_bJetRegRes') ; 
        j2_res = getattr(eTree,'jet_'+str(quad['fgg_idxs'][1])+'_bJetRegRes') ; 
        k1_res = getattr(eTree,'jet_'+str(quad['fgg_idxs'][2])+'_bJetRegRes') ; 
        k2_res = getattr(eTree,'jet_'+str(quad['fgg_idxs'][3])+'_bJetRegRes') ; 
    else:
        j1_res = 1.0  
        j2_res = 1.0  
        k1_res = 1.0  
        k2_res = 1.0  
        
    varDict['h1_dijetSigmaMOverM'] = getSigmaMOverM( LVStore['j1LV'],j1_res , LVStore['j2LV'] , j2_res )
    varDict['h2_dijetSigmaMOverM'] = getSigmaMOverM( LVStore['k1LV'],k1_res , LVStore['k2LV'] , k2_res )
    
    varDict['triHiggs_mass']= LVStore['HHHLV'].M()
    varDict['trihiggs_pt']= LVStore['HHHLV'].Pt()
    if hasattr(eTree,'ttH_MET'):
        varDict['ttH_MET'] = eTree.ttH_MET
    else:
        varDict['ttH_MET'] = 0.0
    return varDict

def fillOtherDerivedVariables(histStore,varDict,wei):
    
    for ky in histStore["miscVars"]:
        if ky not in varDict:
            warnings.warn(ky+" : variable  not defined in varDict ! ",category=UserWarning)
            continue
        histStore["miscVars"][ky].Fill( varDict[ky] , wei)


def addCandidateVars(histStore,candList):
    for ky in candList :
        histStore[ky]={}
        histStore[ky]['nEvents']=ROOT.TH1F("nEvents","",4,-0.5,3.5)
        histStore[ky]["nEvents"].SetCanExtend(ROOT.TH1.kAllAxes);
        histStore[ky]["nEvents"].SetOption("HIST text");

        histStore[ky]['weights']=ROOT.TH1F("weights","",1000,-0.05,0.05)

        histStore[ky]['massHHH']=ROOT.TH1F("mass_hhh","",1500,0.0,1500)
        
        histStore[ky]['massHgg']=ROOT.TH1F("mass_hgg","",1000,0.0,500)
        histStore[ky]['massHbb1']=ROOT.TH1F("mass_h1bb","",2000,0.0,2000)
        histStore[ky]['massHbb2']=ROOT.TH1F("mass_h2bb","",2000,0.0,2000)
        histStore[ky]['h1MassVsh2Mass']=ROOT.TH2F("h1MassVsh2Mass","",1000,0.0,1000,1000,0.0,1000)
        
        histStore[ky]['hhh']={}
        histStore[ky]['hhh']['mass']=ROOT.TH1F("mass","",1500,0.0,1500)
        histStore[ky]['hhh']['pT']=ROOT.TH1F("pT","",1000,0.0,1000)
        histStore[ky]['hhh']['eta']=ROOT.TH1F("eta","",60,-3.0,3.0)
        histStore[ky]['hhh']['phi']=ROOT.TH1F("phi","",64,-3.2,3.2)
        
        histStore[ky]['misc']={}
        histStore[ky]['misc']['sumTheta']=ROOT.TH1F("sumTheta","sumTheta",100,0.0,10.)
        histStore[ky]['misc']["scalarPtSumHHH"]  = ROOT.TH1F("scalarPtSumHHH","",1000,0.0,1000) 
        histStore[ky]['misc']["scalarPtSum4b"]   = ROOT.TH1F("scalarPtSum4b","",1000,0.0,1000)
        histStore[ky]['misc']["scalarPtSum4b2g"] = ROOT.TH1F("scalarPtSum4b2g","",1000,0.0,1000)
        
        for hh in ['h1','h2']:
            histStore[ky][hh]={}
            for bjet in ["lead","sublead"]:
                histStore[ky][hh][bjet]={}
                histStore[ky][hh][bjet]['pT'] =ROOT.TH1F("pT","",500,0.0,2000)
                histStore[ky][hh][bjet]['eta']=ROOT.TH1F("eta","",60,-3.0,3.0)
                histStore[ky][hh][bjet]['phi']=ROOT.TH1F("phi","",62,-3.14,3.14)
                histStore[ky][hh][bjet]['mass'] =ROOT.TH1F("mass","",500,0.0,2000)
                histStore[ky][hh][bjet]['deepCSVScore']=ROOT.TH1F("deepCSVScore","",120,-3.0,3.0)
                histStore[ky][hh][bjet]['deepJetScore']=ROOT.TH1F("deepJetScore","",40,-0.5,1.5)
                histStore[ky][hh][bjet]['puJetIdMVA']=ROOT.TH1F("puJetIdMVA","",60,-1.5,1.5)
                histStore[ky][hh][bjet]['QGL']=ROOT.TH1F("QGL","",40,-0.5,1.5)
        
        histStore[ky]['hgg']={}
        histStore[ky]['hgg']['mass']=ROOT.TH1F("mass","",1000,0.0,500)
        histStore[ky]['hgg']['pT']=ROOT.TH1F("pT","",500,0.0,500)
        histStore[ky]['hgg']['eta']=ROOT.TH1F("eta","",60,-3.0,3.0)
        histStore[ky]['hgg']['phi']=ROOT.TH1F("phi","",62,-3.14,3.14)
        histStore[ky]['hgg']['pT_BB']=ROOT.TH1F("pT_BB","",1500,0.0,1500)
        histStore[ky]['hgg']['pT_BE']=ROOT.TH1F("pT_BE","",1500,0.0,1500)
        histStore[ky]['hgg']['pT_EB']=ROOT.TH1F("pT_EB","",1500,0.0,1500)
        histStore[ky]['hgg']['pT_EE']=ROOT.TH1F("pT_EE","",1500,0.0,1500)
        addFlasggVars(histStore[ky])
        addKinematicVars(histStore[ky])

@ignore_warning(RuntimeWarning)
def addObjectDeltaRValues(histStore): 
    histStore["objDelta"]={'dr':{} ,'dEta':{},'dPhi':{}}
    objs=['g1LV','g2LV','j1LV','j2LV','k1LV','k2LV']
    for i in range(len(objs)):
        for j in range(len(objs)):
            if i>=j: continue
            diObjTag = objs[i]+'To'+objs[j]
            histStore['objDelta']['dr'][diObjTag+'_DR']    = ROOT.TH1F(diObjTag + "_DeltaR","",50,0.0,5.0) 
            histStore['objDelta']['dEta'][diObjTag+'_dEta']  = ROOT.TH1F(diObjTag + "_dEta","",50,0.0,5.0) 
            histStore['objDelta']['dPhi'][diObjTag+'_dPhi']  = ROOT.TH1F(diObjTag + "_dPhi","",32,0.0,3.2)

def fillObjectDeltaRValuesFromLVStore(histStore,LVStore ,wei=1):
    
    objs=['g1LV','g2LV','j1LV','j2LV','k1LV','k2LV']
    for i in range(len(objs)):
        for j in range(len(objs)):
            if i>=j: continue
            diObjTag = objs[i]+'To'+objs[j]
            histStore['objDelta']['dr'][diObjTag+'_DR'].Fill( LVStore[objs[i]].DeltaR(LVStore[objs[j]]) ,  wei )
            histStore['objDelta']['dEta'][diObjTag+'_dEta'].Fill( abs(LVStore[objs[i]].Eta()-LVStore[objs[j]].Eta()) ,  wei )
            histStore['objDelta']['dPhi'][diObjTag+'_dPhi'].Fill( LVStore[objs[i]].DeltaPhi(LVStore[objs[j]]) ,  wei )
    
    
    
@ignore_warning(RuntimeWarning)
def addKinematicVars(histStore,framesToProbe=['CMS_RestFrame','HHH_RestFrame','Hgg_RestFrame']): 
    histStore['kinematicVars']={}
#    framesToProbe=['CMS_RestFrame',
#                   'HHH_RestFrame',
#                   'H1H2_RestFrame',
#                   'H1H3_RestFrame',
#                   'H2H3_RestFrame',
#                   'Hgg_RestFrame',
#                   'H1bb_RestFrame',
#                   'H2bb_RestFrame',
#                    ]
    
    for frame in framesToProbe:
        histStore['kinematicVars'][frame]={}
        
        histStore['kinematicVars'][frame]['H1Pt'] = ROOT.TH1F("H1Pt","",250,0.0,2500)
        histStore['kinematicVars'][frame]['H2Pt'] = ROOT.TH1F("H2Pt","",250,0.0,2500)
        histStore['kinematicVars'][frame]['H3Pt'] = ROOT.TH1F("H3Pt","",250,0.0,2500)
    
        histStore['kinematicVars'][frame]['H1Eta'] = ROOT.TH1F("H1Eta","",100,-5.0,5.0)
        histStore['kinematicVars'][frame]['H2Eta'] = ROOT.TH1F("H2Eta","",100,-5.0,5.0)
        histStore['kinematicVars'][frame]['H3Eta'] = ROOT.TH1F("H3Eta","",100,-5.0,5.0)
    
        histStore['kinematicVars'][frame]['H1y'] = ROOT.TH1F("H1y","",100,-5.0,5.0)
        histStore['kinematicVars'][frame]['H2y'] = ROOT.TH1F("H2y","",100,-5.0,5.0)
        histStore['kinematicVars'][frame]['H3y'] = ROOT.TH1F("H3y","",100,-5.0,5.0)
    
        histStore['kinematicVars'][frame]['H1Phi'] = ROOT.TH1F("H1Phi","",100,-3.14,3.14)
        histStore['kinematicVars'][frame]['H2Phi'] = ROOT.TH1F("H2Phi","",100,-3.14,3.14)
        histStore['kinematicVars'][frame]['H3Phi'] = ROOT.TH1F("H3Phi","",100,-3.14,3.14)

        histStore['kinematicVars'][frame]['H1H2VsH3_CosTheta']  =  ROOT.TH1F("H1H2VsH3_CosTheta","",22,-1.1,1.1)
        histStore['kinematicVars'][frame]['H1H3VsH2_CosTheta']  =  ROOT.TH1F("H1H3VsH2_CosTheta","",22,-1.1,1.1)
        histStore['kinematicVars'][frame]['H2H3VsH1_CosTheta']  =  ROOT.TH1F("H2H3VsH1_CosTheta","",22,-1.1,1.1)
        
        histStore['kinematicVars'][frame]['H1VsH2_CosTheta']  =  ROOT.TH1F("H1VsH2_CosTheta","",22,-1.1,1.1)
        histStore['kinematicVars'][frame]['H1VsH3_CosTheta']  =  ROOT.TH1F("H1VsH3_CosTheta","",22,-1.1,1.1)
        histStore['kinematicVars'][frame]['H2VsH3_CosTheta']  =  ROOT.TH1F("H2VsH3_CosTheta","",22,-1.1,1.1)
        
        histStore['kinematicVars'][frame]['H1H2_DR']  = ROOT.TH1F("H1H2DeltaR","",50,0.0,5.0) 
        histStore['kinematicVars'][frame]['H1H3_DR']  = ROOT.TH1F("H1H3DeltaR","",50,0.0,5.0) 
        histStore['kinematicVars'][frame]['H2H3_DR']  = ROOT.TH1F("H2H3DeltaR","",50,0.0,5.0) 

        histStore['kinematicVars'][frame]['H1H2_dPhi']  =  ROOT.TH1F("H1H2_dPhi","",32,0.0,3.2)
        histStore['kinematicVars'][frame]['H1H3_dPhi']  =  ROOT.TH1F("H1H3_dPhi","",32,0.0,3.2)
        histStore['kinematicVars'][frame]['H2H3_dPhi']  =  ROOT.TH1F("H2H3_dPhi","",32,0.0,3.2)
         
        histStore['kinematicVars'][frame]['H1H2_dEta']  = ROOT.TH1F("H1H2_dEta","",50,0.0,5.0) 
        histStore['kinematicVars'][frame]['H1H3_dEta']  = ROOT.TH1F("H1H2_dEta","",50,0.0,5.0) 
        histStore['kinematicVars'][frame]['H2H3_dEta']  = ROOT.TH1F("H1H2_dEta","",50,0.0,5.0) 

        histStore['kinematicVars'][frame]['H1H2_area']  = ROOT.TH1F("H1H2_area","",200,0.0,20.0) 
        histStore['kinematicVars'][frame]['H1H3_area']  = ROOT.TH1F("H1H3_area","",200,0.0,20.0) 
        histStore['kinematicVars'][frame]['H2H3_area']  = ROOT.TH1F("H2H3_area","",200,0.0,20.0) 
        
        histStore['kinematicVars'][frame]['H1H2H3_volume']  =  ROOT.TH1F("H1H2H3_volume"  ,""  ,40,0.0,2.0)
        histStore['kinematicVars'][frame]['L1L2L3_volume']  =  ROOT.TH1F("L1L2L3_volume"  ,""  ,40,0.0,2.0)
        histStore['kinematicVars'][frame]['sL1sL2sL3_volume']  =   ROOT.TH1F("sL1sL2sL3_volume","",40,0.0,2.0)
        
        addOtherDerivedVariables(histStore['kinematicVars'][frame])
    
def fillKinematicVarsFromLV(LVStore,histStore,wei=1 ,framesToDo=['CMS_RestFrame','HHH_RestFrame','Hgg_RestFrame']):
    
    if 'CMS_RestFrame' in framesToDo:
       fillKinematicVars(histStore['CMS_RestFrame'],LVStore)    

    if 'HHH_RestFrame' in framesToDo:
        boost=-1.0*LVStore['HHHLV'].BoostVector()
        boostedLV=getBoostedLVs(LVStore,boost)
        fillKinematicVars(histStore['HHH_RestFrame'],boostedLV , wei )    
    
    if 'H1H2_RestFrame' in framesToDo:
        boost=-1.0*(LVStore['H1LV'] +LVStore['H2LV']).BoostVector()
        boostedLV=getBoostedLVs(LVStore,boost)
        fillKinematicVars(histStore['H1H2_RestFrame'],boostedLV , wei )    
    
    if 'H1H3_RestFrame' in framesToDo:
        boost=-1.0*(LVStore['H1LV']+LVStore['H2LV']).BoostVector()
        boostedLV=getBoostedLVs(LVStore,boost)
        fillKinematicVars(histStore['H1H3_RestFrame'],boostedLV , wei )    
        
    if 'H2H3_RestFrame' in framesToDo:
        boost=-1.0*(LVStore['H2LV']+LVStore['H3LV']).BoostVector()
        boostedLV=getBoostedLVs(LVStore,boost)
        fillKinematicVars(histStore['H2H3_RestFrame'],boostedLV , wei )    
        
    if 'Hgg_RestFrame' in framesToDo:
        boost=-1.0*LVStore['HggLV'].BoostVector()
        boostedLV=getBoostedLVs(LVStore,boost)
        fillKinematicVars(histStore['Hgg_RestFrame'],boostedLV , wei )    
        
    if 'H1bb_RestFrame' in framesToDo:
        boost=-1.0*LVStore['H1bbLV'].BoostVector()
        boostedLV=getBoostedLVs(LVStore,boost)
        fillKinematicVars(histStore['H1bb_RestFrame'],boostedLV , wei )    

    if 'H2bb_RestFrame' in framesToDo:
        boost=-1.0*LVStore['H2bbLV'].BoostVector()
        boostedLV=getBoostedLVs(LVStore,boost)
        fillKinematicVars(histStore['H2bb_RestFrame'],boostedLV , wei )    

    return


@ignore_warning(RuntimeWarning)
def getRecoHistos():
    histStore={}
    histStore['events']={}
    histStore['events']['nEvents']=ROOT.TH1F("nEvents","",1,-0.5,0.5)
    histStore['events']["nEvents"].SetCanExtend(ROOT.TH1.kAllAxes);
    histStore['events']['allFailIdsEvents']=ROOT.TH1F("allFailedEvents","",1,-0.5,0.5)
    histStore['events']["allFailIdsEvents"].SetCanExtend(ROOT.TH1.kAllAxes);
    histStore['events']['allFailIdsWeights']=ROOT.TH1F("allFailedWeights","",1,-0.5,0.5)
    histStore['events']["allFailIdsWeights"].SetCanExtend(ROOT.TH1.kAllAxes);
    
    histStore['events']['nHiggsMatch']=ROOT.TH1F("nHiggsMatch","",40,-0.5,39.5)
    histStore['events']["nHiggsMatch"].SetCanExtend(ROOT.TH1.kAllAxes);
    
    histStore['diPhotons']={}
    histStore['diPhotons']['mass']=ROOT.TH1F("mass","",1000,0.0,500)
    histStore['diPhotons']['pT']=ROOT.TH1F("pT","",500,0.0,500)
    histStore['diPhotons']['eta']=ROOT.TH1F("eta","",60,-3.0,3.0)
    histStore['diPhotons']['y']=ROOT.TH1F("y","",60,-3.0,3.0)
    histStore['diPhotons']['phi']=ROOT.TH1F("phi","",62,-3.14,3.14)
    histStore['diPhotons']['g1g2_Dr'] =ROOT.TH1F("g1g2_Dr","",50,0.0,5)
    histStore['diPhotons']['gamma1_pT'] =ROOT.TH1F("gamma1_pT","",500,0.0,500)
    histStore['diPhotons']['gamma1_eta']=ROOT.TH1F("gamma1_eta","",60,-3.0,3.0)
    histStore['diPhotons']['gamma1_phi']=ROOT.TH1F("gamma1_phi","",62,-3.14,3.14)
    histStore['diPhotons']['gamma2_pT'] =ROOT.TH1F("gamma2_pT","",500,0.0,500)
    histStore['diPhotons']['gamma2_eta']=ROOT.TH1F("gamma2_eta","",60,-3.0,3.0)
    histStore['diPhotons']['gamma2_phi']=ROOT.TH1F("gamma2_phi","",62,-3.14,3.14)  
    histStore['diPhotons']['detIdx']=ROOT.TH1F("detIdx","",10,-0.5,9.5)    
    
    histStore['h1Tobb']={}
    histStore['h1Tobb']['mass']=ROOT.TH1F("mass","",1000,0.0,2000)
    histStore['h1Tobb']['mass_preReg']=ROOT.TH1F("mass_preReg","",1000,0.0,2000)
    histStore['h1Tobb']['pT']=ROOT.TH1F("pT","",500,0.0,2000)
    histStore['h1Tobb']['eta']=ROOT.TH1F("eta","",60,-3.0,3.0)
    histStore['h1Tobb']['y']=ROOT.TH1F("y","",60,-3.0,3.0)
    histStore['h1Tobb']['phi']=ROOT.TH1F("phi","",62,-3.14,3.14)
    histStore['h1Tobb']['b1b2_Dr'] =ROOT.TH1F("g1g2_Dr","",50,0.0,5)
    histStore['h1Tobb']['bJet1_pT'] =ROOT.TH1F("bJet1_pT","",500,0.0,2000)
    histStore['h1Tobb']['bJet1_pT_preReg'] =ROOT.TH1F("bJet1_pT_preReg","",500,0.0,2000)
    histStore['h1Tobb']['bJet1_eta']=ROOT.TH1F("bJet1_eta","",60,-3.0,3.0)
    histStore['h1Tobb']['bJet1_phi']=ROOT.TH1F("bJet1_phi","",62,-3.14,3.14)
    histStore['h1Tobb']['bJet1_deepCSVScore']=ROOT.TH1F("bJet1_deepCSVScore","",120,-3.0,3.0)
    histStore['h1Tobb']['bJet1_deepJetScore']=ROOT.TH1F("bJet1_deepJetScore","",40,-0.5,1.5)
    histStore['h1Tobb']['bJet1_hFlavour']=ROOT.TH1F("bJet1_hFlavour","",9,-0.5,8.0)
    histStore['h1Tobb']['bJet1_pFlavour']=ROOT.TH1F("bJet1_pFlavour","",21,-10.5,10.5)
    histStore['h1Tobb']['bJet2_pT'] =ROOT.TH1F("bJet2_pT","",500,0.0,2000)
    histStore['h1Tobb']['bJet2_pT_preReg'] =ROOT.TH1F("bJet2_pT_preReg","",500,0.0,2000)
    histStore['h1Tobb']['bJet2_eta']=ROOT.TH1F("bJet2_eta","",60,-3.0,3.0)
    histStore['h1Tobb']['bJet2_phi']=ROOT.TH1F("bJet2_phi","",62,-3.14,3.14)
    histStore['h1Tobb']['bJet2_deepCSVScore']=ROOT.TH1F("bJet2_deepCSVScore","",120,-3.0,3.0)
    histStore['h1Tobb']['bJet2_deepJetScore']=ROOT.TH1F("bJet2_deepJetScore","",40,-0.5,1.5)
    histStore['h1Tobb']['bJet2_hFlavour']=ROOT.TH1F("bJet2_hFlavour","",9,-0.5,8.0)
    histStore['h1Tobb']['bJet2_pFlavour']=ROOT.TH1F("bJet2_pFlavour","",21,-10.5,10.5)
    
    histStore['h2Tobb']={}
    histStore['h2Tobb']['mass']=ROOT.TH1F("mass","",1000,0.0,2000)
    histStore['h2Tobb']['mass_preReg']=ROOT.TH1F("mass_preReg","",1000,0.0,2000)
    histStore['h2Tobb']['pT']=ROOT.TH1F("pT","",500,0.0,2000)
    histStore['h2Tobb']['eta']=ROOT.TH1F("eta","",60,-3.0,3.0)
    histStore['h2Tobb']['y']=ROOT.TH1F("y","",60,-3.0,3.0)
    histStore['h2Tobb']['phi']=ROOT.TH1F("phi","",62,-3.14,3.14)
    histStore['h2Tobb']['b1b2_Dr'] =ROOT.TH1F("b1b2_Dr","",50,0.0,5)
    histStore['h2Tobb']['bJet1_pT'] =ROOT.TH1F("bJet1_pT","",500,0.0,2000)
    histStore['h2Tobb']['bJet1_pT_preReg'] =ROOT.TH1F("bJet1_pT_preReg","",500,0.0,2000)
    histStore['h2Tobb']['bJet1_eta']=ROOT.TH1F("bJet1_eta","",60,-3.0,3.0)
    histStore['h2Tobb']['bJet1_phi']=ROOT.TH1F("bJet1_phi","",62,-3.14,3.14)
    histStore['h2Tobb']['bJet1_deepCSVScore']=ROOT.TH1F("bJet1_deepCSVScore","",120,-3.0,3.0)
    histStore['h2Tobb']['bJet1_deepJetScore']=ROOT.TH1F("bJet1_deepJetScore","",40,-0.5,1.5)
    histStore['h2Tobb']['bJet1_hFlavour']=ROOT.TH1F("bJet1_hFlavour","",9,-0.5,8.0)
    histStore['h2Tobb']['bJet1_pFlavour']=ROOT.TH1F("bJet1_pFlavour","",21,-10.5,10.5)
    histStore['h2Tobb']['bJet2_pT'] =ROOT.TH1F("bJet2_pT","",500,0.0,2000)
    histStore['h2Tobb']['bJet2_pT_preReg'] =ROOT.TH1F("bJet2_pT_preReg","",500,0.0,2000)
    histStore['h2Tobb']['bJet2_eta']=ROOT.TH1F("bJet2_eta","",60,-3.0,3.0)
    histStore['h2Tobb']['bJet2_phi']=ROOT.TH1F("bJet2_phi","",62,-3.14,3.14)
    histStore['h2Tobb']['bJet2_deepCSVScore']=ROOT.TH1F("bJet2_deepCSVScore","",120,-3.0,3.0)
    histStore['h2Tobb']['bJet2_deepJetScore']=ROOT.TH1F("bJet2_deepJetScore","",40,-0.5,1.5)
    histStore['h2Tobb']['bJet2_hFlavour']=ROOT.TH1F("bJet2_hFlavour","",9,-0.5,8.0)
    histStore['h2Tobb']['bJet2_pFlavour']=ROOT.TH1F("bJet2_pFlavour","",21,-10.5,10.5)
    

    histStore['hhTo4b']={}
    histStore['hhTo4b']['pT']=ROOT.TH1F("pT","",500,0.0,2000)
    histStore['hhTo4b']['mass']=ROOT.TH1F("mass","",1000,0.0,2000)
    histStore['hhTo4b']['mass_preReg']=ROOT.TH1F("mass_preReg","",1000,0.0,2000)
    histStore['hhTo4b']['eta']=ROOT.TH1F("eta","",60,-3.0,3.0)
    histStore['hhTo4b']['y']=ROOT.TH1F("y","",60,-3.0,3.0)
    histStore['hhTo4b']['phi']=ROOT.TH1F("phi","",62,-3.14,3.14)
    
    histStore['hhhTo4b2gamma']={}
    histStore['hhhTo4b2gamma']['mass']=ROOT.TH1F("mass","",1000,0.0,2000)
    histStore['hhhTo4b2gamma']['h1MassVsh2Mass']=ROOT.TH2F("h1MassVsh2Mass","",1000,0.0,1000,1000,0.0,1000)
    histStore['hhhTo4b2gamma']['massXbar']=ROOT.TH1F("massXbar","",1000,0.0,2000)
    histStore['hhhTo4b2gamma']['mass_preReg']=ROOT.TH1F("mass_preReg","",1000,0.0,2000)
    histStore['hhhTo4b2gamma']['pT']=ROOT.TH1F("pT","",500,0.0,2000)
    histStore['hhhTo4b2gamma']['eta']=ROOT.TH1F("eta","",60,-3.0,3.0)
    histStore['hhhTo4b2gamma']['y']=ROOT.TH1F("y","",60,-3.0,3.0)
    histStore['hhhTo4b2gamma']['phi']=ROOT.TH1F("phi","",62,-3.14,3.14)
    histStore['hhhTo4b2gamma']['h1MassVsh2Mass_preReg']=ROOT.TH2F("h1MassVsh2Mass_preReg","",1000,0.0,1000,1000,0.0,1000)
    histStore['hhhTo4b2gamma']['eventJetMultiplicity']=ROOT.TH1F("eventJetMultiplicity","",40,-0.5,39.5)

    histStore['hhhTo4b2gamma']['massX0']=ROOT.TH1F("massX0","",1000,0.0,2000)
    histStore['hhhTo4b2gamma']['massX1']=ROOT.TH1F("massX1","",1000,0.0,2000)
    histStore['hhhTo4b2gamma']['massX2']=ROOT.TH1F("massX2","",1000,0.0,2000)
    histStore['hhhTo4b2gamma']['massX3']=ROOT.TH1F("massX3","",1000,0.0,2000)
    histStore["hhhTo4b2gamma"]["massX0X1"]=ROOT.TH2F("h1MassX0X1","h1MassX0X1",1000,0.0,1000,1000,0.0,1000)
    histStore["hhhTo4b2gamma"]["massX2X3"]=ROOT.TH2F("h1MassX2X3","h1MassX2X3",1000,0.0,1000,1000,0.0,1000)
    histStore["hhhTo4b2gamma"]["massX0X2"]=ROOT.TH2F("h1MassX0X2","h1MassX0X2",1000,0.0,1000,1000,0.0,1000)
    histStore["hhhTo4b2gamma"]["massX1X3"]=ROOT.TH2F("h1MassX1X3","h1MassX1X3",1000,0.0,1000,1000,0.0,1000)
    histStore["hhhTo4b2gamma"]["massX0X3"]=ROOT.TH2F("h1MassX0X3","h1MassX0X3",1000,0.0,1000,1000,0.0,1000)
    histStore["hhhTo4b2gamma"]["massX1X2"]=ROOT.TH2F("h1MassX1X2","h1MassX1X2",1000,0.0,1000,1000,0.0,1000)
    
    histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_tightID"]= histStore["hhhTo4b2gamma"]["h1MassVsh2Mass"].Clone()
    histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_tightID"].SetName("h1MassVsh2Mass_tightID")
    
    histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_mediumID"]= histStore["hhhTo4b2gamma"]["h1MassVsh2Mass"].Clone()
    histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_mediumID"].SetName("h1MassVsh2Mass_mediumID")
    
    histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_looseID"]= histStore["hhhTo4b2gamma"]["h1MassVsh2Mass"].Clone()
    histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_looseID"].SetName("h1MassVsh2Mass_looseID")
    
    histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_looseID"]= histStore["hhhTo4b2gamma"]["h1MassVsh2Mass"].Clone()
    histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_looseID"].SetName("h1MassVsh2Mass_looseID")
    
#    histStore["h1Tobb"]["mass_tightID" ] =histStore["h1Tobb"]["mass" ].Clone()
#    histStore["h2Tobb"]["mass_tightID" ] =histStore["h2Tobb"]["mass" ].Clone()
#    histStore["h1Tobb"]["mass_tightID" ].SetName("mass_tightID" )
#    histStore["h2Tobb"]["mass_tightID" ].SetName("mass_tightID" )
#    
#    histStore["h1Tobb"]["mass_mediumID" ] =histStore["h1Tobb"]["mass" ].Clone()
#    histStore["h2Tobb"]["mass_mediumID" ] =histStore["h2Tobb"]["mass" ].Clone()
#    histStore["h1Tobb"]["mass_mediumID" ].SetName("mass_mediumID" )
#    histStore["h2Tobb"]["mass_mediumID" ].SetName("mass_mediumID")
#    
#    histStore["h1Tobb"]["mass_looseID" ] =histStore["h1Tobb"]["mass" ].Clone()
#    histStore["h2Tobb"]["mass_looseID" ] =histStore["h2Tobb"]["mass" ].Clone()
#    histStore["h1Tobb"]["mass_looseID" ].SetName("mass_looseID" )
#    histStore["h2Tobb"]["mass_looseID" ].SetName("mass_looseID" )



    addKinematicVars(histStore)

    return histStore

def fillBJetHistos():
    return

def getCosthetaVars(eTree,LVStore):
    drMin=1e5
    for g in ['g1LV','g2LV']:
        for j in ['j1LV','j2LV','k1LV','k2LV']:
            dr= LVStore[g].DeltaR(LVStore[j])
            if dr < drMin:
                gg=g
                jj=j
                drMin=dr
    drMin2=1e5
    for g in ['g1LV','g2LV']:
        if g==gg:
            continue
        for j in ['j1LV','j2LV','k1LV','k2LV']:
            if j==jj:
                continue
            dr= LVStore[g].DeltaR(LVStore[j])
            if dr < drMin2:
                drMin2=dr
    x=copy.deepcopy(LVStore)          
    x['j1LV'].Boost(-1.0* x['H1bbLV'].BoostVector())
    x['k1LV'].Boost(-1.0* x['H2bbLV'].BoostVector())
    x['HggLV'].Boost(-1.0* x['HHHLV'].BoostVector())
    
    return x['j1LV'].CosTheta() , x['k1LV'].CosTheta() , x['HggLV'].CosTheta() , drMin , drMin2


def getLVStoreFromTreeAndQuad( eTree ,quad):
    
    LVStore={}
    LVStore['j1LV']=ROOT.TLorentzVector();
    LVStore['j2LV']=ROOT.TLorentzVector();
    LVStore['k1LV']=ROOT.TLorentzVector();
    LVStore['k2LV']=ROOT.TLorentzVector();
    
    LVStore['g1LV']=ROOT.TLorentzVector()
    LVStore['g1LV'].SetPtEtaPhiM(eTree.leadingPhoton_pt,eTree.leadingPhoton_eta,eTree.leadingPhoton_phi,0.0)

    LVStore['g2LV']=ROOT.TLorentzVector()
    LVStore['g2LV'].SetPtEtaPhiM(eTree.subleadingPhoton_pt,eTree.subleadingPhoton_eta,eTree.subleadingPhoton_phi,0.0)
    
    kys=['j1LV','j2LV','k1LV','k2LV']
    for i in range(4):
        ix=quad['fgg_idxs'][i]
        pt = getattr(eTree,'jet_'+str(ix)+'_pt')  
        eta= getattr(eTree,'jet_'+str(ix)+'_eta') 
        phi= getattr(eTree,'jet_'+str(ix)+'_phi') 
        mass= getattr(eTree,'jet_'+str(ix)+'_mass')  
        LVStore[kys[i]].SetPtEtaPhiM(pt,eta,phi,mass)

    LVStore['HggLV'] =ROOT.TLorentzVector()
    LVStore['HggLV'].SetPtEtaPhiM( eTree.diphoton_pt , eTree.diphoton_eta , eTree.diphoton_phi , eTree.CMS_hgg_mass )
    LVStore['H1bbLV']= LVStore['j1LV'] + LVStore['j2LV']
    LVStore['H2bbLV']= LVStore['k1LV'] + LVStore['k2LV']
    
    LVStore['HH4bLV'] =LVStore['H1bbLV']+LVStore['H2bbLV']
    LVStore['HHHLV']  =LVStore['HggLV']+LVStore['H1bbLV']+LVStore['H2bbLV']
    
    LVStore['H1LV']=LVStore['HggLV']
    LVStore['H2LV']=LVStore['HggLV']
    LVStore['H3LV']=LVStore['HggLV']
    
    if LVStore['HggLV'].Pt() < LVStore['H1bbLV'].Pt():
        LVStore['H1LV']=LVStore['H1bbLV']
        if LVStore['HggLV'].Pt() < LVStore['H2bbLV'].Pt():
            LVStore['H2LV']=LVStore['H2bbLV']
        else:
            LVStore['H3LV']=LVStore['H2bbLV']
    else:
        LVStore['H2LV']=LVStore['H1bbLV']
        LVStore['H3LV']=LVStore['H2bbLV']
 
    return LVStore




def getLVStore( eTree ):
    
    LVStore={}
    LVStore['j1LV']=ROOT.TLorentzVector();
    LVStore['j2LV']=ROOT.TLorentzVector();
    LVStore['k1LV']=ROOT.TLorentzVector();
    LVStore['k2LV']=ROOT.TLorentzVector();
    
    LVStore['g1LV']=ROOT.TLorentzVector()
    LVStore['g1LV'].SetPtEtaPhiM(eTree.leadingPhoton_pt,eTree.leadingPhoton_eta,eTree.leadingPhoton_phi,0.0)

    LVStore['g2LV']=ROOT.TLorentzVector()
    LVStore['g2LV'].SetPtEtaPhiM(eTree.subleadingPhoton_pt,eTree.subleadingPhoton_eta,eTree.subleadingPhoton_phi,0.0)
    
    LVStore['j1LV'].SetPtEtaPhiM(eTree.h1LeadingJet_pt,eTree.h1LeadingJet_eta,eTree.h1LeadingJet_phi,eTree.h1LeadingJet_mass);
    LVStore['j2LV'].SetPtEtaPhiM(eTree.h1SubleadingJet_pt,eTree.h1SubleadingJet_eta,eTree.h1SubleadingJet_phi,eTree.h1SubleadingJet_mass);
    LVStore['k1LV'].SetPtEtaPhiM(eTree.h2LeadingJet_pt,eTree.h2LeadingJet_eta,eTree.h2LeadingJet_phi,eTree.h2LeadingJet_mass);
    LVStore['k2LV'].SetPtEtaPhiM(eTree.h2SubleadingJet_pt,eTree.h2SubleadingJet_eta,eTree.h2SubleadingJet_phi,eTree.h2SubleadingJet_mass);
    
    LVStore['HggLV'] =ROOT.TLorentzVector()
    LVStore['HggLV'].SetPtEtaPhiM( eTree.diphoton_pt , eTree.diphoton_eta , eTree.diphoton_phi , eTree.CMS_hgg_mass )
    LVStore['H1bbLV']= LVStore['j1LV'] + LVStore['j2LV']
    LVStore['H2bbLV']= LVStore['k1LV'] + LVStore['k2LV']
    
    LVStore['HHHLV']=LVStore['HggLV']+LVStore['H1bbLV']+LVStore['H2bbLV']
    
    LVStore['H1LV']=LVStore['HggLV']
    LVStore['H2LV']=LVStore['HggLV']
    LVStore['H3LV']=LVStore['HggLV']
    
    if LVStore['HggLV'].Pt() < LVStore['H1bbLV'].Pt():
        LVStore['H1LV']=LVStore['H1bbLV']
        if LVStore['HggLV'].Pt() < LVStore['H2bbLV'].Pt():
            LVStore['H2LV']=LVStore['H2bbLV']
        else:
            LVStore['H3LV']=LVStore['H2bbLV']
    else:
        LVStore['H2LV']=LVStore['H1bbLV']
        LVStore['H3LV']=LVStore['H2bbLV']
 
    return LVStore

def getDijetResolutions(eTree , LVStore):
 
    dijet1SigmaMOverM = np.sqrt(pow(eTree.h1LeadingJet_bRegNNResolution,2)*pow( LVStore['j1LV'].M()*LVStore['j1LV'].M() +  LVStore['j1LV'].Dot(LVStore['j2LV']),2)  +
                                         pow(eTree.h1SubleadingJet_bRegNNResolution,2)*pow( LVStore['j2LV'].M()*LVStore['j2LV'].M()+ LVStore['j1LV'].Dot(LVStore['j2LV']),2))/pow(eTree.M1jj,2)
    diket2SigmaMOverM = np.sqrt(pow(eTree.h2LeadingJet_bRegNNResolution,2)*pow( LVStore['k1LV'].M()*LVStore['k1LV'].M() +  LVStore['k1LV'].Dot(LVStore['k2LV']),2)+
                                         pow(eTree.h2SubleadingJet_bRegNNResolution,2)*pow( LVStore['k2LV'].M()*LVStore['k2LV'].M()+ LVStore['k1LV'].Dot(LVStore['k2LV']),2))/pow(eTree.M2jj,2)
    
    return dijet1SigmaMOverM,diket2SigmaMOverM;

def fillFlashggVars(histStore,eTree,LVStore,wei=None):
    if wei ==None:
        wei=eTree.weight   
    histStore["flashggVars"]['HHHCosThetaH1'].Fill( eTree.absCosThetaStar_CS , wei )
    
    a,b,c,d,e=getCosthetaVars(eTree,LVStore)
    
    histStore["flashggVars"]['HHHCosThetaHgg'].Fill( c , wei ) #TODO
    histStore["flashggVars"]['HggCosThetaLeadGamma'].Fill( eTree.absCosTheta_gg , wei )
    
    histStore["flashggVars"]['pTgg_overMgg'].Fill( eTree.diphotonCandidatePtOverdiHiggsM , wei )
    
    histStore["flashggVars"]['pTleadG_overMgg'].Fill( eTree.leadingPhoton_pt / eTree.CMS_hgg_mass , wei )

    histStore["flashggVars"]['pTsubleadG_overMgg'].Fill( eTree.subleadingPhoton_pt / eTree.CMS_hgg_mass , wei )
    histStore["flashggVars"]['customLeadingPhotonIDMVA'].Fill( eTree.customLeadingPhotonIDMVA , wei )
    histStore["flashggVars"]['customSubLeadingPhotonIDMVA'].Fill( eTree.customSubLeadingPhotonIDMVA , wei )
    
    histStore["flashggVars"]['h1leadG_SigOverE'].Fill( eTree.leadingPhotonSigOverE , wei )
    histStore["flashggVars"]['h1subleadG_SigOverE'].Fill( eTree.subleadingPhotonSigOverE , wei )
    histStore["flashggVars"]['hgg_SigMOverM'].Fill( eTree.sigmaMOverM , wei )
    histStore["flashggVars"]['hgg_SigMOverMDecorrelated'].Fill( eTree.sigmaMOverMDecorr , wei )
    
    histStore["flashggVars"]['rho'].Fill( eTree.rho , wei )



def fillCandidateHistograms(histStore,eTree,wei):
    
    histStore["weights"].Fill(wei)
    histStore["nEvents"].Fill("total",1)
    histStore["nEvents"].Fill("totalWeighted",wei)
    
    histStore['massHHH'].Fill(eTree.triHiggs_mass,wei)

    histStore['massHgg'].Fill(eTree.CMS_hgg_mass,wei)
    histStore['massHbb1'].Fill(eTree.M1jj,wei)
    histStore['massHbb2'].Fill(eTree.M2jj,wei)
    histStore['h1MassVsh2Mass'].Fill(eTree.M1jj,eTree.M2jj,wei)

    
    histStore['hhh']['pT'].Fill( eTree.triHiggs_pt , wei )
    histStore['hhh']['eta'].Fill( eTree.triHiggs_eta , wei )
    histStore['hhh']['phi'].Fill( eTree.triHiggs_phi , wei )

    LVStore=getLVStore(eTree)
    sumTheta=LVStore['H1LV'].Angle(LVStore['H2LV'].Vect()) +LVStore['H2LV'].Angle(LVStore['H3LV'].Vect()) +LVStore['H3LV'].Angle(LVStore['H1LV'].Vect())
    
    histStore['misc']['sumTheta'].Fill( sumTheta , wei )
    
    histStore['misc']["scalarPtSumHHH"].Fill(   LVStore["H1LV"].Pt()+LVStore["H2LV"].Pt()+LVStore["H3LV"].Pt() )
    histStore['misc']["scalarPtSum4b"].Fill(    LVStore["j1LV"].Pt() + LVStore["j2LV"].Pt() + LVStore["k1LV"].Pt() + LVStore["k2LV"].Pt() )
    histStore['misc']["scalarPtSum4b2g"].Fill(  LVStore["j1LV"].Pt() + LVStore["j2LV"].Pt() + LVStore["k1LV"].Pt() + LVStore["k2LV"].Pt()  + LVStore["g1LV"].Pt() + LVStore["g2LV"].Pt() )

    if abs(eTree.leadingPhoton_eta) < 1.44 and abs(eTree.subleadingPhoton_eta) < 1.44:
        histStore['hgg']['pT_BB'].Fill( eTree.diphoton_pt , wei )
    if abs(eTree.leadingPhoton_eta) < 1.44 and abs(eTree.subleadingPhoton_eta) > 1.567:
        histStore['hgg']['pT_BE'].Fill( eTree.diphoton_pt , wei )
    if abs(eTree.leadingPhoton_eta) > 1.567 and abs(eTree.subleadingPhoton_eta) < 1.44:
        histStore['hgg']['pT_EB'].Fill( eTree.diphoton_pt , wei )
    if abs(eTree.leadingPhoton_eta) > 1.567 and abs(eTree.subleadingPhoton_eta) < 1.567:
        histStore['hgg']['pT_EE'].Fill( eTree.diphoton_pt , wei )
        
    histStore['hgg']['eta'].Fill( eTree.diphoton_eta , wei )
    histStore['hgg']['phi'].Fill( eTree.diphoton_phi , wei )
    histStore['hgg']['mass'].Fill( eTree.CMS_hgg_mass , wei )
    
    histStore['h1']["lead"]['pT'].Fill( eTree.h1LeadingJet_pt , wei )
    histStore['h1']["lead"]['eta'].Fill( eTree.h1LeadingJet_eta , wei )
    histStore['h1']["lead"]['phi'].Fill( eTree.h1LeadingJet_phi , wei )
    histStore['h1']["lead"]['mass'].Fill( eTree.h1LeadingJet_mass , wei )
    histStore['h1']["lead"]['deepCSVScore'].Fill( eTree.h1LeadingJet_DeepCSV , wei )
    histStore['h1']["lead"]['deepJetScore'].Fill( eTree.h1LeadingJet_DeepFlavour , wei )
    histStore['h1']["lead"]['puJetIdMVA'].Fill( eTree.h1LeadingJet_puJetIdMVA , wei )
    histStore['h1']["lead"]['QGL'].Fill( eTree.h1LeadingJet_QGL , wei )
    
    histStore['h1']["sublead"]['pT'].Fill( eTree.h1SubleadingJet_pt , wei )
    histStore['h1']["sublead"]['eta'].Fill( eTree.h1SubleadingJet_eta , wei )
    histStore['h1']["sublead"]['phi'].Fill( eTree.h1SubleadingJet_phi , wei )
    histStore['h1']["sublead"]['mass'].Fill( eTree.h1SubleadingJet_mass , wei )
    histStore['h1']["sublead"]['deepCSVScore'].Fill( eTree.h1SubleadingJet_DeepCSV , wei )
    histStore['h1']["sublead"]['deepJetScore'].Fill( eTree.h1SubleadingJet_DeepFlavour , wei )
    histStore['h1']["sublead"]['puJetIdMVA'].Fill( eTree.h1SubleadingJet_puJetIdMVA , wei )
    histStore['h1']["sublead"]['QGL'].Fill( eTree.h1SubleadingJet_QGL , wei )
    
    
    histStore['h2']["lead"]['pT'].Fill( eTree.h2LeadingJet_pt , wei )
    histStore['h2']["lead"]['eta'].Fill( eTree.h2LeadingJet_eta , wei )
    histStore['h2']["lead"]['phi'].Fill( eTree.h2LeadingJet_phi , wei )
    histStore['h2']["lead"]['mass'].Fill( eTree.h2LeadingJet_mass , wei )
    histStore['h2']["lead"]['deepCSVScore'].Fill( eTree.h2LeadingJet_DeepCSV , wei )
    histStore['h2']["lead"]['deepJetScore'].Fill( eTree.h2LeadingJet_DeepFlavour , wei )
    histStore['h2']["lead"]['puJetIdMVA'].Fill( eTree.h2LeadingJet_puJetIdMVA , wei )
    histStore['h2']["lead"]['QGL'].Fill( eTree.h2LeadingJet_QGL , wei )
    
    histStore['h2']["sublead"]['pT'].Fill( eTree.h2SubleadingJet_pt , wei )
    histStore['h2']["sublead"]['eta'].Fill( eTree.h2SubleadingJet_eta , wei )
    histStore['h2']["sublead"]['phi'].Fill( eTree.h2SubleadingJet_phi , wei )
    histStore['h2']["sublead"]['mass'].Fill( eTree.h2SubleadingJet_mass , wei )
    histStore['h2']["sublead"]['deepCSVScore'].Fill( eTree.h2SubleadingJet_DeepCSV , wei )
    histStore['h2']["sublead"]['deepJetScore'].Fill( eTree.h2SubleadingJet_DeepFlavour , wei )
    histStore['h2']["sublead"]['puJetIdMVA'].Fill( eTree.h2SubleadingJet_puJetIdMVA , wei )
    histStore['h2']["sublead"]['QGL'].Fill( eTree.h2SubleadingJet_QGL , wei )
    
    fillFlashggVars(histStore["flashggVars"],eTree,wei)
    LVStore=getLVStore(eTree)
    fillKinematicVarsFromLV(LVStore,histStore["kinematicVars"],wei)
 
def getBoostedLVs(eventCanditateLVs,boostVec):
    
    localLVs=copy.deepcopy(eventCanditateLVs)
    
    for ky in localLVs:
        localLVs[ky].Boost(boostVec)

    return localLVs
    
def fillKinematicVars( histCollection , LVStore , wei=1):
    
    histCollection['H1Pt'].Fill( LVStore['H1LV'].Pt() , wei ) 
    histCollection['H2Pt'].Fill( LVStore['H2LV'].Pt() , wei ) 
    histCollection['H3Pt'].Fill( LVStore['H3LV'].Pt() , wei ) 
    
    histCollection['H1Eta'].Fill( LVStore['H1LV'].Eta() , wei ) 
    histCollection['H2Eta'].Fill( LVStore['H2LV'].Eta() , wei ) 
    histCollection['H3Eta'].Fill( LVStore['H3LV'].Eta() , wei ) 
    
    histCollection['H1y'].Fill( LVStore['H1LV'].Rapidity() , wei ) 
    histCollection['H2y'].Fill( LVStore['H2LV'].Rapidity() , wei ) 
    histCollection['H3y'].Fill( LVStore['H3LV'].Rapidity() , wei ) 
    
    histCollection['H1Phi'].Fill( LVStore['H1LV'].Phi() , wei ) 
    histCollection['H2Phi'].Fill( LVStore['H2LV'].Phi() , wei ) 
    histCollection['H3Phi'].Fill( LVStore['H3LV'].Phi() , wei ) 
    
    histCollection['H1H2_DR'].Fill( LVStore['H1LV'].DeltaR(LVStore['H2LV']) ,  wei )
    histCollection['H1H3_DR'].Fill( LVStore['H1LV'].DeltaR(LVStore['H3LV']) ,  wei )
    histCollection['H2H3_DR'].Fill( LVStore['H2LV'].DeltaR(LVStore['H3LV']) ,  wei )
    
    histCollection['H1H2VsH3_CosTheta'].Fill(np.cos( (LVStore['H1LV'] + LVStore['H2LV']).Angle(LVStore['H3LV'].Vect())) , wei )
    histCollection['H1H3VsH2_CosTheta'].Fill(np.cos( (LVStore['H1LV'] + LVStore['H3LV']).Angle(LVStore['H2LV'].Vect())) , wei )
    histCollection['H2H3VsH1_CosTheta'].Fill(np.cos( (LVStore['H2LV'] + LVStore['H3LV']).Angle(LVStore['H1LV'].Vect())) , wei )
 
    histCollection['H1VsH2_CosTheta'].Fill(np.cos( LVStore['H1LV'].Angle(LVStore['H2LV'].Vect()) ), wei )
    histCollection['H1VsH3_CosTheta'].Fill(np.cos( LVStore['H1LV'].Angle(LVStore['H3LV'].Vect()) ), wei )
    histCollection['H2VsH3_CosTheta'].Fill(np.cos( LVStore['H2LV'].Angle(LVStore['H3LV'].Vect()) ), wei )
    
    histCollection['H1H2_dPhi'].Fill( dPhi( LVStore['H1LV'].Phi(),LVStore['H2LV'].Phi()) , wei )
    histCollection['H1H3_dPhi'].Fill( dPhi( LVStore['H1LV'].Phi(),LVStore['H3LV'].Phi()) , wei )
    histCollection['H2H3_dPhi'].Fill( dPhi( LVStore['H2LV'].Phi(),LVStore['H3LV'].Phi()) , wei )
     
    histCollection['H1H2_dEta'].Fill( dEta( LVStore['H1LV'].Eta(),LVStore['H2LV'].Eta()) , wei )
    histCollection['H1H3_dEta'].Fill( dEta( LVStore['H1LV'].Eta(),LVStore['H3LV'].Eta()) , wei )
    histCollection['H2H3_dEta'].Fill( dEta( LVStore['H2LV'].Eta(),LVStore['H3LV'].Eta()) , wei )

    histCollection['H1H2_area'].Fill( abs(0.5*LVStore['H1LV'].Vect().Cross(LVStore['H2LV'].Vect()).Mag()/125/125.0 ), wei )
    histCollection['H1H3_area'].Fill( abs(0.5*LVStore['H1LV'].Vect().Cross(LVStore['H3LV'].Vect()).Mag()/125/125.0 ), wei )
    histCollection['H2H3_area'].Fill( abs(0.5*LVStore['H2LV'].Vect().Cross(LVStore['H3LV'].Vect()).Mag()/125/125.0 ), wei )
    
    histCollection['H1H2H3_volume'].Fill( LVStore['H1LV'].Vect().Dot(LVStore['H2LV'].Vect().Cross(LVStore['H3LV'].Vect()))/125.0/125.0/125.0 , wei )
    histCollection['L1L2L3_volume'].Fill( LVStore['g1LV'].Vect().Dot(LVStore['j1LV'].Vect().Cross(LVStore['k1LV'].Vect()))/125.0/125.0/125.0 , wei )
    histCollection['sL1sL2sL3_volume'].Fill( LVStore['g2LV'].Vect().Dot(LVStore['j2LV'].Vect().Cross(LVStore['k2LV'].Vect()))/125.0/125.0/125.0 , wei )
    varDict=getOtherDerivedVariables(eTree=None,LVStore=LVStore,quad=None)
    fillOtherDerivedVariables(histCollection,varDict,wei)

def getHHHJethistos():

    histStore={'jets':{'nEvents':{}}}
    histStore['jets']['nEvents']=ROOT.TH1F("nEvents","",40,-0.5,39.5)
    histStore['jets']["nEvents"].SetCanExtend(ROOT.TH1.kAllAxes);

    
    for jet in['h1j1','h1j2','h2j1','h2j2']:
        histStore['jets'][jet] ={}
        histStore['jets'][jet]['NHF']       = ROOT.TH1F("NHF"  ,"",100,-2.5,2.5)
        histStore['jets'][jet]['NEMF']      = ROOT.TH1F("NEMF" ,"",100,-2.5,2.5)
        histStore['jets'][jet]['CHF']       = ROOT.TH1F("CHF"  ,"",100,-2.5,2.5)
        histStore['jets'][jet]['MUF']       = ROOT.TH1F("MUF"  ,"",100,-2.5,2.5)
        histStore['jets'][jet]['CEMF']      = ROOT.TH1F("CEMF" ,"",100,-2.5,2.5)
        histStore['jets'][jet]['NumConst']  = ROOT.TH1F("NumConst","",200,-0.5,199.5)
        histStore['jets'][jet]['NumNeutralParticles']   = ROOT.TH1F("NumNeutralParticles","",200,-0.5,199.5)
        histStore['jets'][jet]['CHM']       = ROOT.TH1F("CHM","",100,-2.5,2.5)
        histStore['jets'][jet]['csvScore']       = ROOT.TH1F("csvScore","",100,-2.5,2.5)
        histStore['jets'][jet]['deepCSVScore']       = ROOT.TH1F("deepCSVScore","",100,-2.5,2.5)
        histStore['jets'][jet]['deepJetScore']       = ROOT.TH1F("deepJetScore","",100,-2.5,2.5)
        histStore['jets'][jet]['flavour']        = ROOT.TH1F("hFlav","",60,-30.5,29.5)
        histStore['jets'][jet]['pFlavour']       = ROOT.TH1F("pFlav","",60,-30.5,29.5)
        histStore['jets'][jet]['particleNetAK4_B']       = ROOT.TH1F("particleNetAK4_B"   ,"",100,-2.5,2.5)
        histStore['jets'][jet]['particleNetAK4_CvsL']    = ROOT.TH1F("particleNetAK4_CvsL","",100,-2.5,2.5) 
        histStore['jets'][jet]['particleNetAK4_CvsB']    = ROOT.TH1F("particleNetAK4_CvsB","",100,-2.5,2.5) 
        histStore['jets'][jet]['particleNetAK4_QvsG']    = ROOT.TH1F("particleNetAK4_QvsG","",100,-2.5,2.5) 
        histStore['jets'][jet]['particleNetAK4_puIdDisc']= ROOT.TH1F("particleNetAK4_puIdDisc","",100,-2.5,2.5) 
        for tg in ['Loose' , 'Tight' ,'Tight2017']:
            histStore['jets'][jet]['jetID_'+tg]= ROOT.TH1F("jetID_"+tg,"",5,-2.5,2.5) 

    return histStore

def fillHHHJetHisto(histDict,eTree,idx):
    histDict['NHF'].Fill(eTree.jets_NHF[idx])
    histDict['NEMF'].Fill(eTree.jets_NEMF[idx])
    histDict['CHF'].Fill(eTree.jets_CHF[idx])
    histDict['MUF'].Fill(eTree.jets_MUF[idx])
    histDict['CEMF'].Fill(eTree.jets_CEMF[idx])
    histDict['NumConst'].Fill(eTree.jets_NumConst[idx])
    histDict['NumNeutralParticles'].Fill(eTree.jets_NumNeutralParticles[idx])
    histDict['CHM'].Fill(eTree.jets_CHM[idx])
    histDict['csvScore'].Fill(eTree.jets_csvScore[idx])
    histDict['deepCSVScore'].Fill(eTree.jets_deepCSVScore[idx])
    histDict['deepJetScore'].Fill(eTree.jets_deepJetScore[idx])
    histDict['flavour'].Fill(eTree.jets_flavour[idx])
    histDict['pFlavour'].Fill(eTree.jets_pFlavour[idx])
    histDict['particleNetAK4_B'].Fill(eTree.jets_particleNetAK4_B[idx])
    histDict['particleNetAK4_CvsL'].Fill(eTree.jets_particleNetAK4_CvsL[idx])
    histDict['particleNetAK4_CvsB'].Fill(eTree.jets_particleNetAK4_CvsB[idx])
    histDict['particleNetAK4_QvsG'].Fill(eTree.jets_particleNetAK4_QvsG[idx])
    histDict['particleNetAK4_puIdDisc'].Fill(eTree.jets_particleNetAK4_puIdDisc[idx])
    for tg in ['Loose' , 'Tight' ,'Tight2017']:
        histDict['jetID_'+tg].Fill(jetID(eTree,idx,tg) )


def jetID( eTree , idx, case='Tight2017'  ):

    eta = eTree.jets_eta[idx]
    
    NHF = eTree.jets_NHF[idx]
    NEMF = eTree.jets_NEMF[idx]
    NumConst = eTree.jets_NumConst[idx]
    CHF = eTree.jets_CHF[idx]
    CEMF = eTree.jets_CEMF[idx]
    CHM = eTree.jets_CHM[idx]
    NumNeutralParticles = eTree.jets_NumNeutralParticles[idx]



    jetID_barrel_loose  =  (  NHF<0.99 and NEMF<0.99 and NumConst>1) and ((abs(eta)<=2.4 and CHF>0 and CHM>0 and CEMF<0.99) or abs(eta)>2.4) and abs(eta)<=2.7;
    jetID_barrel_tight  =  (  NHF<0.90 and NEMF<0.90 and NumConst>1) and ((abs(eta)<=2.4 and CHF>0 and CHM>0 and CEMF<0.99) or abs(eta)>2.4) and abs(eta)<=2.7;
    jetID_transition    =  (  NEMF>0.01 and NHF<0.98 and NumNeutralParticles>2 and abs(eta)>2.7 and abs(eta)<3.0);
    jetID_forward       =  (  NEMF<0.90 and NumNeutralParticles >10 and abs(eta)>3.0 );

    jetID_2017_27 = (NHF < 0.9 and NEMF < 0.9 and NumConst > 1);
    jetID_2017_24 = jetID_2017_27 and (CHF > 0. and CHM > 0 );
    jetID_2017_30 = (NEMF > 0.02 and NEMF < 0.99 and NumNeutralParticles > 2);
    jetID_2017_forward = (NEMF < 0.9 and NHF > 0.02 and NumNeutralParticles > 10);


    if case=='Loose':
        if(abs(eta)<=2.7 ):
            return jetID_barrel_loose;
        if(abs(eta)<=3.0 ):
            return jetID_transition;
        if(abs(eta)> 3.0 ): return jetID_forward;
    elif case=='Tight':
        if(abs(eta)<=2.7 ): return jetID_barrel_tight;
        if(abs(eta)<=3.0 ): return jetID_transition;
        if(abs(eta)> 3.0 ): return jetID_forward;
    elif case=='Tight2017':
        if(abs(eta)<2.4 ): return jetID_2017_24;
        if(abs(eta)<2.7 ): return jetID_2017_27;
        if(abs(eta)<3.0 ): return jetID_2017_30;
        return jetID_2017_forward;
    else:
        print("error:: wrong level !!")
        return -1;


#float TrippleHTag::getSigmaM1OverMJets() const
#{
#    float dijetSigmaMOverM = 1./pow(dijet1().M(),2)*sqrt(
#                                         pow(h1LeadJet().userFloat("bRegNNResolution"),2)*pow(pow(h1LeadJet().p4().M(),2) + h1LeadJet().p4().Dot(h1SubleadJet().p4()) ,2)  +
#                                         pow(h1SubleadJet().userFloat("bRegNNResolution"),2)*pow( pow(h1SubleadJet().p4().M(),2) + h1SubleadJet().p4().Dot(h1LeadJet().p4()),2)
#                                                     );
#    return dijetSigmaMOverM;
#}

def getSigmaMOverM( obj1_p4,  obj1_res , obj2_p4, obj2_res ):
    M=(obj1_p4 +obj2_p4).M()
    sMoM = np.sqrt( pow(obj1_res,2)*pow(pow(obj1_p4.M(),2) + obj1_p4.Dot(obj2_p4) ,2) + pow( obj2_res, 2 )*pow( pow(obj2_p4.M(),2) + obj2_p4.Dot(obj1_p4) ,2) ) / M**2
    return sMoM


def fillTrippleHRecoVariables(eTree,histStore, LVStore, quad):
    
    histStore["diPhotons"]["pT"].Fill(  LVStore['HggLV'].Pt())
    histStore["diPhotons"]["eta"].Fill( LVStore['HggLV'].Eta())
    histStore["diPhotons"]["y"].Fill(   LVStore['HggLV'].Rapidity())
    histStore["diPhotons"]["phi"].Fill( LVStore['HggLV'].Phi())
    histStore["diPhotons"]["mass"].Fill(  eTree.CMS_hgg_mass)
    
    detId=0
    if(  abs(LVStore['g1LV'].Eta()) < 1.4 and abs(LVStore['g2LV'].Eta()) < 1.4 ): detId=0
    elif(abs(LVStore['g1LV'].Eta()) < 1.4 and abs(LVStore['g2LV'].Eta()) > 1.4 ) : detId=1
    elif(abs(LVStore['g1LV'].Eta()) > 1.4 and abs(LVStore['g2LV'].Eta()) < 1.4 ):  detId=2
    elif(abs(LVStore['g1LV'].Eta()) > 1.4 and abs(LVStore['g2LV'].Eta()) > 1.4 ) : detId=3
    else:    detId=4

    histStore["diPhotons"]["g1g2_Dr"].Fill( LVStore['g2LV'].DeltaR(LVStore['g1LV']) )
    histStore["diPhotons"]["detIdx"].Fill(  detId)
    histStore["diPhotons"]["gamma1_pT"].Fill(  eTree.leadingPhoton_pt)
    histStore["diPhotons"]["gamma1_eta"].Fill( eTree.leadingPhoton_eta)
    histStore["diPhotons"]["gamma1_phi"].Fill( eTree.leadingPhoton_phi)
    histStore["diPhotons"]["gamma2_pT"].Fill(  eTree.subleadingPhoton_pt)
    histStore["diPhotons"]["gamma2_eta"].Fill( eTree.subleadingPhoton_eta)
    histStore["diPhotons"]["gamma2_phi"].Fill( eTree.subleadingPhoton_phi)

    i1,i2,i3,i4=quad['fgg_idxs']
    
    histStore["h1Tobb"]["mass"].Fill( quad["m1"] )
    histStore["h1Tobb"]["mass_preReg"].Fill( quad["m1_preReg"] )
    histStore["h1Tobb"]["pT"].Fill( quad["p4_h1"].Pt() )
    histStore["h1Tobb"]["eta"].Fill( quad["p4_h1"].Eta() )
    histStore["h1Tobb"]["y"].Fill( quad["p4_h1"].Rapidity() )
    histStore["h1Tobb"]["phi"].Fill( quad["p4_h1"].Phi() )
    
    histStore["h2Tobb"]["mass"].Fill( quad["m2"] )
    histStore["h2Tobb"]["mass_preReg"].Fill( quad["m2_preReg"] )
    histStore["h2Tobb"]["pT"].Fill( quad["p4_h2"].Pt() )
    histStore["h2Tobb"]["eta"].Fill( quad["p4_h2"].Eta() )
    histStore["h2Tobb"]["y"].Fill( quad["p4_h2"].Rapidity() )
    histStore["h2Tobb"]["phi"].Fill( quad["p4_h2"].Phi() )

    for htag,hJets in zip(['1','2'] , [[i1,i2],[i3,i4]]):
        for jtag, idx in zip(['1','2'] , hJets): 
            sidx=str(idx)
            histStore["h"+htag+"Tobb"]["bJet"+jtag+"_pT_preReg"].Fill( getattr(eTree,'jet_'+sidx+'_pt') )
            histStore["h"+htag+"Tobb"]["bJet"+jtag+"_pT"].Fill(        getattr(eTree,'jet_'+sidx+'_pt')*getattr(eTree,'jet_'+sidx+'_bJetRegCorr') )
            histStore["h"+htag+"Tobb"]["bJet"+jtag+"_eta"].Fill(       getattr(eTree,'jet_'+sidx+'_eta'))
            histStore["h"+htag+"Tobb"]["bJet"+jtag+"_phi"].Fill(       getattr(eTree,'jet_'+sidx+'_phi'))
            histStore["h"+htag+"Tobb"]["bJet"+jtag+"_deepCSVScore"].Fill( getattr(eTree,'jet_'+sidx+'_deepCSVScore') );
            histStore["h"+htag+"Tobb"]["bJet"+jtag+"_deepJetScore"].Fill( getattr(eTree,'jet_'+sidx+'_deepJetScore') );
            histStore["h"+htag+"Tobb"]["bJet"+jtag+"_hFlavour"].Fill( getattr(eTree,'jet_'+sidx+'_flavour' ) );
            histStore["h"+htag+"Tobb"]["bJet"+jtag+"_pFlavour"].Fill( getattr(eTree,'jet_'+sidx+'_pFlavour') );
    
    histStore["h1Tobb"]["b1b2_Dr"].Fill(deltaR( getattr(eTree,'jet_'+str(i1)+'_eta') , getattr(eTree,'jet_'+str(i1)+'_phi') , 
                                                getattr(eTree,'jet_'+str(i2)+'_eta') , getattr(eTree,'jet_'+str(i2)+'_phi') ) 
                                                )
    histStore["h2Tobb"]["b1b2_Dr"].Fill(deltaR( getattr(eTree,'jet_'+str(i3)+'_eta') , getattr(eTree,'jet_'+str(i3)+'_phi') , 
                                                getattr(eTree,'jet_'+str(i4)+'_eta') , getattr(eTree,'jet_'+str(i4)+'_phi') ) 
                                                )
 

    idxs=[]
    
    #hhTo4b
    histStore["hhTo4b"]["mass"].Fill( quad["mass"] )
    histStore["hhTo4b"]["mass_preReg"].Fill( quad["mass_preReg"] )
    histStore["hhTo4b"]["pT"].Fill( quad["pT"] )
    histStore["hhTo4b"]["eta"].Fill( quad["eta"] )
    histStore["hhTo4b"]["y"].Fill( quad["y"] )
    histStore["hhTo4b"]["phi"].Fill( quad["phi"] )
    
    #hhhTo4b2gamma
    allDauP4 = (quad['p4_h1']+quad['p4_h2']+LVStore['HggLV'])
    allDauP4_preReg = (quad['p4_h1_preReg']+quad['p4_h2_preReg']+LVStore['HggLV'])
    

    histStore["hhhTo4b2gamma"]["mass"].Fill( allDauP4.M() )
    histStore["hhhTo4b2gamma"]["mass_preReg"].Fill( allDauP4_preReg.M()  )
    histStore["hhhTo4b2gamma"]["pT"].Fill( allDauP4.Pt() )
    histStore["hhhTo4b2gamma"]["eta"].Fill( allDauP4.Eta() )
    histStore["hhhTo4b2gamma"]["y"].Fill( allDauP4.Rapidity() )
    histStore["hhhTo4b2gamma"]["phi"].Fill( allDauP4.Phi() )
    histStore["hhhTo4b2gamma"]["h1MassVsh2Mass"].Fill(quad['m1'],quad['m2'] )
    histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_preReg"].Fill(quad['m1_preReg'],quad['m2_preReg'] )
    
    mx0 =  allDauP4.M() - ( quad['p4_h1'].M()+quad['p4_h2'].M()+LVStore['HggLV'].M() - 375.0 )
    mx1 = (quad['p4_h1']+quad['p4_h2']).M() - ( quad['p4_h1'].M()+quad['p4_h2'].M() - 250.0 ) 
    mx2 = (quad['p4_h1']+LVStore['HggLV']).M()-( quad['p4_h1'].M()+LVStore['HggLV'].M() - 250.0 ) 
    mx3 = (quad['p4_h2']+LVStore['HggLV']).M()-( quad['p4_h2'].M()+LVStore['HggLV'].M() - 250.0 ) 

    histStore["hhhTo4b2gamma"]["massX0"].Fill( mx0 )
    histStore["hhhTo4b2gamma"]["massX1"].Fill( mx1 )
    histStore["hhhTo4b2gamma"]["massX2"].Fill( mx2 )
    histStore["hhhTo4b2gamma"]["massX3"].Fill( mx3 )

    histStore["hhhTo4b2gamma"]["massX0X1"].Fill( mx0,mx1  )
    histStore["hhhTo4b2gamma"]["massX2X3"].Fill( mx2,mx3  )
    histStore["hhhTo4b2gamma"]["massX0X2"].Fill( mx0,mx2  )
    histStore["hhhTo4b2gamma"]["massX1X3"].Fill( mx1,mx3  )
    histStore["hhhTo4b2gamma"]["massX0X3"].Fill( mx0,mx3  )
    histStore["hhhTo4b2gamma"]["massX1X2"].Fill( mx1,mx2  )

def visualizeEvents(eTree ,text=None,outputFname=None ,jetMask=[]) :
    
    artists=[]
    legends=[]
    fig, ax = plt.subplots(figsize=(11,8))
    jetX=[]
    jetY=[]
    if jetMask==[]:
        jetMask=[True for i in range(N_JET_MAX) ]
    for i in range(8):
        if abs(getattr(eTree,'jet_'+str(i)+'_isValid') )  < 0.5:
            continue
        eta,phi=getattr(eTree,'jet_'+str(i)+'_eta') , getattr(eTree,'jet_'+str(i)+'_phi')
        jetX.append(eta)
        jetY.append(phi)
        cir=plt.Circle((eta,phi),0.4,color='b',fill='b',alpha=0.9-i*0.1)
        artists.append( cir )
        lbl='j'+str(i)+' : '+str(np.round(getattr(eTree,'jet_'+str(i)+'_deepJetScore') ,2))
        lbl+=' , '+str(np.round(getattr(eTree,'jet_'+str(i)+'_pt') ,1))+' GeV'
        lbl+=' ['+str(np.round(getattr(eTree,'jet_'+str(i)+'_flavour') ,0))+']'
        if jetMask[i]:
            lbl=u'\u2713 ' +lbl
        else:
            lbl=u'\u2717 ' +lbl

        legends.append( lbl )
        ax.text(eta+0.45,phi,"j"+str(i))
    eta,phi=eTree.leadingPhoton_eta,eTree.leadingPhoton_phi
    cir=plt.Circle((eta,phi),0.5,color='darkorange',fill='darkorange',alpha=0.4,label='leadG')
    artists.append( cir ); legends.append('leadG  '+str(np.round(eTree.leadingPhoton_pt,1))+' GeV')

    
    eta,phi=eTree.subleadingPhoton_eta,eTree.subleadingPhoton_phi
    cir=plt.Circle( (eta,phi),0.5,color='orange',fill='orange',alpha=0.2,label='subLeadG ' )
    artists.append( cir ); legends.append( 'subleadG '+str(np.round(eTree.subleadingPhoton_pt,1))+' GeV' )
    
    bDaus=hhhSelector.getHiggsDauP4s(eTree,5)
    gDaus=hhhSelector.getHiggsDauP4s(eTree,22)
    
    cir=plt.Circle((gDaus[0].Eta(),gDaus[0].Phi()),0.1,color='forestgreen',alpha=0.8); legends.append('geb $\gamma_{1}$ '+str(np.round(gDaus[0].Pt(),1))+' GeV');artists.append( cir ) ;
    cir=plt.Circle((gDaus[1].Eta(),gDaus[1].Phi()),0.1,color='limegreen',alpha=0.6); legends.append('gen $\gamma_{2}$ '+str(np.round(gDaus[1].Pt(),1))+' GeV');artists.append( cir ) ; 
    cir=plt.Circle((bDaus[0].Eta(),bDaus[0].Phi()),0.1,color='magenta' ,alpha=0.8);legends.append('gen b1 ' + str(np.round(bDaus[0].Pt(),1))+' GeV');    artists.append( cir )
    cir=plt.Circle((bDaus[1].Eta(),bDaus[1].Phi()),0.1,color='plum' ,alpha=0.8);legends.append('gen b2 ' + str(np.round(bDaus[1].Pt(),1))+' GeV');    artists.append( cir )
    cir=plt.Circle((bDaus[2].Eta(),bDaus[2].Phi()),0.1,color='orangered',alpha=0.8);legends.append('gen b3 ' + str(np.round(bDaus[2].Pt(),1))+' GeV');    artists.append( cir )
    cir=plt.Circle((bDaus[3].Eta(),bDaus[3].Phi()),0.1,color='coral'    ,alpha=0.8);legends.append('gen b4 ' + str(np.round(bDaus[3].Pt(),1))+' GeV');    artists.append( cir )
       
    ax.set_aspect( 1 )
    ax.set_xlim([-3.0,6.6])
    ax.set_ylim([-3.8,3.2])
    if text:   
        ax.text(-2.7,-3.6,text,zorder=100,backgroundcolor='w')
        
    for art in artists:
        ax.add_artist(art)
    plt.scatter(jetX,jetY,c='k',marker='x',label='reco jet centers')
    ax.legend(artists,legends,loc='upper right')
    ax.axhline( np.pi,color='k',linestyle='dashed')
    ax.axhline( -np.pi,color='k',linestyle='dashed')
    ax.axvline( 2.5,color='k',linestyle='dashed')
    ax.axvline(-2.5,color='k',linestyle='dashed')
    if outputFname:
        #fig.savefig(outputFname)
        fig.savefig(outputFname,bbox_inches='tight')
    plt.close(fig)

def printEventInfo(eTree,jetMask=[]):
    if jetMask==[]:
        jetMask=[True  for i in range(N_JET_MAX)]
    
    print("-"*20)
    print("Event ID : ",eTree.event)
    
    print("lead photon [pt/eta/phi ]   : ", np.round(eTree.leadingPhoton_pt ,1),
                                         np.round(eTree.leadingPhoton_eta,2),
                                         np.round(eTree.leadingPhoton_phi,2))

    print("sublead photon [pt/eta/phi] : ", np.round(eTree.subleadingPhoton_pt ,1),
                                         np.round(eTree.subleadingPhoton_eta,2),
                                         np.round(eTree.subleadingPhoton_phi,2))

    for i in range(N_JET_MAX ):
        if (getattr(eTree,'jet_'+str(i)+'_isValid') < 0.25):
            continue
        isMasked='V' 
        if not jetMask[i]:
            isMasked='X' 
            
        print("Jet ",i,"| ",isMasked," |  : ", np.round(getattr(eTree,'jet_'+str(i)+'_pt') ,1),
                              np.round(getattr(eTree,'jet_'+str(i)+'_eta'),2),
                              np.round(getattr(eTree,'jet_'+str(i)+'_phi'),2) )
 
def printEventInfoFromLVStore(LVStore,eid=None):
    print("-"*20)
    if eid:
        print("Event ID : ",eid)
    
    print("lead photon [pt/eta/phi ]   : ", np.round(LVStore['g1LV'].Pt(),1),
                                            np.round(LVStore['g1LV'].Eta(),2),
                                            np.round(LVStore['g1LV'].Phi(),2))

    print("sublead photon [pt/eta/phi ]   : ", np.round(LVStore['g2LV'].Pt(),1),
                                               np.round(LVStore['g2LV'].Eta(),2),
                                               np.round(LVStore['g2LV'].Phi(),2))

    print("b1 [pt/eta/phi] : ",np.round(LVStore['j1LV'].Pt(),1),
                               np.round(LVStore['j1LV'].Eta(),2),
                               np.round(LVStore['j1LV'].Phi(),2))
    print("b2 [pt/eta/phi] : ",np.round(LVStore['j2LV'].Pt(),1),
                               np.round(LVStore['j2LV'].Eta(),2),
                               np.round(LVStore['j2LV'].Phi(),2))
    print("b3 [pt/eta/phi] : ",np.round(LVStore['k1LV'].Pt(),1),
                               np.round(LVStore['k1LV'].Eta(),2),
                               np.round(LVStore['k1LV'].Phi(),2))
    print("b4 [pt/eta/phi] : ",np.round(LVStore['k2LV'].Pt(),1),
                               np.round(LVStore['k2LV'].Eta(),2),
                               np.round(LVStore['k2LV'].Phi(),2))
def getDataTag(fname):
    
    txt=fname.lower()
    dataType='data'
    if 'diphotonjetsbox1bjet' in txt:
        dataType='ggJets1b'
    elif 'diphotonjetsbox2bjet' in txt:
        dataType='ggJets2b'
    elif 'diphotonjetsbox' in txt:
        dataType='ggJets'

    return dataType

def getNBsFromQuad(eTree,quad):
    bCount=0
    if  getattr(eTree,'jet_'+str(quad['fgg_idxs'][0])+'_flavour')==5:
        bCount+=1
    if  getattr(eTree,'jet_'+str(quad['fgg_idxs'][1])+'_flavour')==5:
        bCount+=1
    if  getattr(eTree,'jet_'+str(quad['fgg_idxs'][2])+'_flavour')==5:
        bCount+=1
    if  getattr(eTree,'jet_'+str(quad['fgg_idxs'][3])+'_flavour')==5:
        bCount+=1

    return bCount



