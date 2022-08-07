from __future__ import print_function
import ROOT ,copy
import numpy as np
import itertools as itrTools
from Util import *



class mcScaler:
    def __init__(self):
        self.dataHist=None
        self.mcHist=None
        self.binEdges=[]
        self.scaleFactors=[]
        self.def_scaleFactor=0
        self.underflow_x=0.0
        self.overflow_x=1.0
        self.underflow_scaleFactor=1.0
        self.overflow_scaleFactor=1.0
        self.scaleFactorHist=0
        self.scaleFactorHistSet=False
        self.File=None
    def setSFHist(self,hist):
        self.scaleFactorHist=hist.Clone()
        self.underflow_x= self.scaleFactorHist.GetBinLowEdge(1)
        self.overflow_x= self.scaleFactorHist.GetBinLowEdge(self.scaleFactorHist.GetNbinsX())
        self.underflow_scaleFactor = self.scaleFactorHist.GetBinContent(1)
        self.overflow_scaleFactor = self.scaleFactorHist.GetBinContent(self.scaleFactorHist.GetNbinsX()-1)
        print("Scale factors defined for x = [ ",self.underflow_x,self.overflow_x," ]")
        print("Under flow Scale factor     =  ",self.underflow_scaleFactor)
        print("Over flow Scale factor      = ",self.overflow_scaleFactor)

        self.scaleFactorHistSet=True
    def setSFHistFromFile(self,fname,histName):
        print("Setting SF hist ",histName," from ",fname,)
        self.File=ROOT.TFile(fname,"READ")
        hist=self.File.Get(histName)
        self.setSFHist(hist)
    
    def getSFForX(self,x):
        if not self.scaleFactorHistSet:
            print("\tWARNING ! : scaleFactor Hist not set , returning -1.0 ")
            return -1.0
        if x < self.underflow_x:
            return self.underflow_scaleFactor
        if x > self.overflow_x:
            return self.overflow_x
        bid=self.scaleFactorHist.FindBin(x)
        #print("bid = ",bid," for x ",x, " Bin edges --> ",
        #         self.scaleFactorHist.GetBinLowEdge(bid), 
        #         self.scaleFactorHist.GetBinLowEdge(bid)+self.scaleFactorHist.GetBinWidth(bid))
        return self.scaleFactorHist.GetBinContent(bid)
        
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
    genPt=[1e4 for i in range(4)]
    genEta=[1e4 for i in range(4)]
    genPhi=[1e4 for i in range(4)]
    genE=[1e4 for i in range(4)]
    k=0
    if abs(eTree.gen_H1_dau1_pdgId)==pdgId:
        genPt[k]=eTree.gen_H1_dau1_pt
        genEta[k]=eTree.gen_H1_dau1_eta
        genPhi[k]=eTree.gen_H1_dau1_phi
        genE[k]=eTree.gen_H1_dau1_e
        k+=1
        genPt[k]=eTree.gen_H1_dau2_pt
        genEta[k]=eTree.gen_H1_dau2_eta
        genPhi[k]=eTree.gen_H1_dau2_phi
        genE[k]=eTree.gen_H1_dau2_e
        k+=1
    if abs(eTree.gen_H2_dau1_pdgId)==pdgId:
        genPt[k]=eTree.gen_H2_dau1_pt
        genEta[k]=eTree.gen_H2_dau1_eta
        genPhi[k]=eTree.gen_H2_dau1_phi
        genE[k]=eTree.gen_H2_dau1_e
        k+=1
        genPt[k]=eTree.gen_H2_dau2_pt
        genEta[k]=eTree.gen_H2_dau2_eta
        genPhi[k]=eTree.gen_H2_dau2_phi
        genE[k]=eTree.gen_H2_dau2_e
        k+=1

    if abs(eTree.gen_H3_dau1_pdgId)==pdgId:
        genPt[k]=eTree.gen_H3_dau1_pt
        genEta[k]=eTree.gen_H3_dau1_eta
        genPhi[k]=eTree.gen_H3_dau1_phi
        genE[k]=eTree.gen_H3_dau1_e
        k+=1
        genPt[k]=eTree.gen_H3_dau2_pt
        genEta[k]=eTree.gen_H3_dau2_eta
        genPhi[k]=eTree.gen_H3_dau2_phi
        genE[k]=eTree.gen_H3_dau2_e
        k+=1
    allP4=[]
    for i in range(k):
      lv=ROOT.TLorentzVector()
      lv.SetPtEtaPhiE(genPt[i],genEta[i],genPhi[i],genE[i])
      allP4.append(lv)
    return allP4

def addFlasggVars(histStore): 
    
    histStore["flashggVars"]={}

    histStore["flashggVars"]['HHHCosThetaH1']= ROOT.TH1F("HHHCosThetaH1","",60,-0.1,1.1)
    histStore["flashggVars"]['HHHCosThetaHgg']= ROOT.TH1F("HHHCosThetaHgg","",60,-0.1,1.1)
    histStore["flashggVars"]['H4bCosThetaLeadJet']= ROOT.TH1F("H4bCosThetaLeadJet","",60,-0.1,1.1)
    histStore["flashggVars"]['HggCosThetaLeadGamma']= ROOT.TH1F("HggCosThetaLeadGamma","",60,-0.1,1.1)
    histStore["flashggVars"]['H1bbCosThetaLeadJet'] = ROOT.TH1F("H1bbCosThetaLeadJet"  ,"",60,-0.1,1.1)
    histStore["flashggVars"]['H2bbCosThetaLeadJet'] = ROOT.TH1F("H2bbCosThetaLeadJet"  ,"",60,-0.1,1.1)
    
    histStore["flashggVars"]['pTgg_overMgg'] = ROOT.TH1F("pTgg_overMgg"  ,"",100,0.0,10.0)
    histStore["flashggVars"]['pTh1jj_overMh1'] = ROOT.TH1F("pTh1jj_overMh1"  ,"",100,0.0,10.0)
    histStore["flashggVars"]['pTh2jj_overMh2'] = ROOT.TH1F("pTh2jj_overMh2"  ,"",100,0.0,10.0)
    
    histStore["flashggVars"]['pTleadG_overMgg'] = ROOT.TH1F("pTleadG_overMgg"  ,"",100,0.0,10.0)
    histStore["flashggVars"]['pTh1leadJ_overMh1'] = ROOT.TH1F("pTh1leadJ_overMh1"  ,"",100,0.0,10.0)
    histStore["flashggVars"]['pTh2leadJ_overMh2'] = ROOT.TH1F("pTh2leadJ_overMh2"  ,"",100,0.0,10.0)

    histStore["flashggVars"]['pTsubleadG_overMgg'] = ROOT.TH1F("pTsubleadG_overMgg"  ,"",100,0.0,10.0)
    histStore["flashggVars"]['pTh1subleadJ_overMh1'] = ROOT.TH1F("pTh1subleadJ_overMh1"  ,"",100,0.0,10.0)
    histStore["flashggVars"]['pTh2subleadJ_overMh2'] = ROOT.TH1F("pTh2subleadJ_overMh2"  ,"",100,0.0,10.0)

    histStore["flashggVars"]['h1leadJ_deepjetScore'] = ROOT.TH1F("h1leadJ_deepjetScore"  ,"",12,-0.1,1.1)
    histStore["flashggVars"]['h1subleadJ_deepjetScore'] = ROOT.TH1F("h1subleadJ_deepjetScore"  ,"",12,-0.1,1.1)
    histStore["flashggVars"]['h2leadJ_deepjetScore'] = ROOT.TH1F("h2leadJ_deepjetScore"  ,"",12,-0.1,1.1)
    histStore["flashggVars"]['h2subleadJ_deepjetScore'] = ROOT.TH1F("h2subleadJ_deepjetScore"  ,"",12,-0.1,1.1)

    histStore["flashggVars"]['h1leadG_customMVA'] = ROOT.TH1F("h1leadG_customMVA"  ,"",48,-0.1,1.1)
    histStore["flashggVars"]['h1subleadG_customMVA'] = ROOT.TH1F("h1subleadG_customMVA"  ,"",48,-0.1,1.1)
    
    histStore["flashggVars"]['h1leadG_SigOverE'] = ROOT.TH1F("h1leadG_SigOverE"  ,"",50,0.0,0.1)
    histStore["flashggVars"]['h1subleadG_SigOverE'] = ROOT.TH1F("h1subleadG_SigOverE"  ,"",50,0.0,0.1)
    histStore["flashggVars"]['hgg_SigMOverM'] = ROOT.TH1F("hgg_SigMOverM"  ,"",20,0.0,0.1)
    histStore["flashggVars"]['hgg_SigMOverMDecorrelated'] = ROOT.TH1F("hgg_SigMOverMDecorrelated"  ,"",20,0.0,0.1)
    histStore["flashggVars"]['h1_dijetSigmaMOverM'] = ROOT.TH1F("h1_dijetSigmaMOverM"  ,"",50,0.0,0.5)
    histStore["flashggVars"]['h2_dijetSigmaMOverM'] = ROOT.TH1F("h2_dijetSigmaMOverM"  ,"",50,0.0,0.5)
    
    histStore["flashggVars"]['rho'] = ROOT.TH1F("rho","",50,0.0,50.0)
    histStore["flashggVars"]['phoJetMinDr'] = ROOT.TH1F("phoJetMinDr"  ,"",60,0.0,6.0)
    histStore["flashggVars"]['otherphoJetMinDr'] = ROOT.TH1F("otherphoJetMinDr"  ,"",60,0.0,6.0)

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
        histStore[ky]['hhh']['phi']=ROOT.TH1F("phi","",62,-3.14,3.14)

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
        addFlasggVars(histStore[ky])
        addKinematicVars(histStore[ky])


def addKinematicVars(histStore): 
    histStore['kinematicVars']={}
    framesToProbe=['CMS_RestFrame',
                   'HHH_RestFrame',
                   'H1H2_RestFrame',
                   'H1H3_RestFrame',
                   'H2H3_RestFrame',
                   'Hgg_RestFrame',
                   'H1bb_RestFrame',
                   'H2bb_RestFrame',
                    ]
    
    for frame in framesToProbe:
        histStore['kinematicVars'][frame]={}
        
        histStore['kinematicVars'][frame]['H1H2_DR']  = ROOT.TH1F("H1H2DeltaR","",50,0.0,5.0) 
        histStore['kinematicVars'][frame]['H1H3_DR']  = ROOT.TH1F("H1H3DeltaR","",50,0.0,5.0) 
        histStore['kinematicVars'][frame]['H2H3_DR']  = ROOT.TH1F("H2H3DeltaR","",50,0.0,5.0) 

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
    
def fillKinematicVarsFromLV(LVStore,histStore,wei=1 ):

        
    fillKinematicVars(histStore['CMS_RestFrame'],LVStore)    

    boost=-1.0*LVStore['HHHLV'].BoostVector()
    boostedLV=getBoostedLVs(LVStore,boost)
    fillKinematicVars(histStore['HHH_RestFrame'],boostedLV , wei )    
    
    boost=-1.0*(LVStore['H1LV'] +LVStore['H2LV']).BoostVector()
    boostedLV=getBoostedLVs(LVStore,boost)
    fillKinematicVars(histStore['H1H2_RestFrame'],boostedLV , wei )    
    
    boost=-1.0*(LVStore['H1LV']+LVStore['H2LV']).BoostVector()
    boostedLV=getBoostedLVs(LVStore,boost)
    fillKinematicVars(histStore['H1H3_RestFrame'],boostedLV , wei )    
    
    boost=-1.0*(LVStore['H2LV']+LVStore['H3LV']).BoostVector()
    boostedLV=getBoostedLVs(LVStore,boost)
    fillKinematicVars(histStore['H2H3_RestFrame'],boostedLV , wei )    
    
    boost=-1.0*LVStore['HggLV'].BoostVector()
    boostedLV=getBoostedLVs(LVStore,boost)
    fillKinematicVars(histStore['Hgg_RestFrame'],boostedLV , wei )    
    
    boost=-1.0*LVStore['H1bbLV'].BoostVector()
    boostedLV=getBoostedLVs(LVStore,boost)
    fillKinematicVars(histStore['H1bb_RestFrame'],boostedLV , wei )    

    
    boost=-1.0*LVStore['H2bbLV'].BoostVector()
    boostedLV=getBoostedLVs(LVStore,boost)
    fillKinematicVars(histStore['H2bb_RestFrame'],boostedLV , wei )    

    return



def getRecoHistos():
    histStore={}
    histStore['events']={}
    histStore['events']['nEvents']=ROOT.TH1F("nEvents","",40,-0.5,39.5)
    histStore['events']["nEvents"].SetCanExtend(ROOT.TH1.kAllAxes);
    
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
    
    histStore["h1Tobb"]["mass_tightID" ] =histStore["h1Tobb"]["mass" ].Clone()
    histStore["h2Tobb"]["mass_tightID" ] =histStore["h2Tobb"]["mass" ].Clone()
    histStore["h1Tobb"]["mass_tightID" ].SetName("mass_tightID" )
    histStore["h2Tobb"]["mass_tightID" ].SetName("mass_tightID" )
    
    histStore["h1Tobb"]["mass_mediumID" ] =histStore["h1Tobb"]["mass" ].Clone()
    histStore["h2Tobb"]["mass_mediumID" ] =histStore["h2Tobb"]["mass" ].Clone()
    histStore["h1Tobb"]["mass_mediumID" ].SetName("mass_mediumID" )
    histStore["h2Tobb"]["mass_mediumID" ].SetName("mass_mediumID")
    
    histStore["h1Tobb"]["mass_looseID" ] =histStore["h1Tobb"]["mass" ].Clone()
    histStore["h2Tobb"]["mass_looseID" ] =histStore["h2Tobb"]["mass" ].Clone()
    histStore["h1Tobb"]["mass_looseID" ].SetName("mass_looseID" )
    histStore["h2Tobb"]["mass_looseID" ].SetName("mass_looseID" )



    for frame in ['vars_v1']:
        histStore[frame]={}
   
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
    LVStore['k1LV'].SetPtEtaPhiM(eTree.h2LeadingJet_pt,eTree.h1LeadingJet_eta,eTree.h1LeadingJet_phi,eTree.h1LeadingJet_mass);
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

def fillFlashggVars(histStore,eTree,wei):
    
    histStore['HHHCosThetaH1'].Fill( eTree.absCosThetaStar_CS , wei )
    
    LVStore = getLVStore(eTree)
    a,b,c,d,e=getCosthetaVars(eTree,LVStore)
    
    histStore['HHHCosThetaHgg'].Fill( c , wei ) #TODO
    histStore['H4bCosThetaLeadJet'].Fill( eTree.absCosTheta_bb , wei )
    
    histStore['HggCosThetaLeadGamma'].Fill( eTree.absCosTheta_gg , wei )
    histStore['H1bbCosThetaLeadJet'].Fill( a , wei ) #TODO
    histStore['H2bbCosThetaLeadJet'].Fill( b , wei ) #TODO
    
    histStore['pTgg_overMgg'].Fill( eTree.diphotonCandidatePtOverdiHiggsM , wei )
    histStore['pTh1jj_overMh1'].Fill( eTree.dije1CandidatePtOverdiHiggsM , wei )
    histStore['pTh2jj_overMh2'].Fill( eTree.dije2CandidatePtOverdiHiggsM , wei )
    
    histStore['pTleadG_overMgg'].Fill( eTree.leadingPhoton_pt / eTree.CMS_hgg_mass , wei )
    histStore['pTh1leadJ_overMh1'].Fill( eTree.h1LeadingJet_pt / eTree.M1jj , wei )
    histStore['pTh2leadJ_overMh2'].Fill( eTree.h2LeadingJet_pt / eTree.M2jj , wei )

    histStore['pTsubleadG_overMgg'].Fill( eTree.subleadingPhoton_pt / eTree.CMS_hgg_mass , wei )
    histStore['pTh1subleadJ_overMh1'].Fill( eTree.h1SubleadingJet_pt / eTree.M1jj , wei )
    histStore['pTh2subleadJ_overMh2'].Fill( eTree.h2SubleadingJet_pt / eTree.M2jj , wei )


    histStore['h1leadJ_deepjetScore'].Fill( eTree.h1LeadingJet_DeepFlavour , wei )
    histStore['h2leadJ_deepjetScore'].Fill( eTree.h2LeadingJet_DeepFlavour , wei )
    histStore['h1subleadJ_deepjetScore'].Fill( eTree.h1SubleadingJet_DeepFlavour , wei )
    histStore['h2subleadJ_deepjetScore'].Fill( eTree.h2SubleadingJet_DeepFlavour , wei )

    histStore['h1leadG_customMVA'].Fill( eTree.customLeadingPhotonIDMVA , wei )
    histStore['h1subleadG_customMVA'].Fill( eTree.customSubLeadingPhotonIDMVA , wei )
    
    histStore['h1leadG_SigOverE'].Fill( eTree.leadingPhotonSigOverE , wei )
    histStore['h1subleadG_SigOverE'].Fill( eTree.subleadingPhotonSigOverE , wei )
    histStore['hgg_SigMOverM'].Fill( eTree.sigmaMOverM , wei )
    histStore['hgg_SigMOverMDecorrelated'].Fill( eTree.sigmaMOverMDecorr , wei )
    #a,b = getDijetResolutions(eTree , LVStore)
    histStore['h1_dijetSigmaMOverM'].Fill( eTree.sigmaM1Jets , wei ) #TODO
    histStore['h2_dijetSigmaMOverM'].Fill( eTree.sigmaM2Jets , wei ) #TODO
    
    histStore['rho'].Fill( eTree.rho , wei )
    histStore['phoJetMinDr'].Fill( eTree.PhoJetOtherDr , wei )
    histStore['otherphoJetMinDr'].Fill( e , wei )



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
    histStore['hhh']['mass'].Fill( eTree.triHiggs_mass , wei )
    
    histStore['hgg']['pT'].Fill( eTree.diphoton_pt , wei )
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

def getGenHistos():
    histStore={'vars':{}}
    
    histStore['vars']['H1Pt'] = ROOT.TH1F("H1Pt","",500,0.0,1000)
    histStore['vars']['H2Pt'] = ROOT.TH1F("H2Pt","",500,0.0,1000)
    histStore['vars']['H3Pt'] = ROOT.TH1F("H3Pt","",500,0.0,1000)

    histStore['vars']['H1Eta'] = ROOT.TH1F("H1Eta","",100,-5.0,5.0)
    histStore['vars']['H2Eta'] = ROOT.TH1F("H2Eta","",100,-5.0,5.0)
    histStore['vars']['H3Eta'] = ROOT.TH1F("H3Eta","",100,-5.0,5.0)

    histStore['vars']['H1y'] = ROOT.TH1F("H1y","",100,-5.0,5.0)
    histStore['vars']['H2y'] = ROOT.TH1F("H2y","",100,-5.0,5.0)
    histStore['vars']['H3y'] = ROOT.TH1F("H3y","",100,-5.0,5.0)

    histStore['vars']['H1Phi'] = ROOT.TH1F("H1Phi","",100,-3.14,3.14)
    histStore['vars']['H2Phi'] = ROOT.TH1F("H2Phi","",100,-3.14,3.14)
    histStore['vars']['H3Phi'] = ROOT.TH1F("H3Phi","",100,-3.14,3.14)

    histStore['vars']['B1Pt'] = ROOT.TH1F("B1Pt","",500,0.0,1000)
    histStore['vars']['B2Pt'] = ROOT.TH1F("B2Pt","",500,0.0,1000)
    histStore['vars']['B3Pt'] = ROOT.TH1F("B3Pt","",500,0.0,1000)
    histStore['vars']['B4Pt'] = ROOT.TH1F("B4Pt","",500,0.0,1000)

    histStore['vars']['B1Eta'] = ROOT.TH1F("B1Eta","",100,-5.0,5.0)
    histStore['vars']['B2Eta'] = ROOT.TH1F("B2Eta","",100,-5.0,5.0)
    histStore['vars']['B3Eta'] = ROOT.TH1F("B3Eta","",100,-5.0,5.0)
    histStore['vars']['B4Eta'] = ROOT.TH1F("B4Eta","",100,-5.0,5.0)

    histStore['vars']['B1Y'] = ROOT.TH1F("B1Y","",100,-5.0,5.0)
    histStore['vars']['B2Y'] = ROOT.TH1F("B2Y","",100,-5.0,5.0)
    histStore['vars']['B3Y'] = ROOT.TH1F("B3Y","",100,-5.0,5.0)
    histStore['vars']['B4Y'] = ROOT.TH1F("B4Y","",100,-5.0,5.0)

    histStore['vars']['B1Phi'] = ROOT.TH1F("B1Phi","",100,-3.14,3.14)
    histStore['vars']['B2Phi'] = ROOT.TH1F("B2Phi","",100,-3.14,3.14)
    histStore['vars']['B3Phi'] = ROOT.TH1F("B3Phi","",100,-3.14,3.14)
    histStore['vars']['B4Phi'] = ROOT.TH1F("B4Phi","",100,-3.14,3.14)



    histStore['vars']['nHiggsAbove300'] = ROOT.TH1F("nHiggsAbove300","",3,0.0,3.0)
    histStore['vars']['nHiggsAbove300'].SetCanExtend(ROOT.TH1.kAllAxes);
    
    histStore['vars']['gamma1Pt'] = ROOT.TH1F("gamma1Pt","",500,0.0,1000)
    histStore['vars']['gamma2Pt'] = ROOT.TH1F("gamma2Pt","",500,0.0,1000)

    histStore['vars']['gamma1Eta'] = ROOT.TH1F("gamma1Eta","",100,-5.0,5.0)
    histStore['vars']['gamma2Eta'] = ROOT.TH1F("gamma2Eta","",100,-5.0,5.0)

    histStore['vars']['gamma1y'] = ROOT.TH1F("gamma1y","",100,-5.0,5.0)
    histStore['vars']['gamma2y'] = ROOT.TH1F("gamma2y","",100,-5.0,5.0)

    histStore['vars']['gamma1Phi'] = ROOT.TH1F("gamma1Phi","",100,-3.14,3.14)
    histStore['vars']['gamma2Phi'] = ROOT.TH1F("gamma2Phi","",100,-3.14,3.14)



    histStore['vars']['H1H2DeltaR'] = ROOT.TH1F("H1H2DeltaR","",50,0.0,5.0)
    histStore['vars']['H2H3DeltaR'] = ROOT.TH1F("H2H3DeltaR","",50,0.0,5.0)
    histStore['vars']['H3H1DeltaR'] = ROOT.TH1F("H3H1DeltaR","",50,0.0,5.0)


    histStore['vars']['H1H2DeltaEta'] = ROOT.TH1F("H1H2DeltaEta","",50,0.0,5.0)
    histStore['vars']['H2H3DeltaEta'] = ROOT.TH1F("H2H3DeltaEta","",50,0.0,5.0)
    histStore['vars']['H3H1DeltaEta'] = ROOT.TH1F("H3H1DeltaEta","",50,0.0,5.0)


    histStore['vars']['H1H2DeltaPhi'] = ROOT.TH1F("H1H2DeltaPhi","",50,0.0,5.0)
    histStore['vars']['H2H3DeltaPhi'] = ROOT.TH1F("H2H3DeltaPhi","",50,0.0,5.0)
    histStore['vars']['H3H1DeltaPhi'] = ROOT.TH1F("H3H1DeltaPhi","",50,0.0,5.0)
    
    
    histStore['vars']['mass_X0'] = ROOT.TH1F("mass_X0","",500,0.0,2000.0)
    histStore['vars']['mass_X1'] = ROOT.TH1F("mass_X1","",500,0.0,2000.0)
    histStore['vars']['mass_X2'] = ROOT.TH1F("mass_X2","",500,0.0,2000.0)
    histStore['vars']['mass_X3'] = ROOT.TH1F("mass_X3","",500,0.0,2000.0)

    histStore['vars']['mass_4b'] = ROOT.TH1F("mass_4b","",500,0.0,2000.0)
    histStore['vars']['mass_2gamma'] = ROOT.TH1F("mass_2gamma","",500,0.0,2000.0)
    histStore['vars']['mass_2bb'] = ROOT.TH1F("mass_2bb","",500,0.0,2000.0)
    histStore['vars']['mass_4b2gamma'] = ROOT.TH1F("mass_4b2gamma","",500,0.0,2000.0)
    
    histStore['kinematicVars']={}

    for frame in ['vars_v1']:
        histStore[frame]={}
        histStore[frame]['HHHCosThetaHgg']= ROOT.TH1F("HHHCosThetaHgg","",44,-1.1,1.1)
        histStore[frame]['HggCosThetaLeadGamma']= ROOT.TH1F("HggCosThetaLeadGamma","",44,-1.1,1.1)
        histStore[frame]['H1bbCosThetaLeadJet'] = ROOT.TH1F("H1bbCosThetaLeadJet"  ,"",44,-1.1,1.1)
        histStore[frame]['H2bbCosThetaLeadJet'] = ROOT.TH1F("H2bbCosThetaLeadJet"  ,"",44,-1.1,1.1)

    addKinematicVars(histStore)

#for frame in ['HHH_RestFrame']:
    #    histStore[frame]={}
    #    histStore[frame]['cosThetaLeadGamma']= ROOT.TH1F("cosThetaLeadGamma","",44,-1.1,1.1)


    return histStore

 
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


def getHHHJethistos():

    histStore={'jets':{'nEvents':{}}}
    histStore['jets']['nEvents']=ROOT.TH1F("nEvents","",40,-0.5,39.5)
    histStore['jets']["nEvents"].SetCanExtend(ROOT.TH1.kAllAxes);

    
    for jet in['h1j1','h1j2','h2j1','h2j2']:
        histStore['jets'][jet] ={}
        histStore['jets'][jet]['NHF']       = ROOT.TH1F("NHF","",-2.5,2.5,100)
        histStore['jets'][jet]['NEMF']      = ROOT.TH1F("NEMF","",-2.5,2.5,100)
        histStore['jets'][jet]['CHF']       = ROOT.TH1F("CHF","",-2.5,2.5,100)
        histStore['jets'][jet]['MUF']       = ROOT.TH1F("MUF","",-2.5,2.5,100)
        histStore['jets'][jet]['CEMF']      = ROOT.TH1F("CEMF","",-2.5,2.5,100)
        histStore['jets'][jet]['NumConst']  = ROOT.TH1F("NumConst","",-0.5,199.5,200)
        histStore['jets'][jet]['NumNeutralParticles']   = ROOT.TH1F("NumNeutralParticles","",-0.5,199.5,200)
        histStore['jets'][jet]['CHM']       = ROOT.TH1F("CHM","",-2.5,2.5,100)
        histStore['jets'][jet]['csvScore']       = ROOT.TH1F("csvScore","",-2.5,2.5,100)
        histStore['jets'][jet]['deepCSVScore']       = ROOT.TH1F("deepCSVScore","",-2.5,2.5,100)
        histStore['jets'][jet]['deepJetScore']       = ROOT.TH1F("deepJetScore","",-2.5,2.5,100)
        histStore['jets'][jet]['flavour']        = ROOT.TH1F("hFlav","",-30.5,29.5,60)
        histStore['jets'][jet]['pFlavour']       = ROOT.TH1F("pFlav","",-30.5,29.5,60)
        histStore['jets'][jet]['particleNetAK4_B']       = ROOT.TH1F("particleNetAK4_B","",-2.5,2.5,100)
        histStore['jets'][jet]['particleNetAK4_CvsL']    = ROOT.TH1F("particleNetAK4_CvsL","",-2.5,2.5,100) 
        histStore['jets'][jet]['particleNetAK4_CvsB']    = ROOT.TH1F("particleNetAK4_CvsB","",-2.5,2.5,100) 
        histStore['jets'][jet]['particleNetAK4_QvsG']    = ROOT.TH1F("particleNetAK4_QvsG","",-2.5,2.5,100) 
        histStore['jets'][jet]['particleNetAK4_puIdDisc']= ROOT.TH1F("particleNetAK4_puIdDisc","",-2.5,2.5,100) 

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


