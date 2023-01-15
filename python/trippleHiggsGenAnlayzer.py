from __future__ import print_function
import ROOT 
import numpy as np
import itertools as itrTools
from Util import *
import trippleHiggsUtils as hhhUtil
import trippleHiggsSelector as hhhSelector


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

    hhhUtil.addKinematicVars(histStore)

#for frame in ['HHH_RestFrame']:
    #    histStore[frame]={}
    #    histStore[frame]['cosThetaLeadGamma']= ROOT.TH1F("cosThetaLeadGamma","",44,-1.1,1.1)


    return histStore
def printGenInfo(eTree):
    bLVs=hhhSelector.getHiggsDauP4s(eTree,5)       
    gammaLVs=hhhSelector.getHiggsDauP4s(eTree,22)       
    print('  g1 pt/eta/phi : ',np.round(gammaLVs[0].Pt(),2 ), np.round(gammaLVs[0].Eta(),2 ), np.round(gammaLVs[0].Phi(),2 ))
    print('  g2 pt/eta/phi : ',np.round(gammaLVs[1].Pt(),2 ), np.round(gammaLVs[1].Eta(),2 ), np.round(gammaLVs[1].Phi(),2 ))
    print('  b1 pt/eta/phi : ',np.round(bLVs[0].Pt(),2 ), np.round(bLVs[0].Eta(),2 ), np.round(bLVs[0].Phi(),2 ))
    print('  b2 pt/eta/phi : ',np.round(bLVs[1].Pt(),2 ), np.round(bLVs[1].Eta(),2 ), np.round(bLVs[1].Phi(),2 ))
    print('  b3 pt/eta/phi : ',np.round(bLVs[2].Pt(),2 ), np.round(bLVs[2].Eta(),2 ), np.round(bLVs[2].Phi(),2 ))
    print('  b4 pt/eta/phi : ',np.round(bLVs[3].Pt(),2 ), np.round(bLVs[3].Eta(),2 ), np.round(bLVs[3].Phi(),2 ))

def getLVStoreFromGen(eTree):
    bLVs=hhhUtil.getHiggsDauP4s(eTree,5)       
    gammaLVs=hhhUtil.getHiggsDauP4s(eTree,22)       
    LVStore={}
    LVStore['j1LV']=bLVs[0]
    LVStore['j2LV']=bLVs[0]
    LVStore['k1LV']=bLVs[2]
    LVStore['k2LV']=bLVs[2]

    if LVStore['j2LV'].Pt() < bLVs[1].Pt():
        LVStore['j1LV']=bLVs[1]
    else:
        LVStore['j2LV']=bLVs[1]

    if LVStore['k2LV'].Pt() < bLVs[3].Pt():
        LVStore['k1LV']=bLVs[3]
    else:
        LVStore['k2LV']=bLVs[3]
        
    
    LVStore['g1LV']=gammaLVs[0]
    LVStore['g2LV']=gammaLVs[1]
    if LVStore['g2LV'].Pt() < gammaLVs[1].Pt():
        LVStore['g1LV']=gammaLVs[1]
    else:
        LVStore['g2LV']=gammaLVs[1]

    LVStore['HggLV'] = LVStore['g1LV'] + LVStore['g2LV']
    LVStore['H1bbLV']= LVStore['j1LV'] + LVStore['j2LV']
    LVStore['H2bbLV']= LVStore['k1LV'] + LVStore['k2LV']
    LVStore['HH4bLV']= LVStore['H2bbLV'] + LVStore['H1bbLV']
    
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
    
    LVStore['HHHLV']=LVStore['H1LV']+LVStore['H2LV']+LVStore['H3LV']
    
    return LVStore


def fillTrippleHGenVariablesFromLVStore(histStore,LVStore):
    

    hgg_=ROOT.TLorentzVector(LVStore['HggLV']) ;  hgg_.Boost(-1.0*LVStore['HHHLV'].BoostVector())
    histStore['vars_v1']['HHHCosThetaHgg'].Fill( hgg_.CosTheta())
    lv=ROOT.TLorentzVector(LVStore['g1LV']) ;    lv.Boost(-1*LVStore['HggLV'].BoostVector())
    histStore['vars_v1']['HggCosThetaLeadGamma'].Fill(lv.CosTheta())
    lv=ROOT.TLorentzVector(LVStore['j1LV']) ;    lv.Boost(-1*LVStore['H1bbLV'].BoostVector())
    histStore['vars_v1']['H1bbCosThetaLeadJet'].Fill(lv.CosTheta())
    lv=ROOT.TLorentzVector(LVStore['k1LV']) ;    lv.Boost(-1*LVStore['H2bbLV'].BoostVector())
    histStore['vars_v1']['H2bbCosThetaLeadJet'].Fill(lv.CosTheta())
    
    
    hhhUtil.fillKinematicVarsFromLV(LVStore,histStore['kinematicVars'])

 #   fillKinematicVars(histStore['kimeaticVars'],LVStore)    

 #   boost=-1.0*LVStore['HHHLV'].BoostVector()
 #   boostedLV=getBoostedLVs(LVStore,boost)
 #   fillKinematicVars(histStore['kimeaticVars']['HHH_RestFrame'],boostedLV)    
 #   
 #   boost=-1.0*(LVStore['H1LV'] +LVStore['H2LV']).BoostVector()
 #   boostedLV=getBoostedLVs(LVStore,boost)
 #   fillKinematicVars(histStore['kimeaticVars']['H1H2_RestFrame'],boostedLV)    
 #   
 #   boost=-1.0*(LVStore['H1LV']+LVStore['H2LV']).BoostVector()
 #   boostedLV=getBoostedLVs(LVStore,boost)
 #   fillKinematicVars(histStore['kimeaticVars']['H1H3_RestFrame'],boostedLV)    
 #   
 #   boost=-1.0*(LVStore['H2LV']+LVStore['H3LV']).BoostVector()
 #   boostedLV=getBoostedLVs(LVStore,boost)
 #   fillKinematicVars(histStore['kimeaticVars']['H2H3_RestFrame'],boostedLV)    
 #   
 #   boost=-1.0*LVStore['HggLV'].BoostVector()
 #   boostedLV=getBoostedLVs(LVStore,boost)
 #   fillKinematicVars(histStore['kimeaticVars']['Hgg_RestFrame'],boostedLV)    
 #   
 #   boost=-1.0*LVStore['H1bbLV'].BoostVector()
 #   boostedLV=getBoostedLVs(LVStore,boost)
 #   fillKinematicVars(histStore['kimeaticVars']['H1bb_RestFrame'],boostedLV)    

 #   
 #   boost=-1.0*LVStore['H2bbLV'].BoostVector()
 #   boostedLV=getBoostedLVs(LVStore,boost)
 #   fillKinematicVars(histStore['kimeaticVars']['H2bb_RestFrame'],boostedLV)    

    return
