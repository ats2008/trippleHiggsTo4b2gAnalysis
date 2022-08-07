from __future__ import print_function
import ROOT 
import numpy as np
import itertools as itrTools
from Util import *
from trippleHiggsUtils import *



def fillTrippleHGenVariables(eTree,histStore,bLVs,gammaLVs):
    
    LVStore={}
    LVStore['j1LV']=bLVs[0]
    LVStore['j2LV']=bLVs[0]
    LVStore['k1LV']=bLVs[2]
    LVStore['k2LV']=bLVs[1]

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
    
    hgg_=ROOT.TLorentzVector(LVStore['HggLV']) ;  hgg_.Boost(-1.0*LVStore['HHHLV'].BoostVector())
    histStore['vars_v1']['HHHCosThetaHgg'].Fill( hgg_.CosTheta())
    
    lv=ROOT.TLorentzVector(LVStore['g1LV']) ;    lv.Boost(-1*LVStore['HggLV'].BoostVector())
    histStore['vars_v1']['HggCosThetaLeadGamma'].Fill(lv.CosTheta())
    lv=ROOT.TLorentzVector(LVStore['j1LV']) ;    lv.Boost(-1*LVStore['H1bbLV'].BoostVector())
    histStore['vars_v1']['H1bbCosThetaLeadJet'].Fill(lv.CosTheta())
    lv=ROOT.TLorentzVector(LVStore['k1LV']) ;    lv.Boost(-1*LVStore['H2bbLV'].BoostVector())
    histStore['vars_v1']['H2bbCosThetaLeadJet'].Fill(lv.CosTheta())
    
    
    fillKinematicVars(histStore['kimeaticVars']['CMS_RestFrame'],LVStore)    

    boost=-1.0*LVStore['HHHLV'].BoostVector()
    boostedLV=getBoostedLVs(LVStore,boost)
    fillKinematicVars(histStore['kimeaticVars']['HHH_RestFrame'],boostedLV)    
    
    boost=-1.0*(LVStore['H1LV'] +LVStore['H2LV']).BoostVector()
    boostedLV=getBoostedLVs(LVStore,boost)
    fillKinematicVars(histStore['kimeaticVars']['H1H2_RestFrame'],boostedLV)    
    
    boost=-1.0*(LVStore['H1LV']+LVStore['H2LV']).BoostVector()
    boostedLV=getBoostedLVs(LVStore,boost)
    fillKinematicVars(histStore['kimeaticVars']['H1H3_RestFrame'],boostedLV)    
    
    boost=-1.0*(LVStore['H2LV']+LVStore['H3LV']).BoostVector()
    boostedLV=getBoostedLVs(LVStore,boost)
    fillKinematicVars(histStore['kimeaticVars']['H2H3_RestFrame'],boostedLV)    
    
    boost=-1.0*LVStore['HggLV'].BoostVector()
    boostedLV=getBoostedLVs(LVStore,boost)
    fillKinematicVars(histStore['kimeaticVars']['Hgg_RestFrame'],boostedLV)    
    
    boost=-1.0*LVStore['H1bbLV'].BoostVector()
    boostedLV=getBoostedLVs(LVStore,boost)
    fillKinematicVars(histStore['kimeaticVars']['H1bb_RestFrame'],boostedLV)    

    
    boost=-1.0*LVStore['H2bbLV'].BoostVector()
    boostedLV=getBoostedLVs(LVStore,boost)
    fillKinematicVars(histStore['kimeaticVars']['H2bb_RestFrame'],boostedLV)    

    return
