import matplotlib.pyplot as plt
import mplhep as hep
import uproot3 as urt
import ROOT
import json
import numpy as np
import os,argparse
from Util import lumiMap
import trippleHiggsPlotter as plotterUtil

if __name__=='__main__':
    

    parser = argparse.ArgumentParser()
    parser.add_argument("--unblind", help="Unblind the Mgg spectrum", action='store_true' )
    args=parser.parse_args()

    unblind= args.unblind

    prefixBase='/home/aravind/cernbox/work/trippleHiggs/hhhTo4b2gamma/genAnalysis/python/analysis/results/plots/analysis/jan23/variables_v0Scaled/geometric/'
    fileDict={}
    with open('workarea/data/bdtNtuples/filelistToUse.json') as f:
        fileDict=json.load(f)
    yearsToProcess=['2018','2017','2016PreVFP','2016PostVFP','run2','2016']
    bkgToProcess=[ 'ggBox1Bjet','ggBox2Bjet', 'ggBox','gJet20To40','gJet40ToInf']

    rdataFrames={}
    for yr  in yearsToProcess:
        rdataFrames[yr]={'sig':{},'bkg':{},'data':{}}
    
        fileName=fileDict[yr]['sig']['ggHHH']
        treeName = "trees/ggHHH_125_13TeV"
        rdataFrames[yr]['sig']['ggHHH']= ROOT.RDataFrame(treeName, fileName)
        print("Registering datset : ggHHH , ",yr," withh tree",treeName)
    
        ky=list(fileDict[yr]['data'].keys())[0]
        fileName=fileDict[yr]['data'][ky]
        treeName = "trees/Data_13TeV_TrippleHTag_0"
        if unblind:
            rdataFrames[yr]['data']['data']= ROOT.RDataFrame(treeName, fileName)
            prefixBase+='/unblinded/'
        else:
            rdataFrames[yr]['data']['data']= ROOT.RDataFrame(treeName, fileName).Filter(' CMS_hgg_mass < 115.0  || CMS_hgg_mass >135.0')
        print("Registering datset : data , ",yr," withh tree",treeName)
    
        rdataFrames[yr]['bkg']={}
        treeName = "trees/bkg_13TeV_TrippleHTag_0"
        for bkg in bkgToProcess:
            if bkg not in fileDict[yr]['bkg']:
                print()
                print("FILE NOT FOUNG FOR : ",yr," background ",bkg)
                print()
                continue
            fileName = fileDict[yr]['bkg'][bkg]
            rdataFrames[yr]['bkg'][bkg]=ROOT.RDataFrame(treeName, fileName).Filter(' CMS_hgg_mass < 115.0  || CMS_hgg_mass >135.0')
            print("Registering datset : ",bkg," , ",yr," withh tree",treeName)

    varToBinMap={  
                   "H1H2JetAbsCosThetaMax"    : ("","",16,0.0,1.0),
                   "H1H2JetAbsCosThetaMin"    : ("","",16,0.0,1.0),
                   "H1bbCosTheta_hhhF"        : ("","",16,0.0,1.0),
                   "H1bbToH2bbAbsCosTheta"    : ("","",16,0.0,1.0),
                   "H2bbCosThetaLeadJet"      : ("","",16,0.0,1.0),                   
                   "H2bbCosTheta_hhhF"        : ("","",16,0.0,1.0),
                   "HH4bCosThetaLeadJet"      : ("","",16,0.0,1.0),
                   "HH4bCosThetaLeadJet_hhhF" : ("","",16,0.0,1.0),
                   "HH4bCosTheta_hhhF"        : ("","",16,0.0,1.0),
                   "HggCosTheta_hhhF"         : ("","",16,0.0,1.0),
                   "HggTo4bAbsCosTheta"       : ("","",16,0.0,1.0),
                   "LeadJetAbsCosThetaMax"    : ("","",16,0.0,1.0),
                   "LeadJetAbsCosThetaMin"    : ("","",16,0.0,1.0),
                   "CosThetaH1_hhhF"          : ("","",16,0.0,1.0),
                   "h1bbCosThetaLeadJet"      : ("","",16,0.0,1.0),
                   "absCosThetaH4bHgg"        : ("","",16,0.0,1.0),
                   "absCosThetaH4bHgg_hhhF"   : ("","",16,0.0,1.0),
                   "D_HH"                         : ("","",20,0.0,200.0), 
                   "r_HH"                         : ("","",20,0.0,300.0),
                   "H1H2JetDrMax"                 : ("","",45,0.0,6.0),
                   "H1H2JetDrMin"                 : ("","",45,0.0,6.0),
                   "LeadJetDrMaxWithOtherJets"    : ("","",45,0.0,6.0),
                   "LeadJetDrMinWithOtherJets"    : ("","",45,0.0,6.0),
                   "PhoJetMaxDr"                  : ("","",45,0.0,6.0),
                   "PhoJetMaxDrOther"             : ("","",45,0.0,6.0),
                   "PhoJetMinDr"                  : ("","",45,0.0,6.0),
                   "PhoJetMinDrOther"     : ("","",45,0.0,6.0)
                }

    varToProcess=[
                       "H1H2JetAbsCosThetaMax"   , 
                       "H1H2JetAbsCosThetaMin"   , 
                       "H1bbCosTheta_hhhF"       , 
                       "H1bbToH2bbAbsCosTheta"   , 
                       "H2bbCosThetaLeadJet"     , 
                       "H2bbCosTheta_hhhF"       , 
                       "HH4bCosThetaLeadJet"     , 
                       "HH4bCosThetaLeadJet_hhhF", 
                       "HH4bCosTheta_hhhF"       , 
                       "HggCosTheta_hhhF"        , 
                       "HggTo4bAbsCosTheta"      , 
                       "LeadJetAbsCosThetaMax"   , 
                       "LeadJetAbsCosThetaMin"   , 
                       "CosThetaH1_hhhF"         , 
                       "h1bbCosThetaLeadJet"     , 
                       "absCosThetaH4bHgg"       , 
                       "absCosThetaH4bHgg_hhhF"  , 
                ]
    allHistoDict={}
    for yr  in yearsToProcess:
        allHistoDict[yr]={}
        print("Procssing for year  : ",yr)
        for var in varToProcess:
            saveBase=prefixBase+'/'+yr+'/'
            os.system('mkdir -p '+saveBase)
            histStore={'sig':{},'bkg':{},'data':{}}
            print(" Var : ",var)
            print("   Year : ",yr,' Signal')
            for ky in rdataFrames[yr]['sig']:
                histStore['sig'][ky]=rdataFrames[yr]['sig'][ky].Histo1D(varToBinMap[var],var,'weight')
            print("   Year : ",yr,' Background')
            for ky in rdataFrames[yr]['bkg']:
                histStore['bkg'][ky]=rdataFrames[yr]['bkg'][ky].Histo1D(varToBinMap[var],var,'weight')
            print("   Year : ",yr,' Data')
            for ky in rdataFrames[yr]['data']:
                histStore['data']=rdataFrames[yr]['data'][ky].Histo1D(varToBinMap[var],var,'weight')
            allHistoDict[yr][var]=histStore
            print("Saving into : ",saveBase)
            plotterUtil.plotVariableDistribution(histStore,var,year=yr,lumi=lumiMap[yr],savePrefix=saveBase)




