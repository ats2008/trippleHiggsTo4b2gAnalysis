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

    prefixBase='/home/aravind/cernbox/work/trippleHiggs/hhhTo4b2gamma/genAnalysis/python/analysis/results/plots/analysis/jan19/'
    fileDict={}
    with open('workarea/data/bdtNtuples/v4/filelist.json') as f:
        fileDict=json.load(f)

    yearsToProcess=['2018','2017','2016PreVFP','2016PostVFP','run2']
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
            rdataFrames[yr]['bkg'][bkg]=ROOT.RDataFrame(treeName, fileName)
            print("Registering datset : ",bkg," , ",yr," withh tree",treeName)

    varToBinMap={  
                   'h1bb_mass' : ("","",50,0.0,300.0),
                   'h2bb_mass' : ("","",50,0.0,300.0),
                   'CMS_hgg_mass' :("","",50,100.0,200.0),
                   'pThgg_overMgg' : ("","",20,0.0,4.0),
                   "h1bb_eta" : ("","",30,-3.0,3.0),
                   "h1bb_phi" : ("","",30,-3.2,3.2),
                   "h1bb_pt"  : ("","",100,0.0,600),
                   "h2bb_eta" : ("","",30,-3.0,3.0),
                   "h2bb_phi" : ("","",30,-3.2,3.2),
                   "h2bb_pt"  : ("","",100,0.0,600),
                   "pT_4b"  : ("","",100,0.0,600),
                   "h1_dijetSigmaMOverM"  : ("","",25,0.0,0.25),
                   "h2_dijetSigmaMOverM"  : ("","",25,0.0,0.25),
                    "pTh1leadJ_overMh1"   : ("","",30,0.0,6.0),
                    "pTh1subleadJ_overMh1": ("","",25,0.0,2.5),
                    "pTh2leadJ_overMh2"   : ("","",30,0.0,6.0),
                    "pTh2subleadJ_overMh2": ("","",25,0.0,2.5),
                    "pThgg_overMgg" :  ("","",30,0.0,3.0),
                    "pTleadG_overMgg" : ("","",16,0.0,3.2),
                    "pTsubleadG_overMgg": ("","",16,0.0,1.2),
                    "scalarPtSum4b" : ("","",100,0.0,1000.0),
                    "scalarPtSum4b2g": ("","",100,0.0,1200.0),
                    "scalarPtSumHHH": ("","",100,0.0,1200.0),
                    "sigmaMOverM" : ("","",40,0.0,0.20),
                    "trihiggs_mass" : ("","",140,0.0,2100.0),
                    "trihiggs_pt" : ("","",60,0.0,600.0),
                    "ttH_MET" : ("","",100,0.0,200.0)
                }

    varToProcess=[
           'h1bb_mass' ,
           'h2bb_mass' ,
           'CMS_hgg_mass' ,
           'pThgg_overMgg' ,
           "h1bb_eta" ,
           "h1bb_phi" ,
           "h1bb_pt"  ,
           "h2bb_eta" ,
           "h2bb_phi" ,
           "h2bb_pt"  ,
           "pT_4b"  ,
           "h1_dijetSigmaMOverM"  ,
           "h2_dijetSigmaMOverM"  ,
            "pTh1leadJ_overMh1"   ,
            "pTh1subleadJ_overMh1",
            "pTh2leadJ_overMh2"   ,
            "pTh2subleadJ_overMh2",
            "pThgg_overMgg" ,
            "pTleadG_overMgg" ,
            "pTsubleadG_overMgg",
            "scalarPtSum4b" ,
            "scalarPtSum4b2g",
            "scalarPtSumHHH",
            "sigmaMOverM" ,
            "trihiggs_mass" ,
            "trihiggs_pt" ,
            "ttH_MET" ,
         ]
    allHistoDict={}
    for yr  in yearsToProcess:
        allHistoDict[yr]={}
        print("Procssing for Year  : ",yr)
        for var in varToProcess:
            print("  Procssing for variable  : ",var)
            saveBase=prefixBase+'/'+yr+'/'
            os.system('mkdir -p '+saveBase)
            histStore={'sig':{},'bkg':{},'data':{}}
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




