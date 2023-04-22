import matplotlib.pyplot as plt
import mplhep as hep
import uproot as urt
import ROOT
import json
import numpy as np
import os,argparse

from Util import lumiMap
import Util as utl
import trippleHiggsPlotter as plotterUtil

if __name__=='__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--unblind", help="Unblind the Mgg spectrum", action='store_true' )
    parser.add_argument("-i","--flist", help="File list To Use", default='workarea/data/bdtNtuples/v8p2/filelist.json' )
    parser.add_argument("-o","--dest", help="destination To Use", default='workarea/results/plots/analysis/feb8/' )
    parser.add_argument("-w","--weight", help="weight var to use", default='weight' )
    parser.add_argument('-y',"--year", help="Year",default='2018')
    args=parser.parse_args()

    unblind= args.unblind
    weightVar=args.weight

    prefixBase=args.dest
    prefixBase+='/variables/'
    fileDict={}
    with open(args.flist) as f:
        fileDict=json.load(f)

    yearsToProcess_all=['2018','2017','2016PreVFP','2016PostVFP','run2','2016']
    yearsToProcess=yearsToProcess_all
    
    if 'all' not in args.year:
        yearsToProcess=[]
        for yr in args.year.split(","):
            if yr not in yearsToProcess_all:
                print(yr," not in catalogue. Skipping !! ")
                continue
            yearsToProcess.append(yr)
    
    bkgToProcess=[ 'all','ggBox1Bjet','ggBox2Bjet', 'ggBox','gJet20To40','gJet40ToInf']
    
    rdataFrames={}
    for yr  in yearsToProcess:
        rdataFrames[yr]={'sig':{},'bkg':{},'data':{}}
    
        fileName=fileDict[yr]['sig']['ggHHH']
        treeName = "trees/ggHHH_125_13TeV"
        rdataFrames[yr]['sig']['ggHHH']= ROOT.RDataFrame(treeName, fileName)
        print("Registering datset : ggHHH , ",yr," withh tree",treeName)
        print("\t\tFilename :", fileName)
        ky=list(fileDict[yr]['data'].keys())[0]
        fileName=fileDict[yr]['data'][ky]
        treeName = "trees/Data_13TeV_TrippleHTag_0"
        if unblind:
            rdataFrames[yr]['data']['data']= ROOT.RDataFrame(treeName, fileName)
            prefixBase+='/unblinded/'
        else:
            rdataFrames[yr]['data']['data']= ROOT.RDataFrame(treeName, fileName).Filter(' CMS_hgg_mass < 115.0  || CMS_hgg_mass >135.0')
        print("Registering datset : data , ",yr," withh tree",treeName)
        print("\t\tFilename :", fileName)
    
        rdataFrames[yr]['bkg']={}
        treeName = "trees/bkg_13TeV_TrippleHTag_0"
        for bkg in  fileDict[yr]['bkg']:
            if bkg=='bkg':
                continue
            if ("all" not in bkgToProcess ) and (bkg not in bkgToProcess) :
                print()
                print("FILE NOT FOUNG FOR : ",yr," background ",bkg)
                print()
                continue
            fileName = fileDict[yr]['bkg'][bkg]
            rdataFrames[yr]['bkg'][bkg]=ROOT.RDataFrame(treeName, fileName).Filter(' CMS_hgg_mass < 115.0  || CMS_hgg_mass >135.0')
            print("Registering datset : ",bkg," , ",yr," withh tree",treeName)
            print("\t\tFilename :", fileName)

    varToBinMap={  
                   'h1bb_mass' : ("","",30,0.0,300.0),
                   'h2bb_mass' : ("","",30,0.0,300.0),
                   "diphoton_pt"  : ("","",20,0.0,400),
                   'pThgg_overMgg' : ("","",20,0.0,4.0),
                   "h1bb_eta" : ("","",30,-3.0,3.0),
                   "h1bb_phi" : ("","",30,-3.2,3.2),
                   "h1bb_pt"  : ("","",20,0.0,400),
                   "h2bb_eta" : ("","",30,-3.0,3.0),
                   "h2bb_phi" : ("","",30,-3.2,3.2),
                   "h2bb_pt"  : ("","",20,0.0,400),
                   "pT_4b"  : ("","",20,0.0,400),
                   "leadingPhotonSigOverE"  : ("","",40,0.0,0.20),
                   "subleadingPhotonSigOverE"  : ("","",40,0.0,0.20),
                   "h1_dijetSigmaMOverM"  : ("","",25,0.0,0.25),
                   "h2_dijetSigmaMOverM"  : ("","",25,0.0,0.25),
                   "pT_h1leadJ"  : ("","",20,0.0,400),
                   "pT_h1subleadJ"  : ("","",20,0.0,400),
                   "pT_h2leadJ"  : ("","",20,0.0,400),
                   "pT_h2subleadJ"  : ("","",20,0.0,400),
                   "pTh1subleadJ_overMh1": ("","",25,0.0,2.5),
                   "pTh1leadJ_overMh1":("","",25,2.0,2.5),
                   "pTh2subleadJ_overMh2": ("","",25,0.0,2.5),
                   "pTh2leadJ_overMh2":("","",25,2.0,2.5),
                   "pThgg_overMgg" :  ("","",30,0.0,3.0),
                   "pTleadG_overMgg" : ("","",16,0.0,3.2),
                   "pTsubleadG_overMgg": ("","",16,0.0,1.2),
                   "scalarPtSum4b" : ("","",20,0.0,1000.0),
                   "scalarPtSum4b2g": ("","",20,0.0,1200.0),
                   "scalarPtSumHHH": ("","",20,0.0,1200.0),
                   "sigmaMOverM" : ("","",40,0.0,0.20),
                   "trihiggs_mass" : ("","",20,0.0,2100.0),
                   "trihiggs_pt" : ("","",20,0.0,400.0),
                   "ttH_MET" : ("","",10,0.0,200.0),
                   "rho" : ("","",80,0.0,80.0),
                }
    varToExpressionMap={

    }

    varToProcess=[
         'h1bb_mass' ,
         'h2bb_mass' ,
         "diphoton_pt"  ,
         'pThgg_overMgg' ,
         "h1bb_eta" ,
         "h1bb_phi" ,
         "h1bb_pt"  ,
         "h2bb_eta" ,
         "h2bb_phi" ,
         "h2bb_pt"  ,
         "pT_4b"  ,
         "leadingPhotonSigOverE"  ,
         "subleadingPhotonSigOverE"  ,
         "h1_dijetSigmaMOverM"  ,
         "h2_dijetSigmaMOverM"  ,
         "pT_h1leadJ"  ,
         "pT_h1subleadJ"  ,
         "pT_h2leadJ"  ,
         "pT_h2subleadJ"  ,
         "pTh1subleadJ_overMh1",
         "pTh1leadJ_overMh1",
         "pTh2subleadJ_overMh2",
         "pTh2leadJ_overMh2",
         "pThgg_overMgg" ,
         "pTleadG_overMgg" ,
         "pTsubleadG_overMgg",
         "scalarPtSum4b" ,
         "scalarPtSum4b2g",
         "scalarPtSumHHH",
         "sigmaMOverM" ,
         "trihiggs_mass" ,
         "trihiggs_pt" ,
         "ttH_MET" 
         ]
    varToProcess=["rho"]
    allHistoDict={}
    for yr  in yearsToProcess:
        allHistoDict[yr]={}
        print("Procssing for Year  : ",yr)
        for var in varToProcess:
            print("  Procssing for variable  : ",var)
            saveBase=prefixBase+'/'+yr+'/composit/'
            os.system('mkdir -p '+saveBase)
            histStore={'sig':{},'bkg':{},'data':{}}
            print("   Year : ",yr,' Signal')
            for ky in rdataFrames[yr]['sig']:
                histStore['sig'][ky]=rdataFrames[yr]['sig'][ky].Histo1D(varToBinMap[var],var,weightVar)
            print("   Year : ",yr,' Background')
            for ky in utl.backgroundStackList:
                if ky not in rdataFrames[yr]['bkg']:
                    continue
                histStore['bkg'][ky]=rdataFrames[yr]['bkg'][ky].Histo1D(varToBinMap[var],var,weightVar)
            print("   Year : ",yr,' Data')
            for ky in rdataFrames[yr]['data']:
                histStore['data']=rdataFrames[yr]['data'][ky].Histo1D(varToBinMap[var],var,weightVar)
            allHistoDict[yr][var]=histStore
            print("Saving into : ",saveBase)
            plotterUtil.plotVariableDistribution(histStore,var,year=yr,lumi=lumiMap[yr],savePrefix=saveBase)




