import matplotlib.pyplot as plt
import mplhep as hep
import uproot as urt
import ROOT
import json
import numpy as np
import os,argparse
from Util import lumiMap
import trippleHiggsPlotter as plotterUtil

if __name__=='__main__':
    

    parser = argparse.ArgumentParser()
    parser.add_argument("--unblind", help="Unblind the Mgg spectrum", action='store_true' )
    parser.add_argument("--lumiScale", help="Scale the weights with Lumi", action='store_true' )
    parser.add_argument("-i","--flist", help="File list To Use", default='workarea/data/bdtNtuples/v8p2/filelist.json' )
    parser.add_argument("-o","--dest", help="destination To Use", default='workarea/results/plots/analysis/feb8/' )
    parser.add_argument("-w","--weight", help="weight var to use", default='weight' )
    parser.add_argument('-y',"--year", help="Year",default='2018')
    args=parser.parse_args()

    unblind= args.unblind
    inputfile= args.flist
    prefixBase=args.dest
    weightVar=args.weight
    prefixBase+='/variables/'
    fileDict={}

    with open(inputfile) as f:
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
    bkgToProcess=['all', 'ggBox1Bjet','ggBox2Bjet', 'ggBox','gJet20To40','gJet40ToInf']

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
            #rdataFrames[yr]['bkg'][bkg]=ROOT.RDataFrame(treeName, fileName)
            rdataFrames[yr]['bkg'][bkg]=ROOT.RDataFrame(treeName, fileName).Filter(' CMS_hgg_mass < 115.0  || CMS_hgg_mass >135.0')
            print("Registering datset : ",bkg," , ",yr," withh tree",treeName)
            print("\t\tFilename :", fileName)

    varToBinMap={  
                   'hhhNR_mva_2023_v0' :("","",20,0.0,1.0),
                }

    varToProcess=[
           'hhhNR_mva_2023_v0' 
         ]
    allHistoDict={}
    for yr  in yearsToProcess:
        allHistoDict[yr]={}
        print("Procssing for Year  : ",yr)
        for var in varToProcess:
            print("  Procssing for variable  : ",var)
            saveBase=prefixBase+'/'+yr+'/mva_score/'
            os.system('mkdir -p '+saveBase)
            histStore={'sig':{},'bkg':{},'data':{}}
            print("   Year : ",yr,' Signal')
            for ky in rdataFrames[yr]['sig']:
                histStore['sig'][ky]=rdataFrames[yr]['sig'][ky].Histo1D(varToBinMap[var],var,weightVar)
            print("   Year : ",yr,' Background')
            for ky in rdataFrames[yr]['bkg']:
                if ky=='bkg':
                    continue
                histStore['bkg'][ky]=rdataFrames[yr]['bkg'][ky].Histo1D(varToBinMap[var],var,weightVar)
            print("   Year : ",yr,' Data')
            for ky in rdataFrames[yr]['data']:
                histStore['data']=rdataFrames[yr]['data'][ky].Histo1D(varToBinMap[var],var,weightVar)
            allHistoDict[yr][var]=histStore
            print("Saving into : ",saveBase)
            plotterUtil.plotVariableDistribution(histStore,var,year=yr,lumi=lumiMap[yr],savePrefix=saveBase,lumiScale=args.lumiScale)

