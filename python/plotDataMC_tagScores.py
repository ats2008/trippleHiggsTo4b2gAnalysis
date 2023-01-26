import matplotlib.pyplot as plt
import mplhep as hep
import uproot as urt
import ROOT
import json,os,argparse
import Util as utl
import numpy as np

from Util import lumiMap
import trippleHiggsPlotter as plotterUtil

if __name__=='__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--unblind", help="Unblind the Mgg spectrum", action='store_true' )
    args=parser.parse_args()

    unblind= args.unblind
    
    prefixBase='/home/aravind/cernbox/work/trippleHiggs/hhhTo4b2gamma/genAnalysis/python/analysis/results/plots/analysis/jan25/variables_v0Scaled/tagScores/'
    fileDict={}
    #with open('workarea/data/bdtNtuples/filelistToUse.json') as f:
    with open('workarea/data/analysisNtuples/v1p1/filelist.json') as f:
        fileDict=json.load(f)
    yearsToProcess=['2018','2017','2016PreVFP','2016PostVFP','run2','2016']
    yearsToProcess=['2018','run2']
    bkgToProcess=[ 'ggBox1Bjet','ggBox2Bjet', 'ggBox','gJet20To40','gJet40ToInf']
    #bkgToProcess=[ 'bkg']

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
                   'sumScore_4j' : ("","",20,0.0,4.0) ,
                   'sumScore_3j' : ("","",15,0.0,3.0) ,
                   'quadjet_0_deepJetScore' : ("","",16,0.0,1.0),
                   'quadjet_1_deepJetScore' : ("","",16,0.0,1.0),
                   'quadjet_2_deepJetScore' : ("","",16,0.0,1.0),
                   'quadjet_3_deepJetScore' : ("","",16,0.0,1.0)
                }

    varToProcess=[
                    'sumScore_4j',
                    'quadjet_0_deepJetScore',
                    'quadjet_1_deepJetScore',
                    'quadjet_2_deepJetScore',
                    'quadjet_3_deepJetScore',
                    'sumScore_3j'
                ]

    allHistoDict={}
    for yr  in yearsToProcess:
        allHistoDict[yr]={}
        saveBase=prefixBase+'/'+yr+'/'
        os.system('mkdir -p '+saveBase)
        print(" Procssing for year  : ",yr)
        for var in varToProcess:
            histStore={'sig':{},'bkg':{},'data':{}}
            print("    Procssing for variable  : ",var)
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
            plotterUtil.plotVariableDistribution(histStore,var,year=yr,lumi=lumiMap[yr],savePrefix=saveBase,lumiScale=False)



