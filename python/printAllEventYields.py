from __future__ import print_function
import Util as util
import os
import ROOT
from collections import OrderedDict
import argparse

def getBinVsContent(hist):
    xAxis=hist.GetXaxis()
    nBins=hist.GetNbinsX()
    cutFlow={}
    
    for i in range(nBins):
        binName=xAxis.GetBinLabel(i)
        binContent=hist.GetBinContent(i)
        if binName=="":
            continue
        cutFlow[binName]=binContent#np.round(float(binContent),5)
    return cutFlow


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", help="File to anlyze",default=None)
    parser.add_argument("--hname", help="File to anlyze",default="sumWeights")
    args = parser.parse_args()

    for preReweight in [False]:
        fileDict=OrderedDict()
        if args.file:
            print("File  : ",args.file)
            fileDict["file"] = args.file
        #fileDict[ "signal"          ] =   "out_ggHHH_AnalysisV1_3p1_0.root"
        #fileDict[ "Data"          ] =   "/home/aravind/workarea/trippleHiggs/batch/cutBased_ntupleSkim/EGamma_alesauva-UL2018_0-10_6_4-v0-Run2018-12Nov2019_UL2018-v4.root"
           
        
        binsToPrint=['total',"*"]

        histStore={}
        fileStore={}
        for tag in fileDict:
            
            print("="*20,"   ",tag,"   ","="*20)
            
            fileStore=ROOT.TFile(fileDict[tag],'READ')
            histStore=util.getTheObjectsFromFile(fileStore)
            nHist=histStore['tagsDumper']['trees'][args.hname].Clone()
            counts=getBinVsContent(nHist)
            
            print()

            tw=0
            t=0
            if "*" in binsToPrint:
                for head in counts:
                    print("\t",head, " : ",counts[head])
            else:
                for head in binsToPrint:
                    if head in counts:
                        print("\t",head, " : ",counts[head])
            
            fileStore.Close()
            
            print()

