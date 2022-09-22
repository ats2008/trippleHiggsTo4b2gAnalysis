from PlotUtil import *
from Util import *
import os
from collections import OrderedDict

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
    
    for preReweight in [False]:
        fileDict=OrderedDict()
        fileDict[ "signal"          ] =   "/home/aravind/workarea/trippleHiggs/batch/cutBased_ntupleSkim/ggHHH_M125_13TeV.root"
        fileDict[ "Data"          ] =   "/home/aravind/workarea/trippleHiggs/batch/cutBased_ntupleSkim/EGamma_alesauva-UL2018_0-10_6_4-v0-Run2018-12Nov2019_UL2018-v4.root"
           
        
        binsToPrint=['total',"*"]

        histStore={}
        fileStore={}
        for tag in fileDict:
            print("="*40)
            print(tag,fileDict[tag].split("/")[-1])
            fileStore=ROOT.TFile(fileDict[tag],'READ')
            histStore=getTheObjectsFromFile(fileStore)
            nHist=histStore['tagsDumper']['trees']['sumWeights'].Clone()
            counts=getBinVsContent(nHist)
            print(counts)
            tw=0
            t=0
            print(counts)
            if "*" in binsToPrint:
                for head in counts:
                    print("\t",head, " : ",counts[head])
            else:
                for head in binsToPrint:
                    if head in counts:
                        print("\t",head, " : ",counts[head])
            
            fileStore.Close()
            print()
            print()

