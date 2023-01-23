from __future__ import print_function
import ROOT ,copy
import numpy as np

class mcScaler:
    def __init__(self):
        self.dataHist={}
        self.mcHist={}
        self.binEdges=[]
        self.scaleFactors=[]
        self.def_scaleFactor=0
        self.underflow_x={}
        self.overflow_x={}
        self.underflow_scaleFactor={}
        self.overflow_scaleFactor={}
        self.scaleFactorHist={}
        self.scaleFactorHistSet={}
        self.File=None
    
    def setSFHist(self,hist,cat='def'):
        self.scaleFactorHist[cat]=hist.Clone()
        self.underflow_x[cat]= self.scaleFactorHist[cat].GetBinLowEdge(1)
        self.overflow_x[cat]= self.scaleFactorHist[cat].GetBinLowEdge(self.scaleFactorHist[cat].GetNbinsX()+1)
        self.underflow_scaleFactor[cat] = self.scaleFactorHist[cat].GetBinContent(1)
        self.overflow_scaleFactor[cat] = self.scaleFactorHist[cat].GetBinContent(self.scaleFactorHist[cat].GetNbinsX())
        print("  Setting scale factor for category " , cat  )
        print("     Scale factors defined for x = [ ",self.underflow_x[cat],self.overflow_x[cat]," ]")
        print("     Under flow Scale factor     =  ",self.underflow_scaleFactor[cat])
        print("     Over flow Scale factor      = ",self.overflow_scaleFactor[cat])

        self.scaleFactorHistSet[cat]=True
    
    def setSFHistFromFile(self,fname,histNames={}):
        self.File=ROOT.TFile(fname,"READ")
        print(histNames)
        for cat in histNames:
            print("CAT : ",cat," || Setting SF hist ",histNames[cat]," from ",fname)
            hist=self.File.Get(histNames[cat])
            self.setSFHist(hist,cat)
    
    def getSFForX(self,x,cat='def'):
        if cat not in self.scaleFactorHistSet:
            print("\tWARNING ! : scaleFactor Hist not set , returning -1.0  [ cat : ",cat," available cats : ",self.scaleFactorHist.keys(),"]")
            return -1.0
        if x < self.underflow_x[cat]:
            return self.underflow_scaleFactor[cat]
        if x > self.overflow_x[cat]:
            return self.overflow_scaleFactor[cat]
        bid=self.scaleFactorHist[cat].FindBin(x)
        #print("bid = ",bid," for x ",x, " Bin edges --> ",
        #         self.scaleFactorHist.GetBinLowEdge(bid), 
        #         self.scaleFactorHist.GetBinLowEdge(bid)+self.scaleFactorHist.GetBinWidth(bid))
        return self.scaleFactorHist[cat].GetBinContent(bid)
    def getScaleFactorHist(self,cat='def'):
        return self.scaleFactorHist[cat]

    def setDataHist(self,hist,cat='def'):
        self.dataHist[cat]=hist.Clone()
        self.dataHist[cat].SetName(cat+'_dataHist')
    def setMCHist(self,hist,cat='def'):
        self.mcHist[cat]=hist
        self.mcHist[cat].SetName(cat+'_mcHist')
    
    def getDataHist(self,cat='def'):
        return self.dataHist[cat]
    def getMCHist(self,cat='def'):
        return self.mcHist[cat]

def getPtDependentScaleFactor(name="scaleFactor",
                              dataHist=None,
                              mcHists=None,
                              mcHistsToSubstract=None,
                              binEdges=None,
                              def_scaleFactor=0.0):
    mcHistSum=mcHists[0].Clone()
    mcHistSum.Reset()
    scaleFactorHist=ROOT.TH1F("","",1,0,100.0)
    if binEdges!=None:
#         binEdges.append(binEdges[-1]) # hck for fix, dont know why last bin not taken
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
    for i in range(1,dataHist.GetNbinsX()+1):
        x=dataHist.GetBinCenter(i)
        dValue=dataHist.GetBinContent(i)
        mcValue=mcHistSum.GetBinContent(i)
        dataHTmp.Fill(x,dValue)
        mcHTmp.Fill(x,mcValue)
        
    for i in range(1,dataHTmp.GetNbinsX()+1):
        dValue =dataHTmp.GetBinContent(i)
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
        print(" [ ",scaleFactorHist.GetBinLowEdge(i)," ] ", " -> " ,scaleFactorHist.GetBinContent(i) )
    scaler.binEdges.append(dataHist.GetBinLowEdge(dataHist.GetNbinsX()+1))
    scaler.setSFHist(scaleFactorHist)
    scaler.def_scaleFactor = 0.0
    scaler.setMCHist(mcHTmp)
    scaler.setDataHist(dataHTmp)
    scaler.setSFHist(scaleFactorHist)

    return scaler

