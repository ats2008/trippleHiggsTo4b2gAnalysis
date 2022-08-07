from trippleHiggsUtils import *
from Util  import *

def getPtDependentScaleFactor(name="scaleFactor",dataHist=None,mcHists=None,mcHistsToSubstract=None,binEdges=None,def_scaleFactor=0.0):
    mcHistSum=mcHists[0].Clone()
    mcHistSum.Reset()
    scaleFactorHist=ROOT.TH1F("","",1,0,100.0)
    if binEdges!=None:
        binEdges.append(binEdges[-1]+200.0) # hck for fix, dont know why last bin not taken
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
    for i in range(1,dataHist.GetNbinsX()+1):
        x=dataHist.GetBinCenter(i)
        mcValue=mcHistSum.GetBinContent(i)
        dataHist.Fill(x,-1*mcValue)
        
    for i in range(1,dataHTmp.GetNbinsX()+1):
        dValue=dataHTmp.GetBinContent(i)
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
    scaler.binEdges.append(dataHist.GetBinLowEdge(dataHist.GetNbinsX()+1))
    scaler.mcHist   = mcHTmp
    scaler.mcHist.SetName(name+"_mcHist")
    scaler.dataHist = dataHTmp
    scaler.dataHist.SetName(name+"_dataHist")
    scaler.setSFHist(scaleFactorHist)
    scaler.def_scaleFactor = 0.0
    return scaler


if __name__ == "__main__":
    fListDict={
        'data2018'       : 'results/nonReweited/data2018.root',
        #'signal'         : 'results/nonReweited/signal_hhhTo4b2g.root',
        'ggM80Inc'       : 'results/nonReweited/diphotonInclusive.root',
        'ggM80Jbox1bjet' : 'results/nonReweited/diphotonJetBox1bjet.root',
        'ggM80Jbox2bjet' : 'results/nonReweited/diphotonJetBox2bjet.root',
        'ttgJets'        : 'results/nonReweited/TTGJets.root',
        'ttgg'           : 'results/nonReweited/TTGG.root',
        'tGJets'         : 'results/nonReweited/TGJets.root',
        'gJet20To40'     : 'results/nonReweited/photonInclusive20To40.root',
        'gJet40ToInf'    : 'results/nonReweited/photonInclusive40ToInf.root'
    }
    
    histStore={}
    fileStore={}
    for tag in fListDict:
        print("Adding ",tag,"  : ",fListDict[tag])
        if tag in fileStore:
            fileStore[tag].Close()
        fileStore[tag]=ROOT.TFile(fListDict[tag],'READ')
        histStore[tag]=getTheObjectsFromFile(fileStore[tag])

    rebinC=1
    
    dataHist=histStore['data2018']['controlCands']['hgg']['pT'].Clone()
    dataHist.Rebin(rebinC)
    mcHists=[]
    mcHistsToSubstract=[]
    #for tag in ['ggM80Inc', 'ggM80Jbox1bjet','ggM80Jbox2bjet','gJet20To40','gJet40ToInf']:
    for tag in ['ggM80Inc', 'ggM80Jbox1bjet','ggM80Jbox2bjet']:
        mcHists.append(histStore[tag]['controlCands']['hgg']['pT'].Clone())
        mcHists[-1].Rebin(rebinC)
    for tag in ['ttgJets', 'ttgg', 'tGJets']:
        mcHistsToSubstract.append(histStore[tag]['controlCands']['hgg']['pT'].Clone())
        mcHistsToSubstract[-1].Rebin(rebinC)
    binEdges=[0,25.0,50.0,75.0,100.0,150.0,200.0,300.0,500.0]
    pTDependentScaler=getPtDependentScaleFactor("pTDependent",dataHist,mcHists,mcHistsToSubstract,binEdges,def_scaleFactor=0.0)    
    
    d=ROOT.TCanvas()

    legend = ROOT.TLegend(0.65,0.78,0.8,0.88)
    legend.SetTextFont(42)
    legend.SetFillColor(0)
    
    h=pTDependentScaler.mcHist.Clone()
    g=pTDependentScaler.dataHist.Clone()
    h.SetAxisRange(0.0,1e4, "Y")
    h.SetAxisRange(0.0,900.0, "X")
    h.SetMarkerColor(8)
    h.SetLineColor(4)
    legend.AddEntry(h,"MC","pe")
    h.Draw()
    g.SetLineColor(2)
    legend.AddEntry(g,"Data","pe")
    g.Draw("same pe")
    legend.Draw()
    d.Draw()
    d.SaveAs("dataMC.png")
    k=pTDependentScaler.scaleFactorHist.Clone()
    k.SetAxisRange(0.0,700.0, "X")
    k.Draw()
    d.Draw()
    d.SaveAs("Scalefactor.png")

    f=ROOT.TFile("pTScaleFactorFile.root","RECREATE")
    pTDependentScaler.dataHist.Write()
    pTDependentScaler.mcHist.Write()
    pTDependentScaler.scaleFactorHist.Write()
    f.Close()
