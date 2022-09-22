from PlotUtil import *
from Util import *
import os


if __name__ == "__main__":
    cands='controlCands'
    pltPrefix='results/plots/analysis/controlRegion'
    
    if not os.path.exists(pltPrefix):
        print("Making folder :" ,pltPrefix)
    os.system("mkdir -p "+pltPrefix)    
    
    fListDict={
        'data2018'       : 'results/data2018.root',
        'signal'         : 'results/signal_hhhTo4b2g.root',
        'ggM80Inc'       : 'results/diphotonInclusive.root',
        'ggM80Jbox1bjet' : 'results/diphotonJetBox1bjet.root',
        'ggM80Jbox2bjet' : 'results/diphotonJetBox2bjet.root',
        'ttgJets'        : 'results/TTGJets.root',
        'ttgg'           : 'results/TTGG.root',
        'tGJets'         : 'results/TGJets.root',
        'gJet20To40'     : 'results/photonInclusive20To40.root',
        'gJet40ToInf'    : 'results/photonInclusive40ToInf.root'
    }
    legendDict={
        'data2018'       :'data2018'        ,
        'signal'       :'hhh#rightarrow 4b2#gamma'        ,
        'ggM80Inc'       :'ggM80Inc'        ,
        'ggM80Jbox1bjet' :'ggM80Jbox1bjet'  ,
        'ggM80Jbox2bjet' :'ggM80Jbox2bjet'  ,
        'ttgJets'        :'ttgJets'         ,
        'ttgg'           :'ttgg'          ,
        'tGJets'        : 'tGJets'
    }
    
    histStore={}
    fileStore={}
    for tag in fListDict:
        print("Adding ",tag)
        if tag in fileStore:
            fileStore[tag].Close()
        fileStore[tag]=ROOT.TFile(fListDict[tag],'READ')
        histStore[tag]=getTheObjectsFromFile(fileStore[tag])

    histSum=getSumHistDicts([histStore['ggM80Jbox1bjet'],histStore['ggM80Jbox2bjet'],histStore['ggM80Inc']])
    histStore['ggJetsM80']=histSum
    legendDict['ggJetsM80']='#gamma#gamma+Jets'
    
    histSum=getSumHistDicts([histStore['gJet20To40'],histStore['gJet40ToInf']])
    histStore['gJets']=histSum
    legendDict['gJets']='#gamma+Jets'



    plots=[]
    cPars=getCanvasParams('GEN')
    rBin=2
    sc=1.0
    plot=getDefaultPlot(prefix=pltPrefix,name='XXX',cPars=cPars) 
    plot.legendPosition = (0.6,0.74,0.9,0.90)
    plot.descPosition   = (0.6,0.85)
    plot.desc =[]
    plot.yTitle = "# / 2 GeV"  
    plot.xTitle = "p_{T}^{#gamma#gamma} [GeV]"
    plot.xRange = (100,180.0) 
    plot.yRange = (1.0,10000.0) 
    plots.append(plot)

    canvas = ROOT.TCanvas("c556", "", 700,700)
    canvas.SetGrid()
    plot.setPlotStyle()
    
    h=ROOT.TH1F("xx",'xx',1,plot.xRange[0],plot.xRange[1])
    h.SetAxisRange(plot.yRange[0], plot.yRange[1], "Y")
    h.SetXTitle( plot.xTitle )
    h.SetYTitle( plot.yTitle )
    h.Draw()
    h.Draw()

    CMSbox=plot.getCMSBox()  
    extraTextBox=plot.getExtraTextBox()
    lumibox=plot.getLumiBox()

    selbox=[]
    n=0
    for inx in range(len(plot.desc)):
        selbox.append(ROOT.TLatex  (plot.descPosition[0], plot.descPosition[1] -n*0.04, plot.desc[n]))
        selbox[-1].SetNDC()
        selbox[-1].SetTextSize(0.035)
        selbox[-1].SetTextFont(12)
        selbox[-1].SetTextAlign(13)
        n+=1
    for i in range(len(selbox)):
        selbox[i].Draw()
    CMSbox.Draw()
    extraTextBox.Draw()
    lumibox.Draw()
    
    legend = ROOT.TLegend(plot.legendPosition[0],plot.legendPosition[1],plot.legendPosition[2],plot.legendPosition[3])
    legend.SetTextFont(42)
    legend.SetFillColor(0)
    
    hs = ROOT.THStack("hs","Stacked 1D histograms");
    iCol=2

    
    for ky in ['ttgJets', 'ttgg','tGJets','gJets','ggJetsM80']:
        hist=histStore[ky][cands]['hgg']['pT'].Clone()
        legend.SetTextSize(0.025) 
        hist.RebinX(rBin)
        hist.Scale(sc)
        hist.SetFillColor(col[iCol]);       
        lege=legendDict[ky]
        legend.AddEntry(hist, lege, "f")
        
        iCol+=1
        print("tag  : ",ky," : ",hist.Integral()," ",hist.GetBinWidth(1))
        hs.Add(hist)
    
    hs.Draw("same HIST")
        
    hist=histStore['data2018'][cands]['hgg']['pT'].Clone()
    legend.AddEntry(hist,"Data 2018","pe")
    hist.RebinX(rBin)
    print("tag  : data : ",hist.Integral()," bw : ",hist.GetBinWidth(1))
    hist.SetLineWidth(2)
    hist.Draw("same HIST")

    legend.Draw("same")
    
    for i in range(len(selbox)):
        selbox[i].Draw()
    CMSbox.Draw()
    extraTextBox.Draw()
    lumibox.Draw()

    canvas.SetLogy()
    canvas.Draw()
    canvas.SaveAs(pltPrefix+'/pt_gg.png')


    for tag in fListDict:
        fileStore[tag].Close()
