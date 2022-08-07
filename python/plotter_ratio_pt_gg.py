from PlotUtil import *
from Util import *
import os


if __name__ == "__main__":
    candDict={}
    candDict['controlCands'] ='results/plots/analysis/controlRegion'
    candDict['trippleHCands']='results/plots/analysis/signalRegion'
    candDict['validation_CR']='results/plots/analysis/validationCR'
    candDict['validation_SR']='results/plots/analysis/validationSR'
  
    for preReweight in  [True,False]:
        fListDict={
            'data2018'       : 'results/reweighted/data2018.root',
#           'signal'         : 'results/reweighted/signal_hhhTo4b2g.root',
            'ggM80Inc'       : 'results/reweighted/diphotonInclusive.root',
            'ggM80Jbox1bjet' : 'results/reweighted/diphotonJetBox1bjet.root',
            'ggM80Jbox2bjet' : 'results/reweighted/diphotonJetBox2bjet.root',
            'ttgJets'        : 'results/reweighted/TTGJets.root',
            'ttgg'           : 'results/reweighted/TTGG.root',
            'tGJets'         : 'results/reweighted/TGJets.root',
            'gJet20To40'     : 'results/reweighted/photonInclusive20To40.root',
            'gJet40ToInf'    : 'results/reweighted/photonInclusive40ToInf.root'
        }
        if preReweight:
            fListDict={
                'data2018'       : 'results/nonReweited/data2018.root',
 #              'signal'         : 'results/nonReweited/signal_hhhTo4b2g.root',
                'ggM80Inc'       : 'results/nonReweited/diphotonInclusive.root',
                'ggM80Jbox1bjet' : 'results/nonReweited/diphotonJetBox1bjet.root',
                'ggM80Jbox2bjet' : 'results/nonReweited/diphotonJetBox2bjet.root',
                'ttgJets'        : 'results/nonReweited/TTGJets.root',
                'ttgg'           : 'results/nonReweited/TTGG.root',
                'tGJets'         : 'results/nonReweited/TGJets.root',
                'gJet20To40'     : 'results/nonReweited/photonInclusive20To40.root',
                'gJet40ToInf'    : 'results/nonReweited/photonInclusive40ToInf.root'
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
        
        histSum=getSumHistDicts([histStore['ttgJets'],histStore['ttgg'],histStore['tGJets']])
        histStore['ttX']=histSum
        legendDict['ttX']='tt/t + X'
        
        histSum=getSumHistDicts([histStore['gJet20To40'],histStore['gJet40ToInf']])
        histStore['gJets']=histSum
        legendDict['gJets']='#gamma+Jets'

        for cands in candDict:
            pltPrefix=candDict[cands]
            if not os.path.exists(pltPrefix):
                print("Making folder :" ,pltPrefix)
            os.system("mkdir -p "+pltPrefix)    
        

            plots=[]
            cPars=getCanvasParams('GEN')
            rBin=5
            sc=1.0
            plot=getDefaultPlot(prefix=pltPrefix,name='XXX',cPars=cPars) 
            plot.legendPosition = (0.6,0.74,0.9,0.90)
            plot.descPosition   = (0.6,0.85)
            plot.desc =[]
            plot.yTitle = "# / 2 GeV"  
            plot.xTitle = "p_{T}^{#gamma#gamma} [GeV]"
            plot.xRange = (0.0,300.0) 
            plot.yRange = (1.0,10000.0) 
            plot.yRatioRange = (0.0,2.0) 
            if preReweight:
                plot.yRatioRange = (0.0,10.0) 
            plots.append(plot)

            canvas = ROOT.TCanvas("c556", "", 700,800)
            canvas.SetGrid()
            plot.setPlotStyle()
            canvas.SetBottomMargin(0.1)
            
            pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
            pad1.SetBottomMargin(0)  # joins upper and lower plot
            pad1.SetGridx()
            pad1.Draw()

            #canvas.cd()  # returns to main canvas before defining pad2
            pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
            pad2.SetTopMargin(0)  # joins upper and lower plot
            pad2.SetBottomMargin(0.3)
            pad2.SetGridx()
            pad2.Draw()

            pad1.cd()
            h=ROOT.TH1F("xx",'xx',1,plot.xRange[0],plot.xRange[1])
            h.SetAxisRange(plot.yRange[0], plot.yRange[1], "Y")
            h.SetXTitle( plot.xTitle )
            h.SetYTitle( plot.yTitle )
            h.GetYaxis().SetTitleOffset(1)
            h.Draw()
            h.Draw()
            pad2.cd()
            hR=h.Clone()
            hR.SetAxisRange(plot.yRatioRange[0], plot.yRatioRange[1], "Y")
            hR.SetYTitle( "ratio" )
            y = hR.GetYaxis()
            y.SetTitle("ratio")
            y.SetNdivisions(505)
            y.SetTitleSize(20)
            y.SetTitleFont(43)
            y.SetTitleOffset(1)
            y.SetLabelFont(43)
            y.SetLabelSize(15)

            # Adjust x-axis settings
            x = hR.GetXaxis()
            x.SetTitleSize(20)
            x.SetTitleFont(43)
            x.SetTitleOffset(4.0)
            x.SetLabelFont(43)
            x.SetLabelSize(15)

            hR.Draw()
            l=ROOT.TLine(plot.xRange[0],1.0,plot.xRange[1],1.0);
            l.SetLineColor(2);
            l.Draw("same");

            CMSbox=plot.getCMSBox(pad1)  
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
            pad1.cd()
            CMSbox.Draw()
            extraTextBox.Draw()
            lumibox.Draw()
            
            legend = ROOT.TLegend(plot.legendPosition[0],plot.legendPosition[1],plot.legendPosition[2],plot.legendPosition[3])
            legend.SetTextFont(42)
            legend.SetFillColor(0)
            
            hs = ROOT.THStack("hs","Stacked 1D histograms");
            iCol=2

            histMC=None
            #for ky in ['ttX','ggJetsM80']:
            for ky in ['ttX','ggJetsM80','gJets']:
                hist=histStore[ky][cands]['hgg']['pT'].Clone()
                hist=hist.Rebin(rBin)
                
                if histMC==None:
                    histMC=hist.Clone()
                else:
                    histMC=histMC+hist
                
                legend.SetTextSize(0.025) 
                hist.Scale(sc)
                hist.SetFillColor(col[iCol]);       
                lege=legendDict[ky]
                legend.AddEntry(hist, lege, "f")

                iCol+=1
                print("tag  : ",ky," : ",hist.Integral()," ",hist.GetBinWidth(1))
                hs.Add(hist)
            
                
            dataHist=histStore['data2018'][cands]['hgg']['pT'].Clone()
            dataHist=dataHist.RebinX(rBin)
            legend.AddEntry(dataHist,"Data 2018","pe")
            print("tag  : data : ",dataHist.Integral()," bw : ",dataHist.GetBinWidth(1))
            dataHist.SetLineColor(1)
            dataHist.SetLineWidth(3)
            
            pad1.cd()
            pad1.SetLogy()
            hs.Draw("same HIST")
            dataHist.Draw("same HIST pe")
            legend.Draw("same")

            pad2.cd()
            histRatio=dataHist.Clone()
            histRatio.Divide(histMC)
            print("tag  : ratio : ",histRatio.Integral()," bw : ",histRatio.GetBinWidth(1))
            histRatio.SetLineWidth(2)
            histRatio.Draw("same HIST pe")
            
            axis = ROOT.TGaxis(plot.xRange[0],plot.yRange[0],plot.yRange[0], plot.yRange[1], plot.yRange[0], plot.yRange[1],110, "")
            axis.SetLabelFont(43)
            axis.SetLabelSize(15)
            axis.Draw()
            
            for i in range(len(selbox)):
                selbox[i].Draw()
            CMSbox.Draw()
            extraTextBox.Draw()
            lumibox.Draw()

            #canvas.SetLogy()
            #canvas.Draw()
            oFname=pltPrefix+'/pT_gg_ratio.png'
            if preReweight:
                oFname=pltPrefix+'/pT_preReweighting_gg_ratio.png'
            canvas.SaveAs(oFname)


        for tag in fListDict:
            fileStore[tag].Close()
