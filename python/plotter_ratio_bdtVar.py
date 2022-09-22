from PlotUtil import *
from Util import *
import os

col=[ 0,
      ROOT.kBlack,
      ROOT.kMagenta+1,
      ROOT.kBlue-7,
      ROOT.kCyan-4,
      ROOT.kOrange+10,
      7,
      8,9,11,12,13,14,15  ]


bEdge_dijetSigmaMOverM=[0.0,0.025,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.18,0.2,0.25]
bEdge_SigOverE=[0.01,0.015,0.02,0.025,0.030,0.035,0.04,0.05,0.07,0.1,1.0]
drEdges=[0.0,0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6,4.0,4.4,4.8,5.2,5.6,6.0,7.0,8.0]
pTgg_overMggEdged=[0.05*i for i in range(22)]
pTg_overMgg=[0.0,0.25,0.33,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,4.0]
pTX_overMX=[0.0,0.33,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,4.0]
varList={
            'hgg_SigMOverM' : {} ,
            'h1leadJ_deepjetScore' : {} ,
            'h1subleadJ_deepjetScore' : {} ,
            'h2leadJ_deepjetScore' : {} ,
            'h2subleadJ_deepjetScore' : {} ,
            'h1leadG_customMVA' : {'rebin':{'binCount':4,'binEdges':None}} ,
            'h1subleadG_customMVA' : {'rebin':{'binCount':4,'binEdges':None}} ,
            'h1leadG_SigOverE' :{'xRange':(0.0,0.1) , 'rebin':{'binCount':len(bEdge_SigOverE)-1,'binEdges':bEdge_SigOverE}} ,
            'h1subleadG_SigOverE' : {'xRange':(0.0,0.1) , 'rebin':{'binCount':len(bEdge_SigOverE)-1,'binEdges':bEdge_SigOverE}} ,
            'H1bbCosThetaLeadJet' : {'rebin':{'binCount':4,'binEdges':None}} ,
            'H2bbCosThetaLeadJet' :  {'rebin':{'binCount':4,'binEdges':None}},
            'H1bbCosThetaLeadJet' : {'xRange':(0.0,1.0),'rebin':{'binCount':4,'binEdges':None}} ,
            'H2bbCosThetaLeadJet' : {'xRange':(0.0,1.0),'rebin':{'binCount':4,'binEdges':None}},
            'HggCosThetaLeadGamma' : {'xRange':(0.0,1.0),'rebin':{'binCount':4,'binEdges':None}} ,
            'H4bCosThetaLeadJet' : {'xRange':(0.0,1.0),'rebin':{'binCount':4,'binEdges':None}} ,
            'HHHCosThetaH1' : {'xRange':(0.0,1.0),'rebin':{'binCount':4,'binEdges':None}} ,
            'HHHCosThetaHgg' : {'xRange':(0.0,1.0),'rebin':{'binCount':4,'binEdges':None}} ,
            'phoJetMinDr' :{'xRange':(0.0,8.0) , 'rebin':{'binCount':len(drEdges)-1,'binEdges':drEdges}} ,
            'otherphoJetMinDr' :{'xRange':(0.0,8.0) , 'rebin':{'binCount':len(drEdges)-1,'binEdges':drEdges}} ,
            'pTh1jj_overMh1' :{'xRange':(0.0,1.0) } ,
            'pTgg_overMgg' :{'xRange':(0.0,1.0)   } ,
            'pTh2jj_overMh2' :{'xRange':(0.0,1.0) } ,
            'pTleadG_overMgg'   : {'xRange':(0.0,4.0),'rebin':{'binCount':len(pTg_overMgg)-1,'binEdges':pTg_overMgg} } ,
            'pTsubleadG_overMgg'   : {'xRange':(0.0,4.0),'rebin':{'binCount':len(pTg_overMgg)-1,'binEdges':pTg_overMgg} } ,
            'pTh1leadJ_overMh1'   : {'xRange':(0.0,4.0),'rebin':{'binCount':len(pTX_overMX)-1,'binEdges':pTX_overMX} } ,
            'pTh2leadJ_overMh2'   : {'xRange':(0.0,4.0),'rebin':{'binCount':len(pTX_overMX)-1,'binEdges':pTX_overMX} } ,
            'pTh1subleadJ_overMh1'   : {'xRange':(0.0,4.0),'rebin':{'binCount':len(pTX_overMX)-1,'binEdges':pTX_overMX} } ,
            'pTh2subleadJ_overMh2'   : {'xRange':(0.0,4.0),'rebin':{'binCount':len(pTX_overMX)-1,'binEdges':pTX_overMX} } ,
            'h1leadJ_deepjetScore' : {'yRange':(1.0,1e5)} ,
            'h1subleadJ_deepjetScore' : { 'yRange':(1.0,1e5) } ,
            'h2leadJ_deepjetScore' : {'yRange':(1.0,1e5)} ,
            'h2subleadJ_deepjetScore' : {'yRange':(1.0,1e5)} ,
            'rho' :{'xRange':(0.0,70.0) ,'rebin':{'binCount':4 , 'binEdges':None}   } ,
            'h1_dijetSigmaMOverM' : {'xRange':(0.0,0.25) , 'rebin':{'binCount':len(bEdge_dijetSigmaMOverM)-1,'binEdges':bEdge_dijetSigmaMOverM}} ,
            'h2_dijetSigmaMOverM' : {'xRange':(0.0,0.25) , 'rebin':{'binCount':len(bEdge_dijetSigmaMOverM)-1,'binEdges':bEdge_dijetSigmaMOverM}} ,
    }
#varList={
#        'h1subleadG_customMVA' : {'rebin':{'binCount':4,'binEdges':None}} ,
#        'h1leadG_SigOverE' :{'xRange':(0.0,0.1) , 'rebin':{'binCount':len(bEdge_SigOverE)-1,'binEdges':bEdge_SigOverE}} ,
#}
if __name__ == "__main__":
    candDict={}
    candDict['controlCands'] ='results/plots/analysis/controlRegion/dataMC_bdtVars'
    candDict['trippleHCands']='results/plots/analysis/signalRegion/dataMC_bdtVars'
#    candDict['validation_CR']='results/plots/analysis/validationCR'
#    candDict['validation_SR']='results/plots/analysis/validationSR'
    
    for preReweight in [False]:
        fListDict={
            'data2018'       : 'results/reweighted/data2018.root',
            'signal'         : 'results/reweighted/signal_hhhTo4b2g.root',
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
            'signal'         :'hhh#rightarrow 4b2#gamma'        ,
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
            print("Adding ",tag," : ",fListDict[tag])
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


        for var in varList:
            for cands in candDict:
                pltPrefix=candDict[cands]
                if not os.path.exists(pltPrefix):
                    print("Making folder :" ,pltPrefix)
                os.system("mkdir -p "+pltPrefix)    
            
                plots=[]
                cPars=getCanvasParams('GEN')
                rBin=1
                sc=1.0
                plot=getDefaultPlot(prefix=pltPrefix,name='XXX',cPars=cPars) 
                plot.legendPosition = (0.6,0.74,0.9,0.90)
                plot.descPosition   = (0.6,0.85)
                plot.desc =[]
                plot.yTitle = "#"  
                plot.xTitle = var
                plot.xRange = ('auto','auto') 
                plot.yRange = (0.9,7000.0) 
                plot.yRatioRange = (0.0,2.0) 
                if preReweight:
                    plot.yRatioRange = (0.0,15.0) 
                plots.append(plot)
                
                if 'yRange' in varList[var]:
                    plot.yRange=varList[var]['yRange']

                legend = ROOT.TLegend(plot.legendPosition[0],plot.legendPosition[1],plot.legendPosition[2],plot.legendPosition[3])
                legend.SetTextFont(42)
                legend.SetFillColor(0)

                                   
                dataHist=histStore['data2018'][cands]['flashggVars'][var].Clone()
                if 'rebin' in varList[var]:
                    dataHist=rebinTheHistogram( dataHist , varList[var]['rebin']['binCount'] ,varList[var]['rebin']['binEdges']  )
                legend.AddEntry(dataHist,"Data 2018","pe")
                print("tag  : data : ",dataHist.Integral()," bw : ",dataHist.GetBinWidth(1))
                dataHist.SetLineWidth(2)
                dataHist.SetLineColor(ROOT.kBlack)
                
                hDummy=None
                if 'xRange' in varList[var]:
                    hDummy=ROOT.TH1F("dummy","dummy",1,varList[var]['xRange'][0],varList[var]['xRange'][1])
                else:
                    hDummy=dataHist.Clone()
                plot.xRange=(hDummy.GetBinLowEdge(1) , hDummy.GetBinLowEdge(hDummy.GetNbinsX()+1))

                hs = ROOT.THStack("hs","Stacked 1D histograms");
                iCol=2

                histMC=None
                #for ky in ['ttX','ggJetsM80']:
                for ky in ['ttX','gJets','ggJetsM80']:
                    hist=histStore[ky][cands]['flashggVars'][var].Clone()
                    if 'rebin' in varList[var]:
                        hist=rebinTheHistogram( hist , varList[var]['rebin']['binCount'] ,varList[var]['rebin']['binEdges']  )
                    
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
                #h=ROOT.TH1F("xx",'xx',1,plot.xRange[0],plot.xRange[1])
                hDummy.Reset()
                hDummy.SetAxisRange(plot.yRange[0], plot.yRange[1], "Y")
                hDummy.SetXTitle( plot.xTitle )
                hDummy.SetYTitle(plot.yTitle)
                #h.SetYTitle( plot.yTitle )
                hDummy.GetYaxis().SetTitleOffset(1)
                hDummy.Draw()
                hDummy.Draw()
                pad2.cd()
                hR=hDummy.Clone()
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
                oFname=pltPrefix+'/'+var+'_ratio.png'
                if preReweight:
                    oFname=pltPrefix+'/'+var+'_preReweighting_ratio.png'
                canvas.SaveAs(oFname)


        for tag in fListDict:
            fileStore[tag].Close()
