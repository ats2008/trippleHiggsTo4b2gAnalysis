from PlotUtil import *
from Util import *
import os


if __name__ == "__main__":
    candDict={}
    candDict['trippleHCands']='results/plots/dataMC/'
    
    for preReweight in [False]:
        fListDict={
            'data2018'       : 'workarea/batch/fitterSkimsV8_noPeakingMVA/EGamma_alesauva-UL2018_0-10_6_4-v0-Run2018-12Nov2019_UL2018-v2.root',
            'signal'         : 'workarea/batch/fitterSkimsV8_noPeakingMVA/c3_1_c4_1_HHHto4b2gamma_UL17.root',
            'ggM80Inc'       : 'workarea/batch/fitterSkimsV8_noPeakingMVA/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa_v1.root',
            'ggM80Jbox1bjet' : 'workarea/batch/fitterSkimsV8_noPeakingMVA/output_DiPhotonJetsBox1BJet_MGG-80toInf_13TeV-sherpa.root',
            'ggM80Jbox2bjet' : 'workarea/batch/fitterSkimsV8_noPeakingMVA/output_DiPhotonJetsBox2BJets_MGG-80toInf_13TeV-sherpa.root',
            'ttgJets'        : 'workarea/batch/fitterSkimsV8_noPeakingMVA/output_TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root',
            'ttgg'           : 'workarea/batch/fitterSkimsV8_noPeakingMVA/output_TTGG_TuneCP5_13TeV-amcatnlo-pythia8.root',
            'tGJets'         : 'workarea/batch/fitterSkimsV8_noPeakingMVA/output_TGJets_TuneCP5_13TeV-amcatnlo-madspin-pythia8.root',
            'gJet20To40'     : 'workarea/batch/fitterSkimsV8_noPeakingMVA/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root',
            'gJet40ToInf'    : 'workarea/batch/fitterSkimsV8_noPeakingMVA/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root'
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
            histStore[tag]=getTheObjectsFromFile(fileStore[tag])['tagsDumper']['trees']

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
            rBin=2
            sc=58.90
            plot=getDefaultPlot(prefix=pltPrefix,name='XXX',cPars=cPars) 
            plot.legendPosition = (0.6,0.74,0.9,0.90)
            plot.descPosition   = (0.6,0.85)
            plot.desc =[]
            plot.yTitle = "# / BWidth GeV"  
            plot.xTitle = "m_{jj}^{1} [GeV]"
            plot.xRange = (90,160.0) 
            plot.yRange = (1e-3,2000.0) 
            plot.yRatioRange = (0.0,2.0) 
            if preReweight:
                plot.yRatioRange = (0.0,15.0) 
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
                hist=histStore[ky]['cat0_m1jj'].Clone()
                hist=hist.Rebin(rBin)
                hist.Scale(sc)
                
                if histMC==None:
                    histMC=hist.Clone()
                else:
                    histMC=histMC+hist
                
                legend.SetTextSize(0.025) 
                hist.SetFillColor(col[iCol]);       
                lege=legendDict[ky]
                legend.AddEntry(hist, lege, "f")

                iCol+=1
                print("tag  : ",ky," : ",hist.Integral()," ",hist.GetBinWidth(1))
                hs.Add(hist)
            h.SetYTitle(plot.yTitle.replace("BWidth",str(hist.GetBinWidth(1))))

###  ----------------------------------------------
            ky='signal'
            histSig=histStore[ky]['cat0_m1jj'].Clone()
            histSig=histSig.Rebin(rBin)
            histSig.Print()
            histSig.Scale(3e3*sc)
            histSig.Print()
            histSig.SetLineWidth(2)
            histSig.SetLineColor(ROOT.kRed);       
            
            for i in range(histSig.GetNbinsX()):
                print(i,histSig.GetBinContent(i))
            
            legend.SetTextSize(0.025) 
            lege='signal x3e3'
            legend.AddEntry(histSig, lege, "pe")
            iCol+=1
###  ----------------------------------------------

            
                
            dataHist=histStore['data2018']['cat0_m1jj'].Clone()
            dataHist=dataHist.RebinX(rBin)
            legend.AddEntry(dataHist,"Data 2018","pe")
            print("tag  : data : ",dataHist.Integral()," bw : ",dataHist.GetBinWidth(1))
            dataHist.SetLineWidth(2)
            
            pad1.cd()
            pad1.SetLogy()
            hs.Draw("same HIST")
            dataHist.Draw("same HIST pe")
            histSig.Draw("same HIST pe")
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
            oFname=pltPrefix+'/cat0_m1jj_ratio.png'
            if preReweight:
                oFname=pltPrefix+'/cat0_m1jj_preReweighting_ratio.png'
            canvas.SaveAs(oFname)


        for tag in fListDict:
            fileStore[tag].Close()
