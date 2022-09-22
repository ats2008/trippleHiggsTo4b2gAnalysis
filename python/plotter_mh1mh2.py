import * from PlotUtil

if __name__ == "__main__":

    plots=[]
    cPars=getCanvasParams('GEN')
    rBin=5
    plot=getDefaultPlot(prefix=prefix,name='XXX',cPars=cPars) 
    plot.legendPosition = (0.65,0.78,0.9,0.90)
    plot.descPosition   = (0.60,0.85)
    plot.desc =[]
    plot.yTitle = "a.u"  
    plot.xTitle = "mass [GeV]"
    plot.xRange = (0.0,200.0) 
    plot.yRange = (0.0,200.0) 
    plots.append(plot)
            
    canvas = ROOT.TCanvas("c", "", 700,700)
    canvas.SetGrid()
    plot.setPlotStyle()
    
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
    
    hist=histStore['c3_1_c4_1_FS']['hhhTo4b2gamma']['h1MassVsh2Mass'].Clone()
#     hist.RebinX(2)
#     hist.RebinY(2)
    hist.Draw('colz')
    x=hist.GetXaxis()
    x.SetRangeUser(40.0,210.0)
    y=hist.GetYaxis()
    y.SetRangeUser(40.0,210.0)
    hist.SetXTitle('m_{bb}(H1) [ GeV ]')
    hist.SetYTitle('m_{bb}(H2) [ GeV ]')
    
    signalR=ROOT.TEllipse(125.0,125.0,25.0,25.0)
    signalR.SetFillColorAlpha(0,0)
    signalR.SetLineColor(ROOT.kRed)
    signalR.Draw()
    
    for i in range(len(selbox)):
        selbox[i].Draw()
    CMSbox.Draw()
    extraTextBox.Draw()
    lumibox.Draw()
    
    canvas.SetRightMargin(0.13);
    canvas.SetLeftMargin(0.18);
    canvas.Draw()
    canvas.SaveAs(prefix+'/2b2bMassDistribution_dbgd.png')
