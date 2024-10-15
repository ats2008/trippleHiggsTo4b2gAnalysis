#!/usr/bin/env python3
import os
import copy
import ROOT
import PlotUtil as putil
import HistUtil as hutil
import CanvasSetup as cSetup

import dataDesc as ddc
import plotDesc as pdc
colourWheel=[ROOT.kBlack , ROOT.kGray , ROOT.kRed , ROOT.kGreen , ROOT.kBlue , ROOT.kYellow , ROOT.kMagenta , ROOT.kCyan , ROOT.kOrange , ROOT.kSpring , ROOT.kTeal , ROOT.kAzure , ROOT.kViolet , ROOT.kPink ]

#print(histDicts)

tag="Barrel"
fnames={ 
           'Barrel':"/home/aravind/cernbox/work/bs2mumug/photonID/scripts/BDTPlots/MVAScores.root",
           'EndCap' : "/home/aravind/cernbox/work/bs2mumug/photonID/scripts/BDTPlots/MVAScoresECap.root"
       }    

plots=[]

for tag in ['Barrel','EndCap']:
    fname=fnames[tag]    
    inputFile=ROOT.TFile.Open(fname)
    
    cPars=cSetup.getCanvasParams('GEN')
    plotDir="plots/MVAScores"
    if not os.path.exists(plotDir):
        os.system('mkdir -p '+plotDir)    
    
    plots=[]
    plots.append(putil.PlotCanvas(prefix=plotDir,canvasPars=cPars))
    hName='MVA Scores '+tag
    histDicts={}
    plots[-1].name   = hName  
    histDicts[hName]= pdc.getDefaultRatioDictionary(hName)
    plots[-1].xRange =histDicts[hName]['xRange']  
    plots[-1].yRange =histDicts[hName]['yRange']  
    plots[-1].yTitle =histDicts[hName]['yTitle']  
    plots[-1].xTitle =hName.replace(tag,"")
    plots[-1].desc   =histDicts[hName]['desc'] 
    plots[-1].legendPosition = histDicts[hName]['legendPosition'] 
    plots[-1].descPosition   = histDicts[hName]['descPosition']   
    plots[-1].logx = histDicts[hName]['logx']   
    plots[-1].logy = histDicts[hName]['logy']   
    histDicts[hName]['Options']['Normalize']=False 
    histName={}
    i=0
    for mva in ['BDTG']
        histName[mva]='aMVA_'+tag+'_BDTG/Method_BDTG/BDTG/MVA_BDTG_trainingRejBvsS'
        if i > len(colourWheel):
            i=0
        plotParams[mva]={'Legend':mva,'MarkerColor' : colourWheel[i],'LineColor':colourWheel[i],'LineStyle':9,'MarkerStyle':105,'DrawLegend':True ,'ScaleFactor':1.0}
        i+=1

    histDicts[hName]['legendPosition']=(0.35,0.60,0.73,0.75)
    histDicts[hName]['descPosition']   = (0.38,0.80)

    plots[-1].legendPosition = histDicts[hName]['legendPosition'] 
    plots[-1].descPosition   = histDicts[hName]['descPosition']   
    for tag in histName:
        histo0 = inputFile.Get(histName[tag])
        histo0.__class__ = ROOT.TH1F
        aplot = putil.Plot(Name=histName[tag],Histo=histo0,**plotParams[tag],**histDicts[hName]['Options'])
        plots[-1].addPlot(aplot)
    canvas=[]        
    for plot in plots:
        canvas.append(plot.plot())

    inputFile.Close()   
