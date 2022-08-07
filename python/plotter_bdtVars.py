from PlotUtil import *
from Util import *
import os


if __name__ == "__main__":
    candDict={}
    candDict['controlCands'] ='results/plots/analysis/controlRegion/bdtVars'
    candDict['trippleHCands']='results/plots/analysis/signalRegion/bdtVars'
    candDict['validation_CR']='results/plots/analysis/validationCR/bdtVars'
    candDict['validation_SR']='results/plots/analysis/validationSR/bdtVars'
    varList=[
        'HHHCosThetaH1',
        'HHHCosThetaHgg',
        'H4bCosThetaLeadJet',
        'HggCosThetaLeadGamma',
        'H1bbCosThetaLeadJet',
        'H2bbCosThetaLeadJet',
        
        'pTgg_overMgg',
        'pTh1jj_overMh1',
        'pTh2jj_overMh2',
        
        'pTleadG_overMgg',
        'pTh1leadJ_overMh1',
        'pTh2leadJ_overMh2',
        
        'pTsubleadG_overMgg',
        'pTh1subleadJ_overMh1',
        'pTh2subleadJ_overMh2',
        
        'h1leadJ_deepjetScore',
        'h1subleadJ_deepjetScore',
        'h2leadJ_deepjetScore',
        'h2subleadJ_deepjetScore',
        
        'h1leadG_customMVA',
        'h1subleadG_customMVA',
        
        'h1leadG_SigOverE',
        'h1subleadG_SigOverE',
        'hgg_SigMOverM',
        'hgg_SigMOverMDecorrelated',
        'h1_dijetSigmaMOverM',
        'h2_dijetSigmaMOverM',
        
        'rho',
        'phoJetMinDr',
        'otherphoJetMinDr'
       ]

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
        legendDict={
            'data2018'       :'Data 2018'        ,
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



        for cands in candDict:
            pltPrefix=candDict[cands]
            if not os.path.exists(pltPrefix):
                print("Making folder :" ,pltPrefix)
            os.system("mkdir -p "+pltPrefix)    
        
############
            plots=[]
            cPars=getCanvasParams('Run2MC2017')
            for var in varList:
                rBin=1
                sc=1.0
                plot=getDefaultPlot(prefix=pltPrefix,name=var,cPars=cPars) 
                plot.legendPosition = (0.6,0.74,0.9,0.90)
                plot.descPosition   = (0.6,0.85)
                plot.desc =[]
                plot.yTitle = "# / 2 GeV"  
                plot.xTitle = "p_{T}^{#gamma#gamma} [GeV]"
                plot.xRange = (100,180.0) 
                plot.yRange = (1.0,10000.0) 
                plots.append(plot)
           
                for ky in ['siganl','gJets','ggJetsM80']:
                    hist=histStore[ky][cands]['hgg'][var].Clone()
                    hist.Rebin(rBin)
                    pparams=getDefaultPlotParams(col=i,marker=3)
                    pparams['Legend']=legendDict[ky]
                    pparams['MarkerColor']= col[i]; pparams['LineColor'] = col[i]
                    pparams['Options']='HIST '
                    aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
                    aplot.normalize = True ;aplot.scaleFactor = 1.0
                    plots[-1].addPlot(aplot)
############

            #canvas.SetLogy()
            #canvas.Draw()
            canvas = []
            for plot in plots:
                canvas.append(plot.plot())

        for tag in fListDict:
            fileStore[tag].Close()
