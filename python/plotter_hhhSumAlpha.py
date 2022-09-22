from PlotUtil import *
from Util import *
import os

col=[ 0,1,2,4,6,ROOT.kOrange,ROOT.kMagenta+2,9,11,12,13,14,15  ]

if __name__ == "__main__":
    candDict={}
    candDict['allCands'] ='results/plots/analysis/allCands'
    varList=[
        'pT',
       ]
#    varList=[
#        'HHHCosThetaH1',
#        'HHHCosThetaHgg',
#        'H4bCosThetaLeadJet'
#        ]

    for preReweight in [False]:
        fListDict={
            'data2018'       : 'results/unityWeight/data2018.root',
            'signal'         : 'results/unityWeight/signal_hhhTo4b2g.root',
            'ggM80Inc'       : 'results/unityWeight/diphotonInclusive.root',
            'ggM80Jbox1bjet' : 'results/unityWeight/diphotonJetBox1bjet.root',
            'ggM80Jbox2bjet' : 'results/unityWeight/diphotonJetBox2bjet.root',
            'ttgJets'        : 'results/unityWeight/TTGJets.root',
            'ttgg'           : 'results/unityWeight/TTGG.root',
            'tGJets'         : 'results/unityWeight/TGJets.root',
            'gJet20To40'     : 'results/unityWeight/photonInclusive20To40.root',
            'gJet40ToInf'    : 'results/unityWeight/photonInclusive40ToInf.root'
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
                rBin=4
                sc=1.0
                plot=getDefaultPlot(prefix=pltPrefix,name=var,cPars=cPars) 
                plot.legendPosition = (0.6,0.8,0.9,0.95)
                plot.descPosition   = (0.6,0.85)
                plot.desc =[]
                plot.yTitle = "# / 2 GeV"  
                plot.xTitle = var
                plot.xRange = ('auto','auto') 
                plot.yRange = ('auto','auto') 
                plots.append(plot)
                
                i=1
                for ky in ['signal','ggJetsM80']:
                    hist=histStore[ky][cands]['hhh'][var].Clone()
                    hist.Rebin(rBin)
                    pparams=getDefaultPlotParams(col=i,marker=3)
                    pparams['Legend']=legendDict[ky]
                    pparams['MarkerColor']= col[i]; pparams['LineColor'] = col[i]
                    pparams['Options']='HIST e'
                    aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
                    aplot.normalize = True ;aplot.scaleFactor = 1.0
                    plots[-1].addPlot(aplot)
                    i+=1
############

            #canvas.SetLogy()
            #canvas.Draw()
            canvas = []
            for plot in plots:
                canvas.append(plot.plot())
                canvas[-1].Close()

  #      for tag in fListDict:
  #          fileStore[tag].Close()
