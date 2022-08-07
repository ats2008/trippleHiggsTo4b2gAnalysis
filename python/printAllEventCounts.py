from PlotUtil import *
from Util import *
import os

def getBinVsContent(hist):
    xAxis=hist.GetXaxis()
    nBins=hist.GetNbinsX()
    cutFlow={}
    
    for i in range(nBins):
        binName=xAxis.GetBinLabel(i)
        binContent=hist.GetBinContent(i)
        if binName=="":
            continue
        cutFlow[binName]=np.round(float(binContent),3)
    return cutFlow


if __name__ == "__main__":
    candDict={}
    candDict['allCands'] ='results/plots/analysis/allCands'
    candDict['controlCands'] ='results/plots/analysis/controlRegion'
    candDict['trippleHCands']='results/plots/analysis/signalRegion'
    candDict['validation_CR']='results/plots/analysis/validationCR'
    candDict['validation_SR']='results/plots/analysis/validationSR'
    
    for preReweight in [False]:
        fListDict={
            'data2018'       : 'results/reweighted/data2018.root',
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

        for cands in candDict:
            print("Candidate : ",cands)
            for tag in histStore:
                nHist=histStore[tag][cands]['nEvents'].Clone()
                counts=getBinVsContent(nHist)
                tw=0
                t=0
                if "totalWeighted" in counts:
                    tw=counts["totalWeighted"]
                if "total" in counts:
                    t=counts["total"]
                print("\t",tag," : ",tw,"(",t,")")
            print()

        for tag in fListDict:
            fileStore[tag].Close()
