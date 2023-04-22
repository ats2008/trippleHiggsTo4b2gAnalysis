import  Util  as utl
import scaleFactorUtil as sfUtil

if __name__ == "__main__":
    fListDict={
        'data2018'       : 'results/forRewaiting/data2018.root',
        #'signal'         : 'results/forRewaiting/signal_hhhTo4b2g.root',
        'ggM80Inc'       : 'results/forRewaiting/diphotonInclusive.root',
        'ggM80Jbox1bjet' : 'results/forRewaiting/diphotonJetBox1bjet.root',
        'ggM80Jbox2bjet' : 'results/forRewaiting/diphotonJetBox2bjet.root',
        'ttgJets'        : 'results/forRewaiting/TTGJets.root',
        'ttgg'           : 'results/forRewaiting/TTGG.root',
        'tGJets'         : 'results/forRewaiting/TGJets.root',
        'gJet20To40'     : 'results/forRewaiting/photonInclusive20To40.root',
        'gJet40ToInf'    : 'results/forRewaiting/photonInclusive40ToInf.root'
    }
    
    f=ROOT.TFile("pTScaleFactorFile.root","RECREATE")
    
    allAddedHists= ROOT.TH1F("AllAddedHists"," ",1,0.0,1.0)
    allAddedHists.SetCanExtend(ROOT.TH1.kAllAxes);
    allSubsHists = ROOT.TH1F("AllSubstractedHists"," ",1,0.0,1.0)
    allSubsHists.SetCanExtend(ROOT.TH1.kAllAxes);

    histStore={}
    fileStore={}
    for tag in fListDict:
        print("Adding ",tag,"  : ",fListDict[tag])
        if tag in fileStore:
            fileStore[tag].Close()
        fileStore[tag]=ROOT.TFile(fListDict[tag],'READ')
        histStore[tag]=getTheObjectsFromFile(fileStore[tag])

    rebinC=1
    for histName in [ 'pT_BB'  ,'pT_BE'  , 'pT_EB' , 'pT_EE'  ]: 
        dataHist=histStore['data2018']['controlCands']['hgg'][histName].Clone()
        dataHist.Rebin(rebinC)
        mcHists=[]
        mcHistsToSubstract=[]

        for tag in ['ggM80Inc', 'ggM80Jbox1bjet','ggM80Jbox2bjet','gJet20To40','gJet40ToInf']:
            mcHists.append( histStore[tag]['controlCands']['hgg'][histName].Clone() )
            mcHists[-1].Rebin(rebinC)
            allAddedHists.Fill(tag,1)
        
        for tag in ['ttgJets', 'ttgg', 'tGJets']:
            mcHistsToSubstract.append( histStore[tag]['controlCands']['hgg'][histName].Clone() )
            mcHistsToSubstract[-1].Rebin( rebinC )
            allSubsHists.Fill(tag,1)

        binEdges=[0,25.0,50.0,75.0,100.0,150.0,200.0,300.0,400.0,800.0]

        pTDependentScaler=getPtDependentScaleFactor("pTDependent"+histName,dataHist,mcHists,mcHistsToSubstract,binEdges,def_scaleFactor=0.0)    
        
        d=ROOT.TCanvas()

        legend = ROOT.TLegend(0.65,0.78,0.8,0.88)
        legend.SetTextFont(42)
        legend.SetFillColor(0)
        f.cd()
        h=pTDependentScaler.mcHist.Clone()
        g=pTDependentScaler.dataHist.Clone()
        h.SetAxisRange(0.0,1e4, "Y")
        h.SetAxisRange(0.0,1500.0, "X")
        h.SetMarkerColor(8)
        h.SetLineColor(4)
        legend.AddEntry(h,"MC","pe")
        h.Draw()
        g.SetLineColor(2)
        legend.AddEntry(g,"Data","pe")
        g.Draw("same pe")
        legend.Draw()
        d.Draw()
        d.SaveAs(histName+"_dataMC.png")
        for cat in pTDependentScaler.scaleFactorHist:
            k=pTDependentScaler.scaleFactorHist[cat].Clone()
            k.SetAxisRange(0.0,1500.0, "X")
            k.Draw()
            d.Draw()
            d.SaveAs(histName+"_Scalefactor.png")

        pTDependentScaler.dataHist.Write()
        pTDependentScaler.mcHist.Write()
        allAddedHists.Write()
        allSubsHists.Write()
        for cat in pTDependentScaler.scaleFactorHist:
            pTDependentScaler.scaleFactorHist[cat].Write()
    f.Close()
