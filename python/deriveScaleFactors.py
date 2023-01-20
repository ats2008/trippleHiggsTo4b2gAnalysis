import matplotlib.pyplot as plt
import mplhep as hep
import uproot3 as urt
import ROOT
import json,os
import numpy as np
import Util as utl
import scaleFactorUtil as sf


hep.style.use("CMS")


if __name__=='__main__':
    
    prefixBase='/home/aravind/cernbox/work/trippleHiggs/hhhTo4b2gamma/genAnalysis/python/analysis/results/plots/analysis/jan20/scaleFactors'
    
    fileDict={}
    with open('../workarea/data/bdtNtuples/v5p0/filelist.json') as f:
        fileDict=json.load(f)
    
    yearsToProcess=['run2','2018','2016PreVFP','2016PostVFP','2017']
    bkgToProcess=[ 'ggBox1Bjet',
                   'ggBox2Bjet', 
                   'ggBox',
                   'gJet20To40',
                   'gJet40ToInf']    
    binEdges=[0.0,50.,100.0,200.0,300.0,500.0,800.,1200.0]
    bkgsToTake=['ggBox','ggBox1Bjet','ggBox2Bjet']
    ggPt_Binnig=("","",240,0.0,1200)
    categories={
        'inc': 'abs(leadingPhoton_eta) > -1 ',
        'BB' : 'abs(leadingPhoton_eta) < 1.44   && abs(subleadingPhoton_eta) < 1.44 ',
        'BE' : 'abs(leadingPhoton_eta) < 1.44   && abs(subleadingPhoton_eta) > 1.567 ',
        'EB' : 'abs(leadingPhoton_eta) > 1.567  && abs(subleadingPhoton_eta) < 1.44 ',
        'EE' : 'abs(leadingPhoton_eta) > 1.567  && abs(subleadingPhoton_eta) > 1.567 ',
    }
    # We read the tree from the file and create a RDataFrame, a class that
    # allows us to interact with the data contained in the tree.
    rdataFrames={}
    for yr  in yearsToProcess:
        rdataFrames[yr]={'sig':{},'bkg':{},'data':{}}
    
        fileName=fileDict[yr]['sig']['ggHHH']
        treeName = "trees/ggHHH_125_13TeV"
        rdataFrames[yr]['sig']['ggHHH']= ROOT.RDataFrame(treeName, fileName)
        print("Registering datset : ggHHH , ",yr," withh tree",treeName)
    
        ky=list(fileDict[yr]['data'].keys())[0]
        fileName=fileDict[yr]['data'][ky]
        treeName = "trees/Data_13TeV_TrippleHTag_0"
        rdataFrames[yr]['data']['data']= ROOT.RDataFrame(treeName, fileName)
        print("Registering datset : data , ",yr," withh tree",treeName)
    
        rdataFrames[yr]['bkg']={}
        treeName = "trees/bkg_13TeV_TrippleHTag_0"
        for bkg in bkgToProcess:
            if bkg not in fileDict[yr]['bkg']:
                print()
                print("FILE NOT FOUNG FOR : ",yr," background ",bkg)
                print()
                continue
            fileName = fileDict[yr]['bkg'][bkg]
            rdataFrames[yr]['bkg'][bkg]=ROOT.RDataFrame(treeName, fileName)
            print("Registering datset : ",bkg," , ",yr," withh tree",treeName)
    
    
    ### Defining the bare Histograms
    if True:
        allHistsForSF={}
        for yr  in yearsToProcess:
            saveBase=prefixBase+'/'+yr+'/'
            os.system('mkdir -p '+saveBase)
            allHistsForSF[yr]={'data':{},'bkg':{}}
            print("Processing Year : ",yr)
            f=ROOT.TFile(saveBase+'/scaleFactors.root',"RECREATE")
            
            f2=plt.figure(figsize=(6,6))
            ax2=plt.subplot(1,1,1)
            for cat in categories:
                print("-"*10)
                print("Category : ",cat)
                filterStr=categories[cat]
                bkgHistList=[]
                for ky in bkgsToTake:
                    h1=rdataFrames[yr]['bkg'][ky].Filter(filterStr).Histo1D(ggPt_Binnig,"diphoton_pt","weight")
                    bkgHistList.append(h1)
                
                bkgHist=bkgHistList[0].Clone()
                bkgHist.Scale( float(utl.lumiMap[yr] ) )

                for i in range(1,len(bkgHistList)):
                    bkgHist.Add(bkgHistList[i].Clone())
                dataHist=rdataFrames[yr]['data']['data'].Filter(filterStr).Histo1D( ggPt_Binnig, "diphoton_pt", "weight")
        
                scaler=sf.getPtDependentScaleFactor(  name=yr+'_'+cat,
                                              dataHist=dataHist,
                                              mcHists=[bkgHist],
                                              mcHistsToSubstract=None,
                                              binEdges=binEdges,
                                              def_scaleFactor=0.0)
                f1=plt.figure(figsize=(6,6))
                ax1=plt.subplot(1,1,1)
                hdata =scaler.getDataHist()
                hmc   =scaler.getMCHist()
                hep.histplot(hdata,ax=ax1,label='Data',color='b',yerr=False)
                hep.histplot(hmc,ax=ax1,label='MC',color='r',yerr=False)
                ax1.set_xlabel('pT$_{\gamma\gamma}$')
                ax1.set_ylabel('# Events')
                ax1.legend()
                ax1.grid(alpha=0.1,color='k')
                f1.savefig(saveBase+'/'+cat+'_dataMC.png',bbox_inches='tight')

                hsf =scaler.getScaleFactorHist()
                
                f1=plt.figure(figsize=(6,6))
                ax1=plt.subplot(1,1,1)
                hep.histplot(hsf,ax=ax1,label=cat,color='b',yerr=False)
                ax1.set_xlabel('pT$_{\gamma\gamma}$')
                ax1.set_ylabel('# Events')
                ax1.legend()
                ax1.grid(alpha=0.1,color='k')
                f1.savefig(saveBase+'/'+cat+'_scaleFactor.png',bbox_inches='tight')

                hep.histplot(hsf,ax=ax2,label=cat,yerr=False)
                
                hsf.Write()
                hdata.Write()
                hmc.Write()

            ax2.set_xlabel('pT$_{\gamma\gamma}$')
            ax2.set_ylabel('Scale Factor')
            ax2.legend()
            ax2.grid(alpha=0.1,color='k')
            ax2.semilogy()
            f2.savefig(saveBase+'/scalefactor.png',bbox_inches='tight')
            f.Close()
                                   
