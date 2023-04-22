import matplotlib.pyplot as plt
import mplhep as hep
import uproot as urt
import ROOT
import json,os
import numpy as np
import Util as utl
import scaleFactorUtil as sf
import argparse


hep.style.use("CMS")


if __name__=='__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-v',"--version", help="Version of the specific derivation ",default='vX')
    args = parser.parse_args()

    version = args.version
    prefixBase='/home/aravind/cernbox/work/trippleHiggs/hhhTo4b2gamma/genAnalysis/python/analysis/results/plots/analysis/jan23/scaleFactors'
    bkgsToSubtract=['gJet40ToInf','gJet40ToInf']
    bkgsToTake=['ggBox','ggBox1Bjet','ggBox2Bjet']
    if version=='v0':
        prefixBase+='_v0'    
        bkgsToTake=['ggBox','ggBox1Bjet','ggBox2Bjet']
        bkgsToSubtract=['gJet40ToInf','gJet40ToInf']
    elif version=='v1':
        prefixBase+='_v1'    
        bkgsToTake=['ggBox','ggBox1Bjet','ggBox2Bjet','gJet40ToInf']
        bkgsToSubtract=['gJet40ToInf']
    elif version=='v2':
        prefixBase+='_v2'    
        bkgsToTake=['ggBox','ggBox1Bjet','ggBox2Bjet','gJet40ToInf','gJet20To40']
        bkgsToSubtract=[]
    else:   
        print("version type  : ",version,"  not defined" )
        print(" Please enter a valid version")
        exit(1)


    
    fileDict={}
    with open('workarea/data/bdtNtuples/v5p0/filelist.json') as f:
        fileDict=json.load(f)
    yearsToProcess=['run2','2018','2016PreVFP','2016PostVFP','2017','2016']
    #yearsToProcess=['2018']
    yearsToValidate=['run2','2018','2016PreVFP','2016PostVFP','2017','2016']
    #yearsToValidate=['2018']
    varsToValidate=[]
    bkgToProcess=[ 'ggBox1Bjet',
                   'ggBox2Bjet', 
                   'ggBox',
                   'gJet19To40',
                   'gJet40ToInf']    
    
    binEdges=[0.0,20.0,30,40,60.,80,100.0,120.0,140.0,180.0,240.0,300.0,360.,420.,500.0,600.0,800.,1200.0]

    ggPt_Binnig=("","",240,0.0,1200)
    
    yearsToLoad=np.unique( yearsToProcess+yearsToValidate)

    categories={
        'inc': 'abs(leadingPhoton_pt) > -1 ',
    #    'BB' : 'abs(leadingPhoton_eta) < 1.44   && abs(subleadingPhoton_eta) < 1.44 ',
    #    'BE' : 'abs(leadingPhoton_eta) < 1.44   && abs(subleadingPhoton_eta) > 1.567 ',
    #    'EB' : 'abs(leadingPhoton_eta) > 1.567  && abs(subleadingPhoton_eta) < 1.44 ',
    #    'EE' : 'abs(leadingPhoton_eta) > 1.567  && abs(subleadingPhoton_eta) > 1.567 ',
    }
    # We read the tree from the file and create a RDataFrame, a class that
    # allows us to interact with the data contained in the tree.
    

    filesMade=[]
    rdataFrames={}
    for yr  in yearsToLoad:
        rdataFrames[yr]={'sig':{},'bkg':{},'data':{}}
    
        ky=list(fileDict[yr]['data'].keys())[0]
        fileName=fileDict[yr]['data'][ky]
        treeName = "trees/Data_13TeV_TrippleHTag_0"
        rdataFrames[yr]['data']['data']= ROOT.RDataFrame(treeName, fileName).Filter(' CMS_hgg_mass < 115.0  || CMS_hgg_mass >135.0')
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
            rdataFrames[yr]['bkg'][bkg]=ROOT.RDataFrame(treeName, fileName).Filter(' CMS_hgg_mass < 115.0  || CMS_hgg_mass >135.0')
            print("Registering datset : ",bkg," , ",yr," withh tree",treeName)
    
    
    ### Defining the bare Histograms : Total data and MC histograms that needed to be ratio
    if True:
        allSFactors={}
        for year  in yearsToProcess:
            allSFactors[year]={}
            saveBase=prefixBase+'/'+year+'/'
            os.system('mkdir -p '+saveBase)
            print("Processing Year : ",year)
            
            fileout=ROOT.TFile(saveBase+'/scaleFactors.root',"RECREATE")
            filesMade.append(fileout)
            for cat in categories:
                print("-"*10)
                print("Category : ",cat)
                filterStr=categories[cat]
                bkgHistList=[]
                bkgHistToSubsList=[]
                for ky in bkgsToTake:
                    h1=rdataFrames[year]['bkg'][ky].Filter(filterStr).Histo1D(ggPt_Binnig,"diphoton_pt","weight")
                    bkgHistList.append(h1)
                for ky in bkgsToSubtract:
                    h1=rdataFrames[year]['bkg'][ky].Filter(filterStr).Histo1D(ggPt_Binnig,"diphoton_pt","weight")
                    bkgHistToSubsList.append(h1)
                
                bkgHist=bkgHistList[0].Clone()
                for i in range(1,len(bkgHistList)):
                    bkgHist.Add(bkgHistList[i].Clone())
                bkgHistToSubstract=None
                if len(bkgHistToSubsList) > 0:
                    bkgHistToSubstract=bkgHistToSubsList[0].Clone()
                    for i in range(1,len(bkgHistToSubsList)):
                        bkgHistToSubstract.Add(bkgHistToSubsList[i].Clone())
                    bkgHistToSubstract=[bkgHistToSubstract]
                print("\tScaling to baground to luminosity : ",utl.lumiMap[year])
                bkgHist.Scale( float(utl.lumiMap[year] ) )
                
                dataHist=rdataFrames[year]['data']['data'].Filter(filterStr).Histo1D( ggPt_Binnig, "diphoton_pt", "weight")
        
                scaler=sf.getPtDependentScaleFactor(  name=year+'_'+cat,
                                              dataHist=dataHist,
                                              mcHists=[bkgHist],
                                              mcHistsToSubstract=bkgHistToSubstract,
                                              binEdges=binEdges,
                                              def_scaleFactor=0.0,
                                              cat=cat)
                print("       SCALE FACTOR val @ diphoton_pt = ", 0.0," ->",scaler.getSFForX(0.0,cat))
                print("       SCALE FACTOR val @ diphoton_pt = ",50.0," -> ",scaler.getSFForX(50.0,cat))
                print("       SCALE FACTOR val @ diphoton_pt = ",180.0," -> ",scaler.getSFForX(180.0,cat))
                print("       SCALE FACTOR val @ diphoton_pt = ",800.0," -> ",scaler.getSFForX(800.0,cat))
                print("       SCALE FACTOR val @ diphoton_pt = ",1800.0," -> ",scaler.getSFForX(1800.0,cat))

                allSFactors[year][cat]=scaler

                f1=plt.figure(figsize=(6,6))
                ax1=plt.subplot(1,1,1)
                hdata =scaler.getDataHist(cat)
                hmc   =scaler.getMCHist(cat)
                hep.histplot(hdata,ax=ax1,label='Data',color='b',yerr=False)
                hep.histplot(hmc,ax=ax1,label='MC',color='r',yerr=False)
                ax1.set_xlabel('pT$_{\gamma\gamma}$')
                ax1.set_ylabel('# Events')
                ax1.legend()
                ax1.grid(alpha=0.1,color='k')
                f1.savefig(saveBase+'/'+cat+'_dataMC.png',bbox_inches='tight')
                plt.close(f1)

                hsf =scaler.getScaleFactorHist(cat)
                f1=plt.figure(figsize=(6,6))
                ax1=plt.subplot(1,1,1)
                hep.histplot(hsf,ax=ax1,label=cat,color='b',yerr=False)
                ax1.set_xlabel('pT$_{\gamma\gamma}$')
                ax1.set_ylabel('# Events')
                ax1.legend()
                ax1.semilogy(0)
                ax1.grid(alpha=0.1,color='k')
                ax1.set_xlim([0.0,500])
                f1.savefig(saveBase+'/'+cat+'_scaleFactor.png',bbox_inches='tight')
                plt.close(f1)

            
                f2=plt.figure(figsize=(6,6))
                ax2=plt.subplot(1,1,1)
                hep.histplot(hsf,ax=ax2,label=cat,yerr=False)
                ax2.set_xlabel('pT$_{\gamma\gamma}$')
                ax2.set_ylabel('Scale Factor')
                ax2.legend()
                ax2.grid(alpha=0.1,color='k')
                f2.savefig(saveBase+'/'+cat+'_scalefactor.png',bbox_inches='tight')
                plt.close(f2)
                
                hsf.Write()
                hdata.Write()
                hmc.Write()
            
            dataForVal={}
            for yr  in yearsToValidate:
                for cat in allSFactors[year]:
                    scaler=allSFactors[year][cat]
                    filterStr=categories[cat]
                    valBase=saveBase+'/validation/'+yr+'/'+cat+'/'
                    os.system('mkdir -p '+valBase)

                    dataForVal[yr]={'data':{},'bkg':{}}
                    print("\t Validating for Year : ",yr, " [ from scalefactors from : ",year," ] ")
                    lumi=float(utl.lumiMap[yr])                
            
                    # loading scaled weights and variables  of interest
                    for ky in rdataFrames[yr]['bkg']:
                        dataForVal[yr]['bkg'][ky]=rdataFrames[yr]['bkg'][ky].Filter(filterStr).AsNumpy(["diphoton_pt","weight"]+varsToValidate)
                        print("\tProcessing bkg : ",ky)
                        if ky in bkgsToTake:
                            scaleF=[]
                            w2=[]
                            for i in range(len(dataForVal[yr]['bkg'][ky]["weight"])):
                                if i%25000==0:
                                    print("\t\t i = ",i," / " ,len(dataForVal[yr]['bkg'][ky]["weight"]) )
                                sFactor=scaler.getSFForX(dataForVal[yr]['bkg'][ky]["diphoton_pt"][i],cat)
                                scaleF.append(sFactor)
                                w2.append( dataForVal[yr]['bkg'][ky]["weight"][i]*sFactor )
                            dataForVal[yr]['bkg'][ky]["w2"]=np.array(w2)
                            dataForVal[yr]['bkg'][ky]["sf"]=np.array(scaleF)
                        else:
                            dataForVal[yr]['bkg'][ky]["w2"]=dataForVal[yr]['bkg'][ky]["weight"]
                            dataForVal[yr]['bkg'][ky]["sf"]=dataForVal[yr]['bkg'][ky]["weight"]*0+1.0
                    # loading data
                    dataForVal[yr]['data']= rdataFrames[yr]['data']['data'].Filter(filterStr).AsNumpy(["diphoton_pt","weight"])
                    dataForVal[yr]['data']["w2"]=dataForVal[yr]['data']["weight"]
                    dataForVal[yr]['data']["sf"]=dataForVal[yr]['data']["weight"]*0+1.0
                    
                    
                    # ploting weights after the scale factor is applied
                    for bkg in bkgsToTake:               
                        f3,ax=plt.subplots(1,1,figsize=(5,5))
                        hist_d=np.histogram( dataForVal[yr]['bkg'][bkg]["weight"] ) # default weight
                        hep.histplot(hist_d,ax=ax,color='r',label='pre Scaling')
                        hist_d=np.histogram( dataForVal[yr]['bkg'][bkg]["w2"] )
                        hep.histplot(hist_d,ax=ax,color='b',label='post Scaling')
                        ax.annotate(bkg,xy=(0.5,0.5),xycoords='axes fraction')
                        ax.set_xlabel("weight")
                        ax.grid(alpha=0.1,color='k')
                        f3.savefig(valBase+'/'+bkg+'_weights.png',bbox_inches='tight')
                        plt.close(f3)


                    ## Making DATA-MC Histograms
                    #hist_d=np.histogram( dataForVal[yr]['data']["diphoton_pt"],
                    #                     weights=dataForVal[yr]['data']["w2"],bins=binEdges)
                    hist_d=np.histogram( dataForVal[yr]['data']["diphoton_pt"],
                                         weights=dataForVal[yr]['data']["w2"],bins=np.linspace(0.0,400.0,20))
                    bins_=hist_d[1]
                    
                    allhists=[] ; label=[]
                    for bkg in utl.backgroundStackList:
                        if bkg not in dataForVal[yr]['bkg']:
                            continue
                        h=np.histogram( dataForVal[yr]['bkg'][bkg]["diphoton_pt"],
                                         weights=dataForVal[yr]['bkg'][bkg]["w2"]*lumi,
                                         bins=bins_ )
                        allhists.append(utl.getTH1FromNumpHist(h)) ; label.append(bkg)
                    
                    totHist=allhists[0].Clone()
                    for i in range(1,len(allhists)):
                        totHist.Add(allhists[i])

                    ## Ploting Data MC Histogram
                    f=plt.figure(figsize=(12,12))
                    ax=plt.subplot(4,1,(1,3))
                    hep.histplot(allhists,stack=True,label=label,histtype='fill',ax=ax)
                    hep.histplot(hist_d,ax=ax,color='k',label='Data')
                    ax.grid(alpha=0.3,color='k')
                    ax.semilogy()
                    ax.legend()
                    
                    ratioHist=utl.getTH1FromNumpHist(hist_d)
                    ratioHist.Divide(totHist)

                    axRatio=plt.subplot(4,1,4)
                    hep.histplot(ratioHist,ax=axRatio,color='m',label='MC')
                    axRatio.set_ylim(bottom=0.0)
                    axRatio.axhline(1.0,color='k')
                    axRatio.set_xlabel("diphoton_pt")
                    axRatio.grid(alpha=0.1,color='k')
                    f.savefig(valBase+'/'+"diphoton_pt"+'.png',bbox_inches='tight')
                    plt.close(f)
    for f in filesMade:
        f.Close()
