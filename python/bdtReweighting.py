import matplotlib.pyplot as plt
import mplhep as hep
import uproot as urt
import ROOT
import json,sys,os,argparse
import numpy as np
import Util as utl
import pickle

import hep_ml as hepml
from hep_ml.reweight import BinsReweighter, GBReweighter, FoldingReweighter
from hep_ml.preprocessing import IronTransformer


hep.style.use("CMS")
blind='CMS_hgg_mass < 115 || CMS_hgg_mass > 135.0'

def plotDataMCComparison(var,x_mc,x_data,w_data,w_mc,w_mc_binned,w_mc_bdt,bins=None,saveBase=None):

    ### Plotting the Data MC of different variables after reweighting for comparison study
    varID=0
    binEdges=bins
    # binEdges=np.arange(-4.0,4.0,20)
    dolog=True
    
    f=plt.figure(figsize=(24,8))
    
    if True : # base hist
        hd =np.histogram(x_data,weights=w_data,bins=binEdges)
        hmc=np.histogram(x_mc  ,weights=w_mc  ,bins=binEdges)
    
        hData=utl.getTH1FromNumpHist(hd)
        hMC  =utl.getTH1FromNumpHist(hmc)
        hRatio=hData.Clone()
        hRatio.Divide(hMC)
    
        ax1=plt.subplot(4,3,(1,7))
        hep.histplot(hData,ax=ax1,label='Data')
        hep.histplot(hMC,ax=ax1,label='MC')
        #     ax1.set_xticks([])
        ax1.grid()
    
        ax2=plt.subplot(4,3,10)
        hep.histplot(hRatio,ax=ax2,label='MC')
        ax2.axhline(1.0,c='k')
    #     ax2.set_xlabel([0.0,2.0])
        ax2.grid(color='red',which='minor')
        ax1.annotate('Raw  Weights' ,(0.60,0.75),xycoords='axes fraction', fontsize=15)
        ax1.annotate(var,(0.60,0.65),xycoords='axes fraction', fontsize=15)
        ax1.legend()
        if dolog:
            ax1.semilogy()
    
    if True : # reweighted hist with binned reweighting
        hd =np.histogram(x_data,weights=w_data     ,bins=binEdges )
        hmc=np.histogram(x_mc  ,weights=w_mc_binned,bins=binEdges )
    
        hData=utl.getTH1FromNumpHist(hd)
        hMC=utl.getTH1FromNumpHist(hmc)
        hRatio=hData.Clone()
        hRatio.Divide(hMC)
    
        ax3=plt.subplot(4,3,(2,8))
        hep.histplot(hData,ax=ax3,label='Data')
        hep.histplot(hMC,ax=ax3,label='MC')
        #     aax3.set_xticks([])
        ax3.grid()
    
        ax4=plt.subplot(4,3,11)
        hep.histplot(hRatio,ax=ax4,label='MC')
        ax4.axhline(1.0,c='k')
        ax4.set_ylim([0.0,4.0])
        ax4.grid(color='red',which='minor')
        ax3.annotate(var,(0.60,0.65),xycoords='axes fraction',fontsize=15)
        ax3.annotate('Binned Reweighting',(0.60,0.75),xycoords='axes fraction',fontsize=15)
        ax3.legend()
        if dolog:
            ax3.semilogy()
    
    if True : # reweighted hist with bdt reweighting
        hd =np.histogram(x_data,weights=w_data,bins=binEdges)
        hmc=np.histogram(x_mc  ,weights=w_mc_bdt,bins=binEdges)
    
        hData=utl.getTH1FromNumpHist(hd)
        hMC=utl.getTH1FromNumpHist(hmc)
        hRatio=hData.Clone()
        hRatio.Divide(hMC)
    
        ax3=plt.subplot(4,3,(3,9))
        hep.histplot(hData,ax=ax3,label='Data')
        hep.histplot(hMC,ax=ax3,label='MC')
        #     aax3.set_xticks([])
        ax3.grid()
    
        ax4=plt.subplot(4,3,12)
        hep.histplot(hRatio,ax=ax4,label='MC')
        ax4.axhline(1.0,c='k')
        ax4.set_ylim([0.0,4.0])
        ax4.grid(color='red',which='minor')
        ax3.annotate(var,(0.60,0.65),xycoords='axes fraction' , fontsize=15)
        ax3.annotate('BDT Reweighting',(0.60,0.75),xycoords='axes fraction',fontsize=15)
        ax3.legend()
        if dolog:
            ax3.semilogy()
        if saveBase:
            foutname=saveBase+'/'+var+'.jpeg'
            print("Saving file : ",foutname)
            f.savefig(foutname,bbox_inches='tight')
        plt.close(f)

if __name__=='__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-v',"--version", help="Version of the specific derivation ",default='')
    parser.add_argument('-t',"--transformData", help="trnsoform data using IronTransformer ",default=False,action='store_true')
    args = parser.parse_args()
    version = args.version
    doIronTrnsformation=False
    
    fileDict={}
    with open('../workarea/data/analysisNtuples/v2p0/filelist.json') as f:
        fileDict=json.load(f)
    
    prefixBase='/home/aravind/cernbox/work/trippleHiggs/hhhTo4b2gamma/genAnalysis/python/analysis/results/plots/analysis/jan27/bdtReweighting'
    
    varToPlot=[
         'diphoton_pt' ,
         "nonResonantMVA_v0",
         "quadjet_0_deepJetScore" ,
         "quadjet_1_deepJetScore" ,
         "quadjet_2_deepJetScore" ,
         "quadjet_3_deepJetScore" ,
         "sumScore_4j" ,
         "sumScore_3j" ,
         'h1bb_mass' ,
         'h2bb_mass' ,
         "h1bb_eta" ,
         "h1bb_phi" ,
         "h1bb_pt"  ,
         "h2bb_eta" ,
         "h2bb_phi" ,
         "h2bb_pt"  ,
         'diphoton_pt',
         'diphoton_eta',
         'CMS_hgg_mass'
    ]

    ## Making the variables for th reweighter training
    varForReweighting=[
        'leadingPhoton_pt','leadingPhoton_eta',
        'pT_h1leadJ' ,'eta_h1leadJ',
        'pT_h1subleadJ' ,'eta_h1subleadJ',
        'pT_h2leadJ' ,'eta_h2leadJ',
        'pT_h2subleadJ' ,'eta_h2subleadJ'
    ] 
    if args.transformData:
        prefixBase+='_ironTransformed/'    
        doIronTrnsformation=True
    if 'WithScores'==version:
        prefixBase+='_'+version+'/'    
        varForReweighting+=[
                "quadjet_0_deepJetScore" ,
                "quadjet_1_deepJetScore" ,
                "quadjet_2_deepJetScore" ,
                "quadjet_3_deepJetScore" 
              ]
    elif 'WithScoresAndSums'==version:
        prefixBase+='_'+version+'/'    
        varForReweighting+=[
                "quadjet_0_deepJetScore" ,
                "quadjet_1_deepJetScore" ,
                "quadjet_2_deepJetScore" ,
                "quadjet_3_deepJetScore" ,
                "sumScore_4j" ,
                "sumScore_3j" 
              ]
    
    else:
        prefixBase+=version+'/'    


    yearsToProcess=['2018','2017','2016PreVFP','2016PostVFP','run2']
    bkgToProcess=[ 'ggBox1Bjet',
                   'ggBox2Bjet', 
                   'ggBox', 
                   'gJet20To40',
                   'gJet40ToInf'
                   ]
    varToBinMap={}
    with open('data/binMap.json','r') as f:
        varToBinMap=json.load(f)['varToBinMap']
        

    # We read the tree from the file and create a RDataFrame, a class that
    # allows us to interact with the data contained in the tree.
    rdataFrames={}
    for yr  in yearsToProcess:
        rdataFrames[yr]={'sig':{},'bkg':{},'data':{}}
        try:
            fileName=fileDict[yr]['sig']['ggHHH']
        except:
            print(yr)
        treeName = "trees/ggHHH_125_13TeV"
        rdataFrames[yr]['sig']['ggHHH']= ROOT.RDataFrame(treeName, fileName)
        print("Registering datset : ggHHH , ",yr," withh tree",treeName)
        
        ky=list(fileDict[yr]['data'].keys())[0]
        fileName=fileDict[yr]['data'][ky]
        treeName = "trees/Data_13TeV_TrippleHTag_0"
        rdataFrames[yr]['data']['data']= ROOT.RDataFrame(treeName, fileName).Filter(blind)
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
            rdataFrames[yr]['bkg'][bkg]=ROOT.RDataFrame(treeName, fileName).Filter(blind)
            print("Registering datset : ",bkg," , ",yr," withh tree",treeName)
                       
    
    # reading the vars from file    
    oVars=['weight','weight_v0','lumi']
    varsToGet=[ i for i in np.unique(varForReweighting + oVars + varToPlot)]
    allVars = np.unique(varForReweighting+varToPlot)
    
    dataStore={}
    for yr in yearsToProcess:
        models={}
        print("Processing for year : ",yr)
        saveBase=prefixBase+'/'+yr+'/'
        os.system('mkdir -p '+saveBase)
        dataStore[yr]={}
        bkgStack=None
        dataStore[yr]['bkg']={}
        for ky in rdataFrames[yr]['bkg']:
            dataStore[yr]['bkg'][ky]=rdataFrames[yr]['bkg'][ky].AsNumpy(varsToGet)
        dataStore[yr]['data']={}
        dataStore[yr]['data']['data']=rdataFrames[yr]['data']['data'].AsNumpy(varsToGet)
        
        data={
                'data':{'x':None,'weight':None},
                'mc'  :{'x':None,'weight':None}
             }
        
        isFirst=True
        print("\tLoading bkg ! ")
        for ky in dataStore[yr]['bkg']:
            print("\t\t adding : ",ky)
            xMC=np.stack([dataStore[yr]['bkg'][ky][k] for k in varForReweighting])
            if isFirst:
                data['mc']['x']=xMC
                data['mc']['weight']      =dataStore[yr]['bkg'][ky]['weight_v0']*dataStore[yr]['bkg'][ky]['lumi']
                data['mc']['weight_bined']=dataStore[yr]['bkg'][ky]['weight']
                isFirst=False
            else:
                data['mc']['x']=np.concatenate( [data['mc']['x'],xMC ] , axis=-1)
                w=dataStore[yr]['bkg'][ky]['weight_v0']*dataStore[yr]['bkg'][ky]['lumi']
                data['mc']['weight']=np.concatenate( [data['mc']['weight'],
                                                      w ] ,
                                                      axis=-1)
                data['mc']['weight_bined']=np.concatenate( [data['mc']['weight_bined'],
                                                      dataStore[yr]['bkg'][ky]['weight'] ] ,
                                                      axis=-1)
        
        print("\tLoading data ! ")
        isFirst=True
        for ky in dataStore[yr]['data']:
            xData=np.stack([dataStore[yr]['data'][ky][k] for k in varForReweighting])
            if isFirst:
                data['data']['x']=xData
                data['data']['weight']=dataStore[yr]['data'][ky]['weight']
                isFirst=False
            else:
                data['data']['x']=np.concatenate( [data['data']['x'],xData ] , axis=-1)
                data['data']['weight']=np.concatenate( [data['data']['weight'],
                                                        dataStore[yr]['data'][ky]['weight_v0'] ] ,
                                                      axis=-1)
        transformedData={
                'data':{},
                'mc'  :{}
             }
        
        if doIronTrnsformation: 
            print("\tDoing the Iron transormation on data")
            f=plt.figure()
            transformer = IronTransformer().fit(data['data']['x'].T,sample_weight=data['data']['weight'])
            model['dataTransformer'] = transformer
            transformedData['data']['x']     = transformer.transform(data['data']['x'].T).to_numpy().T
            transformedData['data']['weight']= data['data']['weight']
            transformedData['mc']['x']     = transformer.transform(data['mc']['x'].T).to_numpy().T
            transformedData['mc']['weight']= data['mc']['weight']
            transformedData['mc']['weight_bined']= data['mc']['weight_bined']
            
            
            _=plt.hist(transformedData['data']['x'][0],weights=transformedData['data']['weight'],
                       bins=40,histtype='step',label='data')
            _=plt.hist(transformedData['mc']['x'][0]  ,weights=transformedData['mc']['weight'],
                       bins=40,histtype='step',label='mc')
            plt.legend()
            plt.semilogy()
            plt.close(f)
        #  Define and Train the reweighter  
        dataToUse=data #transformedData
        if doIronTrnsformation:
            dataToUse=transformedData

        reweighter_base = GBReweighter(n_estimators=40,min_samples_leaf=300,max_depth=3, gb_args={'subsample': 0.5})
        reweighter = FoldingReweighter(reweighter_base, n_folds=3)
        
        print("\tFitting the reweighter ! ")
        reweighter.fit(original=dataToUse['mc']['x'].T,original_weight=dataToUse['mc']['weight'],
                       target=dataToUse['data']['x'].T,
                       target_weight=dataToUse['data']['weight'])
        model['reweighter']=reweighter
        print("\t Obtainning then Scale factors ! ")
        MC_weights = reweighter.predict_weights(dataToUse['mc']['x'].T)
        #MC_weights = dataToUse['mc']['weight']
        f=plt.figure()    
        plt.hist(dataToUse['mc']['weight'],histtype='step',label='raw',bins=80)
        plt.hist(dataToUse['mc']['weight_bined'],histtype='step',label='binScaled',bins=80)
        plt.hist(MC_weights,histtype='step',label='bdtSclaed',bins=80)
        plt.legend()
        plt.semilogy()
        plt.semilogx()
        plt.savefig(saveBase+'/'+yr+'_weightDistributions',bbox_inches='tight')
        plt.close(f)
    
        w_data      = data['data']['weight']
        w_mc        = data['mc']['weight'] 
        w_mc_binned = data['mc']['weight_bined'] 
        w_mc_bdt    = MC_weights 
        for var in allVars:
            x_mc  =np.concatenate([dataStore[yr]['bkg'][ky][var] for ky in bkgToProcess])
            x_data=np.concatenate([dataStore[yr]['data'][ky][var] for ky in dataStore[yr]['data']])
            binEdges=np.linspace(varToBinMap[var][3] , varToBinMap[var][4] ,varToBinMap[var][2] )
            plotDataMCComparison(var,
                                 x_mc,
                                 x_data,
                                 w_data,w_mc,w_mc_binned,w_mc_bdt,
                                 bins=binEdges,
                                 saveBase=saveBase)
        foutname=saveBase+'/'+'savedModel.pkl'
        with open(foutname,'w') as f:
            pickle.dump(model,out)
            

