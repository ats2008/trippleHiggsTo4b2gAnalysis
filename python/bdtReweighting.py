import matplotlib.pyplot as plt
import mplhep as hep
import uproot as urt
import ROOT
import json,sys,os,argparse
import numpy as np
import Util as utl
import pickle
import prettytable as ptab

import scaleFactorUtil as scl
import statUtil as stat
import hep_ml as hepml
import prettytable


hep.style.use("CMS")
blind='CMS_hgg_mass < 115 || CMS_hgg_mass > 135.0'

def main():
 
    parser = argparse.ArgumentParser()
    parser.add_argument('-v',"--version", help="Version of the specific derivation ",default='')
    parser.add_argument('-i',"--inputFile", help="Input File",default=None)
    parser.add_argument('-y',"--year", help="Year",default='2018')
    parser.add_argument('-t',"--transformData", help="trnsoform data using IronTransformer ",default=False,action='store_true')
    args = parser.parse_args()
    version = args.version
    doIronTrnsformation=False
    inputFile  = args.inputFile
    saveOutput = True
    outDict={}

    if inputFile:
        saveOutput=False
    fileDict={}
    with open('../workarea/data/analysisNtuples/v2p0/filelist.json') as f:
        fileDict=json.load(f)
    
    prefixBase='/home/aravind/cernbox/work/trippleHiggs/hhhTo4b2gamma/genAnalysis/python/analysis/results/plots/analysis/feb1/bdtReweighting'
    
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
                "quadjet_0_deepJetScore" ,"quadjet_1_deepJetScore" ,"quadjet_2_deepJetScore" ,"quadjet_3_deepJetScore" 
              ]
    elif 'WithScoresAndSums'==version:
        prefixBase+='_'+version+'/'    
        varForReweighting+=[
                "quadjet_0_deepJetScore" ,"quadjet_1_deepJetScore" ,"quadjet_2_deepJetScore" ,
                "quadjet_3_deepJetScore" ,"sumScore_4j" , "sumScore_3j" 
              ]
    elif 'AllVars'==version:
        prefixBase+='_'+version+'/'    
        varForReweighting+=[
                "CMS_hgg_mass",  "H1H2JetDrMax",  "H1bbCosTheta_hhhF",  "HH4bCosTheta_hhhF", "HggCosTheta_hhhF",
                "LeadJetAbsCosThetaMax","LeadJetDrMaxWithOtherJets","LeadJetDrMinWithOtherJets","PhoJetMaxDr","PhoJetMaxDrOther",
                "PhoJetMinDr","PhoJetMinDrOther","absCosThetaH4bHgg","customLeadingPhotonIDMVA","customSubLeadingPhotonIDMVA",
                "diphoton_eta","diphoton_pt","eta_h1leadJ","eta_h1subleadJ","eta_h2leadJ",
                "eta_h2subleadJ","h1_dijetSigmaMOverM","h1bb_eta","h1bb_mass","h1bb_phi",
                "h1bb_pt","h2_dijetSigmaMOverM","h2bb_eta","h2bb_mass","h2bb_phi",
                "h2bb_pt","leadingPhoton_eta","leadingPhoton_pt","pT_h1leadJ","pT_h1subleadJ",
                "pT_h2leadJ","pT_h2subleadJ","pTh1leadJ_overMh1","pTh1subleadJ_overMh1","pTh2leadJ_overMh2",
                "pTh2subleadJ_overMh2","pThgg_overMgg","pTleadG_overMgg","pTsubleadG_overMgg","quadjet_0_deepJetScore",
                "quadjet_1_deepJetScore","quadjet_2_deepJetScore","quadjet_3_deepJetScore","r_HH","scalarPtSum4b",
                "scalarPtSum4b2g","scalarPtSumHHH","sumScore_3j","sumScore_4j","trihiggs_mass",
                "trihiggs_pt","ttH_MET" 
              ]
    
    else:
        prefixBase+=version+'/'    
    varForReweighting_=np.unique(varForReweighting)
    varForReweighting=[ i for i in varForReweighting_]
    outDict['varForReweighting']=varForReweighting
    yearsToProcess_all=['2018','2017','2016PreVFP','2016PostVFP','run2']
    
    yearsToProcess=[]
    for yr in args.year.split(","):
        if yr not in yearsToProcess_all:
            print(yr," not in catalogue. Skipping !! ")
            continue
        yearsToProcess.append(yr)

    bkgToProcess=[ 'ggBox1Bjet',
                   'ggBox2Bjet', 
                   'ggBox', 
                   'gJet20To40',
                   'gJet40ToInf'
                   ]
    bkgToFit = [
            'ggBox',
            'ggBox1Bjet',
            'ggBox2Bjet'
    ]
    bkgToProcess = [
            'ggBox',
            'ggBox1Bjet',
            'ggBox2Bjet'
    ]
    varToBinMap={}
    with open('data/binMap.json','r') as f:
        varToBinMap=json.load(f)['varToBinMap']
        
    print("Variables being used for training : ",end="\n\t")
    for i in range(len(varForReweighting)):
        print(varForReweighting[i],end = " , ")
        if i%4==0:
            print(end="\n\t")
    print()
    print("Outputs are stored in ",prefixBase)
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
        print("Registering datset : ggHHH , ",yr," with tree",treeName)
        
        ky=list(fileDict[yr]['data'].keys())[0]
        fileName=fileDict[yr]['data'][ky]
        treeName = "trees/Data_13TeV_TrippleHTag_0"
        rdataFrames[yr]['data']['data']= ROOT.RDataFrame(treeName, fileName).Filter(blind)
        print("Registering datset : data , ",yr," with tree",treeName)
        
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
            print("Registering datset : ",bkg," , ",yr," with tree",treeName)
                       
    
    # reading the vars from file    
    oVars=['weight','weight_v0','lumi']
    varsToGet=[ i for i in np.unique(varForReweighting + oVars + varToPlot)]
    allVars = np.unique(varForReweighting+varToPlot)
    
    for var in allVars:
        if var not in varToBinMap:
            print(var,"  Not in the bin Map")

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
                'data'    :{'x':None,'weight':None},
                'mc'      :{'x':None,'weight':None},
                'mc_fit'  :{'x':None,'weight':None}
             }
        
        isFirst=True
        print("\tLoading bkg ! ")
        for ky in dataStore[yr]['bkg']:
            if ky not in bkgToFit:
                continue
            print("\t\t adding background ",ky)
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
        print()                                                      
        for ky in dataStore[yr]['bkg']:
            if ky not in bkgToFit:
                continue
            print("\t\t adding background for fitting ",ky)
            xMC=np.stack([dataStore[yr]['bkg'][ky][k] for k in varForReweighting])
            if isFirst:
                data['mc_fit']['x']=xMC
                data['mc_fit']['weight']      =dataStore[yr]['bkg'][ky]['weight_v0']*dataStore[yr]['bkg'][ky]['lumi']
                data['mc_fit']['weight_bined']=dataStore[yr]['bkg'][ky]['weight']
                isFirst=False
            else:
                data['mc_fit']['x']=np.concatenate( [data['mc']['x'],xMC ] , axis=-1)
                w=dataStore[yr]['bkg'][ky]['weight_v0']*dataStore[yr]['bkg'][ky]['lumi']
                data['mc_fit']['weight']=np.concatenate( [data['mc']['weight'],
                                                      w ] ,
                                                      axis=-1)
                data['mc_fit']['weight_bined']=np.concatenate( [data['mc']['weight_bined'],
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
        
        bdtReweighter=None
        data_x = data['data']['x']
        mc_x   = data['mc']['x'] 
        data_w = data['data']['weight']
        mc_w   = data['mc']['weight'] 

        print("\t  nData : ",data_x.shape[1])
        print("\t  nMC   : ",mc_x.shape[1])
        outDict['nData'] = str( data_x.shape[1] )
        outDict['nMC'] = str( mc_x.shape[1] )

        if not inputFile:
            bdtReweighter = scl.bdtScaler( balanceDataMCCounts=True)
            bdtReweighter.setModel( n_estimators=40,min_samples_leaf=100,max_depth=2 ) #,gb_args={'subsample': 0.5}
            if doIronTrnsformation: 
                print("\tSetting the Iron transormation ")
                bdtReweighter.setTransformation(True)

            print( "\tFitting the reweighter ! " )
            bdtReweighter.fit( data['data']['x']  , data['data']['weight']  , data['mc']['x']  , data['mc']['weight'] )
            print("\tDone ! ")
            data_x,data_w = bdtReweighter.getData()
            mc_x,mc_w     = bdtReweighter.getMC()
            print("setting the mc_x with ",mc_x.shape)
            bdtReweighter.setClassifier(n_estimators=10, learning_rate=0.8 , max_depth=3)
            print("\tFitting the classifier ! ")
            bdtReweighter.trainClassifierModel()
            print("\tDone ! ")
        else:
            result={}
            with open(inputFile,'rb') as f:
                results = pickle.load(f)
                print(results)
                bdtReweighter=results['model']
                doTransformation=results['doTransformation']
                if doTransformation:
                    data_x = bdtReweighter.transform( data_x )
                    mc_x   = bdtReweighter.transform( mc_x   )
        
        if doIronTrnsformation:
            f=plt.figure()
            _=plt.hist(data_[0],weights=data_w,bins=40,histtype='step',label='data')
            _=plt.hist(mc_x[0] ,weights=mc_w  ,bins=40,histtype='step',label='mc'  )
            plt.legend()
            plt.semilogy()
            f.savefig( saveBase + '/transformedVar_0.png' , bbox_inches='tight' )
            plt.close(f)
        
       
        ###        Plotting te obtained new weights and scale factors ! ##
        print("\t Obtaining then Scale factors for full MC ! ")
        MC_weights  = bdtReweighter.predictWeight(mc_x)
        wMC_total   = np.sum(MC_weights)
        wData_Total = np.sum(data_w)
        normalizationFactor= wData_Total/wMC_total
        print("\t     Normalization factor obtained as : ", normalizationFactor)
        MC_weights= MC_weights * normalizationFactor
        #MC_weights = data['mc']['weight']
        f=plt.figure()    
        plt.hist(data['mc']['weight'],histtype='step',label='raw',bins=80)
        plt.hist(data['mc']['weight_bined'],histtype='step',label='binScaled',bins=80)
        plt.hist(MC_weights,histtype='step',label='bdtSclaed',bins=80)
        plt.legend()
        plt.semilogy()
        plt.semilogx()
        plt.savefig(saveBase+'/'+yr+'_weightDistributions',bbox_inches='tight')
        plt.close(f)
        
        f=plt.figure()    
        plt.hist(data['mc']['weight_bined']/data['mc']['weight'],histtype='step',label='binScaled',bins=80)
        plt.hist(MC_weights/data['mc']['weight'],histtype='step',label='bdtSclaed',bins=80)
        plt.legend()
        plt.semilogy()
        plt.savefig(saveBase+'/'+yr+'_scaleFactorDistributions',bbox_inches='tight')
        plt.close(f)
     
        w_data      = data['data']['weight']
        w_mc        = data['mc']['weight'] 
        w_mc_binned = data['mc']['weight_bined'] 
        w_mc_bdt    = MC_weights 
        results={}
        
        ###                 YIELD Calculation             ###

        print("\t Obtaining the yields post fitting ! ")
        
        outDict['yields']={}
        tabYieldData=ptab.PrettyTable(['Year','Process','Events','Raw Yield','Reweighted v1','Reweighted v2' ] ) 
        wraw   =np.sum(dataStore[yr]['data']['data']['weight_v0'])
        wbinned=np.sum(dataStore[yr]['data']['data']['weight'])
        outDict['yields']['data']          =  {}
        outDict['yields']['data']['raw']   =str(np.round(wraw,3))
        outDict['yields']['data']['binned']=str(wbinned)
        tabYieldData.add_row([yr , 'data',len(dataStore[yr]['data']['data']['weight_v0']) ,
            outDict['yields']['data']['raw'] , outDict['yields']['data']['binned'] , -1 ])
       
        tabYield=ptab.PrettyTable(['Year','Process','Events','Raw Yield','Reweighted v1','Reweighted v2' ] ) 
        sumY={'raw':0.0 ,'binned':0.0 ,'bdt':0.0,'evts':0}
        for ky in dataStore[yr]['bkg']:
            wraw   =np.sum(dataStore[yr]['bkg'][ky]['weight_v0']*dataStore[yr]['bkg'][ky]['lumi'])
            wbinned=np.sum(dataStore[yr]['bkg'][ky]['weight'])
            xMC=np.stack([dataStore[yr]['bkg'][ky][k] for k in varForReweighting])
            wbdt = np.sum(bdtReweighter.predictWeight( xMC ))*normalizationFactor
            outDict['yields'][ky]          =  {}
            outDict['yields'][ky]['raw']   =str(np.round(wraw,3))
            outDict['yields'][ky]['binned']=str(np.round(wbinned,3))
            outDict['yields'][ky]['bdt']   =str(np.round(wbdt,3))
            if ky !='bkg':
                sumY['raw']    += wraw
                sumY['binned'] += wbinned
                sumY['bdt']    += wbdt
                sumY['evts']    += len(dataStore[yr]['bkg'][ky]['weight'])
            tabYield.add_row([yr , ky ,len(dataStore[yr]['bkg'][ky]['weight']),
                                outDict['yields'][ky]['raw'] , outDict['yields'][ky]['binned'] , outDict['yields'][ky]['bdt']   ])
        tabYield.add_row([yr , 'sum' , sumY['evts'] ,np.round(sumY['raw'],3) , np.round(sumY['binned'],3) , np.round(sumY['bdt'],3) ])
        print(tabYieldData)
        print(tabYield)
 
        ###                                    ROC TEST                     ###

        print("\tROC Test : ")
        results['validation_default']=bdtReweighter.validateClassifierModel( data_x , mc_x , w_mc )
        results['validation_binned'] =bdtReweighter.validateClassifierModel( data_x , mc_x , w_mc_binned )
        
        plt.figure()        

        _=plt.hist(results['validation_default']['data_score'],bins=np.linspace(0,1,20),color='k',histtype='step', density=True ,label='data')
        _=plt.hist(results['validation_default']['mc_score'],weights=w_mc       ,bins=np.linspace(0,1,20),color='g',histtype='step',linewidth=2,density=True ,label='mc raw weights')
        _=plt.hist(results['validation_default']['mc_score'],weights=w_mc_binned,bins=np.linspace(0,1,20),color='b',histtype='step',linewidth=2,density=True ,label='mc binned weights')
        _=plt.hist(results['validation_default']['mc_score'],weights=w_mc_bdt   ,bins=np.linspace(0,1,20),color='m',histtype='step',linewidth=2,density=True ,label='mc bdt weights')
        plt.legend()
        plt.savefig( saveBase+'/classifier_scores.png'  , bbox_inches='tight' )



        outDict['auc']={}
        tabROC=ptab.PrettyTable(['Category' , "ROC"] ) 
        f=plt.figure()
        roc=results['validation_default']['roc_raw']
        tabROC.add_row(["ROC with all Event Weight set to 1.0",np.round(results['validation_default']['auc_raw'],3)])
        outDict['auc']['raw'] = str(np.round(results['validation_default']['auc_raw'],3))
        
        roc=results['validation_default']['roc_pre']
        plt.plot( roc['tpr'], roc['fpr'],label='Default [ AUC : '+str(np.round(results['validation_default']['auc_pre'],3))+' ]')
        outDict['auc']['pre'] = str(np.round(results['validation_default']['auc_pre'],3))
        tabROC.add_row(["ROC with default Event Weights",np.round(results['validation_default']['auc_pre'],3)] )
        
        roc=results['validation_binned']['roc_pre']
        plt.plot( roc['tpr'], roc['fpr'],label='Binned Rw. [ AUC :'+str(np.round(results['validation_binned']['auc_pre'],3))+' ]')
        outDict['auc']['binned'] = str(np.round(results['validation_binned']['auc_pre'],3))
        tabROC.add_row(["ROC with Event Weights [ binned Rw. ]",np.round(results['validation_binned']['auc_pre'],3) ] )
        
        roc=results['validation_binned']['roc_post']
        plt.plot( roc['tpr'], roc['fpr'],label='BDT-Rewei. [ '+str(np.round(results['validation_default']['auc_post'],3))+' ]')
        outDict['auc']['bdt'] = str(np.round(results['validation_binned']['auc_post'],3))
        tabROC.add_row(["ROC with Event Weights [ BDT Rw. ]",np.round(results['validation_default']['auc_post'],3) ] )
        
        t=tabROC.get_string()
        print("\t\t"+t.replace("\n","\n\t\t"))

        foutname=saveBase+'/'+'validation_rocs.png'
        plt.legend()
        plt.savefig(foutname,bbox_inches='tight')


        ###                                   KS TEST                     ###

        results['ks_test']={   }
        var=''
        i=0
        makePlot=True
        print("Processing the variables for KS")
        for var in allVars:
            x_mc  =np.concatenate([dataStore[yr]['bkg'][ky][var] for ky in bkgToProcess])
            x_data=np.concatenate([dataStore[yr]['data'][ky][var] for ky in dataStore[yr]['data']])
            binEdges=np.linspace(varToBinMap[var][3] , varToBinMap[var][4] ,varToBinMap[var][2] )
            if makePlot:
                scl.plotDataMCComparison(var,
                                     x_mc,
                                     x_data,
                                     w_data,w_mc,w_mc_binned,w_mc_bdt,
                                     bins=binEdges,
                                     saveBase=saveBase)
            makePlot=False                                     
            results['ks_test'][var] = {
                            'pre' : stat.twoSample_KSTest(x_data, x_mc ,w_data, w_mc),
                            'binReweighted' : stat.twoSample_KSTest(x_data, x_mc ,w_data, w_mc_binned),
                            'bdtReweighted' : stat.twoSample_KSTest(x_data, x_mc ,w_data, w_mc_bdt)
                        }
        
        outDict['ks_test']={}
        if var in results['ks_test']:
            heads=list(results['ks_test'][var].keys())
            tab=ptab.PrettyTable(['Variable']+heads ) 
            for var in results['ks_test']:
               tab.add_row([var]+[ np.round(results['ks_test'][var][h],4) for h in heads  ])
               outDict['ks_test'][var]={ h : str(np.round(results['ks_test'][var][h],4)) for h in heads }
            t=tab.get_string()
            print("\t\t"+t.replace("\n","\n\t\t"))

        foutname=saveBase+'/'+'savedModel.pkl'
        if saveOutput:
            with open(foutname,'wb') as f:
                print("Saved the model to file :  ",foutname)
                bdtReweighter.saveModel(f)
                #pickle.dump(bdtReweighter,f)
        
        foutname=saveBase+'/'+'modelValidation.pkl'
        with open(foutname, 'wb') as f:
            pickle.dump(results,f)
        foutname=saveBase+'/'+'modelValidation.json'
        with open(foutname, 'w') as f:
            json.dump(outDict,f,indent=4)

if __name__=='__main__':
    main( )

