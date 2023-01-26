import matplotlib.pyplot as plt
import mplhep as hep
import uproot3 as urt
import ROOT
import json
import numpy as np
import os,argparse
from Util import lumiMap
import Util as utl
import trippleHiggsPlotter as plotterUtil

if __name__=='__main__':
    

    parser = argparse.ArgumentParser()
    parser.add_argument("--unblind", help="Unblind the Mgg spectrum", action='store_true' )
    args=parser.parse_args()

    unblind= args.unblind

    prefixBase='/home/aravind/cernbox/work/trippleHiggs/hhhTo4b2gamma/genAnalysis/python/analysis/results/plots/analysis/jan23/correlationsX/'
    fileDict={}
    with open('workarea/data/bdtNtuples/v7p2/filelist.json') as f:
        fileDict=json.load(f)
    yearsToProcess=['2018','2017','2016PreVFP','2016PostVFP','run2']
    yearsToProcess=['2018']
    bkgToProcess=[ 'ggBox1Bjet','ggBox2Bjet', 'ggBox','gJet20To40','gJet40ToInf']

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
        if unblind:
            rdataFrames[yr]['data']['data']= ROOT.RDataFrame(treeName, fileName)
            prefixBase+='/unblinded/'
        else:
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
    
    varlist=[
        'sumScore_4j','h1bb_mass','h2bb_mass','sumScore_3j','trihiggs_pt','scalarPtSum4b',
        'quadjet_0_deepJetScore','quadjet_2_deepJetScore','HH4bCosTheta_hhhF',
        'quadjet_1_deepJetScore','pThgg_overMgg','H1H2JetDrMin',
        'pTsubleadG_overMgg','quadjet_3_deepJetScore','pTh1leadJ_overMh1',
        'HggCosTheta_hhhF','H1bbCosTheta_hhhF','PhoJetMaxDr',
        'LeadJetDrMinWithOtherJets','CosThetaH1_hhhF','PhoJetMinDrOther',
        'PhoJetMinDr','pTh1subleadJ_overMh1',
        'pTh2subleadJ_overMh2','customSubLeadingPhotonIDMVA','h1bb_eta',
        'H1H2JetDrMax','h2bb_eta','customLeadingPhotonIDMVA','pTh2leadJ_overMh2',
        'r_HH','pTleadG_overMgg','PhoJetMaxDrOther','ttH_MET','h2_dijetSigmaMOverM',
        'absCosThetaH4bHgg','LeadJetDrMaxWithOtherJets','h1_dijetSigmaMOverM',
        'LeadJetAbsCosThetaMax','scalarPtSum4b2g','scalarPtSumHHH'
    ]
#    varlist=[
#        'LeadJetDrMinWithOtherJets','CosThetaH1_hhhF','PhoJetMinDrOther'
#    ]
    NVARS=len(varlist)
    dataForCorr={}
    for yr  in yearsToProcess:
        print("Processing Year : ",yr)
        dataForCorr[yr]={'sig':{},'bkg':{},'data':{}}
    
        _data=[]
        _weight=[]
        nBkg=0
        for ky in rdataFrames[yr]['bkg']:
            n=rdataFrames[yr]['bkg'][ky].Count().GetValue()
            if n<2000:
                print("\t !! n for bkg:" , ky,rdataFrames[yr]['bkg'][ky].Count().GetValue()  )
                nBkg+=n
            else:
                nBkg+=2000
            _dataDict=rdataFrames[yr]['bkg'][ky].Filter('rdfentry_ < 2000').AsNumpy(columns=['weight']+varlist)
            _data.append([_dataDict[ky] for ky in varlist ])
            _weight.append(_dataDict['weight'])
        dataForCorr[yr]['bkg']['x']=np.concatenate(_data,axis=1)
        dataForCorr[yr]['bkg']['weight']=np.concatenate(_weight)
        nStr=str(nBkg)
        _data=[]
        _weight=[]
        for ky in rdataFrames[yr]['sig']:
            _dataDict=rdataFrames[yr]['sig'][ky].Filter('rdfentry_ < '+nStr).AsNumpy(columns=['weight']+varlist)
            _data.append([_dataDict[ky] for ky in varlist ])
            _weight.append(_dataDict['weight'])
        dataForCorr[yr]['sig']['x']=np.concatenate(_data,axis=1)
        dataForCorr[yr]['sig']['weight']=np.concatenate(_weight)
    
    
        _data=[]
        _weight=[]
        for ky in rdataFrames[yr]['data']:
            _dataDict=rdataFrames[yr]['data'][ky].Filter('rdfentry_ < '+nStr).AsNumpy(columns=['weight']+varlist)
            _data.append([_dataDict[ky] for ky in varlist ])
            _weight.append(_dataDict['weight'])
        dataForCorr[yr]['data']['x']=np.concatenate(_data,axis=1)
        dataForCorr[yr]['data']['weight']=np.concatenate(_weight)

    correlationDictAllYears={}
    correlationMatrixesAllYears={}
    kk=0
    for yr  in yearsToProcess:
        kk+=1
        correlationDict={'sig':{},'bkg':{},'data':{}}
        correlationDict['varlist']={ i : varlist[i] for i in range(NVARS) }
        corrData=np.zeros((NVARS,NVARS))
        corrSig =np.zeros((NVARS,NVARS))
        corrBkg =np.zeros((NVARS,NVARS))
        print("Processing Year ",yr," [ ",kk,'/',len(yearsToProcess)," ]")
        for i,var in enumerate(varlist):
            correlationDict['data'][i]={}
            correlationDict['bkg'][i]={}
            correlationDict['sig'][i]={}
            if var=='LeadJetDrMinWithOtherJets':
                print("Data  : ",dataForCorr[yr]['data']['x'][i])
            for j in range(len(varlist)):
                correlationDict['data'][i][j]=utl.wei_corr(dataForCorr[yr]['data']['x'][i],
                           dataForCorr[yr]['data']['x'][j] ,
                           dataForCorr[yr]['data']['weight'])
                
                correlationDict['bkg'][i][j]=utl.wei_corr(dataForCorr[yr]['bkg']['x'][i],
                           dataForCorr[yr]['bkg']['x'][j] ,
                           dataForCorr[yr]['bkg']['weight'])
                
                correlationDict['sig'][i][j]=utl.wei_corr(dataForCorr[yr]['sig']['x'][i],
                           dataForCorr[yr]['sig']['x'][j] ,
                           dataForCorr[yr]['sig']['weight'])
                corrData[i][j]=correlationDict['data'][i][j]
                corrSig[i][j]=correlationDict['sig'][i][j]
                corrBkg[i][j]=correlationDict['bkg'][i][j]

            correlationDictAllYears[yr]=correlationDict
            correlationMatrixesAllYears[yr]={ 'sig':corrSig , 'bkg':corrBkg ,'data': corrData ,'dataMinusBkg' : corrData - corrBkg}

    
    for yr  in yearsToProcess:
        saveBase=prefixBase+'/'+yr+'/'
        os.system('mkdir -p '+saveBase)
        
        for key in correlationMatrixesAllYears[yr]:
            mat=correlationMatrixesAllYears[yr][key]
            f,ax=plt.subplots(figsize=(18,12))
            c=ax.imshow(mat,vmin=-1.0, vmax=1.0,cmap='tab20c',extent= (0, NVARS, 0, NVARS*4),aspect=0.14)
            varlist_y=[varlist[n]+'  ['+str(n)+'] ' for n in range(NVARS)]
            varlist_x=['['+str(n)+']' for n in range(NVARS)]
            t=ax.set_yticks(np.arange(0.0,NVARS,1)*4,varlist_y)
            t=ax.set_xticks(np.arange(0.0,NVARS,1),varlist_x,ha='left')
            ax2 = f.add_axes([ax.get_position().x1+0.005,ax.get_position().y0,
                               0.02,ax.get_position().height])
            plt.colorbar(c,cax=ax2)
            ax.grid(color='k',alpha=1)
            
            f.savefig(saveBase+'/correlations_'+key+'.pdf',bbox_inches='tight')
            f.savefig(saveBase+'/correlations_'+key+'.png',bbox_inches='tight')
            plt.close(f)
