import matplotlib.pyplot as plt
import mplhep as hep
import uproot as urt
import ROOT
import json,sys,os,argparse
import numpy as np
import Util as utl
import pickle,itertools
import prettytable as ptab

import scaleFactorUtil as scl
import statUtil as stat
import hep_ml as hepml
import prettytable


hep.style.use("CMS")
data_blind='CMS_hgg_mass < 115 || CMS_hgg_mass > 135.0'

def main():
 
    parser = argparse.ArgumentParser()
    parser.add_argument('-v',"--version", help="Version of the specific derivation ",default='')
    parser.add_argument('-i',"--inputFile", help="Input File",default=None)
    parser.add_argument('-y',"--year", help="Year",default='2018')
    parser.add_argument("-o","--dest", help="destination To Use", default='workarea/results/plots/tmp/' )
    parser.add_argument("-c","--cuts", help="list of cuts to be applied", default=None )
    parser.add_argument('-t',"--transformData", help="trnsoform data using IronTransformer ",default=False,action='store_true')
    parser.add_argument("--doSR", help="Do Signal Region",default=False,action='store_true')
    args = parser.parse_args()
    version = args.version
    doIronTrnsformation=False
    inputFile  = args.inputFile
    saveOutput = True
    cutsToApply=[]
    if args.cuts:
        with open(args.cuts,'r') as f:
            txt=f.readlines()
            for l in txt:
                if l[0]=='#':
                    continue
                cutsToApply.append(l[:-1])
        print("Cuts Being applied : ")
        for cut in cutsToApply:
            print("\t -> ",cut)

    prefixBase=args.dest
    outDict={}
    blind=data_blind
    if args.doSR:
        blind='CMS_hgg_mass >= 115 && CMS_hgg_mass <= 135.0'
        prefixBase+='/SR/'
    print("Processing file list : ",args.inputFile)
    if inputFile:
        saveOutput=False
    fileDict={}
    with open(args.inputFile) as f:
        fileDict=json.load(f)
    
   
    yearsToProcess_all=['2018','2017','2016PreVFP','2016PostVFP','run2']
    yearsToProcess=yearsToProcess_all
    
    if 'all' not in args.year:
        yearsToProcess=[]
        
        for yr in args.year.split(","):
            if yr not in yearsToProcess_all:
                print(yr," not in catalogue. Skipping !! ")
                continue
            yearsToProcess.append(yr)

    bkgToProcess=[ 'bkg',
             #      'ggBox1Bjet',
             #      'ggBox2Bjet', 
             #      'ggBox', 
             #      'gJet20To40',
             #      'gJet40ToInf'
                   ]
    varToBinMap={}
    
    print()
    print("Outputs are stored in ",prefixBase)
    print()
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
        rdataFrames[yr]['sig']['ggHHH']= ROOT.RDataFrame(treeName, fileName).Filter(blind)
        for cut in cutsToApply:
            rdataFrames[yr]['sig']['ggHHH'] = rdataFrames[yr]['sig']['ggHHH'].Filter(cut)
        
        print("Registering datset : ggHHH , ",yr," with tree",treeName)
        print("\t\t ",fileName)   
        ky=list(fileDict[yr]['data'].keys())[0]
        fileName=fileDict[yr]['data'][ky]
        treeName = "trees/Data_13TeV_TrippleHTag_0"
        rdataFrames[yr]['data']['data']= ROOT.RDataFrame(treeName, fileName).Filter(data_blind).Filter(blind)
        for cut in cutsToApply:
            rdataFrames[yr]['data']['data'] = rdataFrames[yr]['data']['data'].Filter(cut)
            
        print("Registering datset : data , ",yr," with tree",treeName)
        print("\t\t ",fileName)   
        
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
            for cut in cutsToApply:
                rdataFrames[yr]['bkg'][bkg] = rdataFrames[yr]['bkg'][bkg].Filter(cut)
            
            print("Registering datset : ",bkg," , ",yr," with tree",treeName)
            print("\t\t ",fileName)   
                       
    
    varsToGet=['weight_bdt','weight_binned','weight_v0','lumi','genRecoCategory']
    varsToGet+=[
            'quadjet_0_isBtagLoose',
            'quadjet_1_isBtagLoose',
            'quadjet_2_isBtagLoose',
            'quadjet_3_isBtagLoose',
            'quadjet_0_isBtagMedium',
            'quadjet_1_isBtagMedium',
            'quadjet_2_isBtagMedium',
            'quadjet_3_isBtagMedium',
            'quadjet_0_isBtagTight',
            'quadjet_1_isBtagTight',
            'quadjet_2_isBtagTight',
            'quadjet_3_isBtagTight',
        ]
    
    allProducts=itertools.product(['L','M','T'],repeat=4)

    keys=[[],[]]
    allprds=[]
    nKeys=0
    for  prd in allProducts:
        keyStr= "".join(prd)
        allprds.append(keyStr)
        keys[0].append(keyStr[0:2])
        keys[1].append(keyStr[2:])
        nKeys+=1 
    nKeys= int(np.sqrt(nKeys))
    keys[0]=np.array(keys[0]).reshape(nKeys,nKeys)
    keys[1]=np.array(keys[1]).reshape(nKeys,nKeys)
    for yr in yearsToProcess:
        saveBase=prefixBase+'/'+yr+'/bTagEfficiencies/'
        os.system('mkdir -p '+saveBase)
        outDict={}
        scoreMaps={}
        for dset in ['bkg','data','sig']:
            outDict[dset]={}
            tabYieldData=ptab.PrettyTable(['cut','yield'] ) 
            for ky in rdataFrames[yr][dset]:
                print("Processing ",ky)
                varData=rdataFrames[yr][dset][ky].AsNumpy(varsToGet)

                idMap={'b1':{},'b2':{},'b3':{},'b4':{}}
                idMap['b1']['T'] = varData['quadjet_0_isBtagTight'] >0.2
                idMap['b2']['T'] = varData['quadjet_1_isBtagTight'] >0.2
                idMap['b3']['T'] = varData['quadjet_2_isBtagTight'] >0.2
                idMap['b4']['T'] = varData['quadjet_3_isBtagTight'] >0.2
                idMap['b1']['M'] = varData['quadjet_0_isBtagMedium'] >0.2
                idMap['b2']['M'] = varData['quadjet_1_isBtagMedium'] >0.2
                idMap['b3']['M'] = varData['quadjet_2_isBtagMedium'] >0.2
                idMap['b4']['M'] = varData['quadjet_3_isBtagMedium'] >0.2
                idMap['b1']['L'] = varData['quadjet_0_isBtagLoose'] >0.2
                idMap['b2']['L'] = varData['quadjet_1_isBtagLoose'] >0.2
                idMap['b3']['L'] = varData['quadjet_2_isBtagLoose'] >0.2
                idMap['b4']['L'] = varData['quadjet_3_isBtagLoose'] >0.2
                vals=[]

                for keyStr in allprds:
                    mask=np.ones(idMap['b1'][keyStr[0]].shape[0],dtype=bool)
                    mask=np.logical_and(mask,idMap['b1'][keyStr[0]])
                    mask=np.logical_and(mask,idMap['b2'][keyStr[1]])
                    mask=np.logical_and(mask,idMap['b3'][keyStr[2]])
                    mask=np.logical_and(mask,idMap['b4'][keyStr[3]])
                    vals.append(np.sum(varData['weight_bdt'][mask]))
                #print("Selcted numbers : ",vals)
                norm=np.sum(varData['weight_bdt'])
                scoreMaps[dset]=np.array(vals).reshape(nKeys,nKeys)
                #print(scoreMaps[dset])
                vals=scoreMaps[dset]/norm
                #print(scoreMaps[dset])
                outDict[dset]['evnts']={}
                outDict[dset]['fractions']={}
                for i in range(vals.shape[0]):
                    for j in range(vals.shape[0]):
                        outDict[dset]['evnts'][keys[0][i,0]+keys[1][0,j]] = str(np.round(scoreMaps[dset][i][j] , 4  ))
                        outDict[dset]['fractions'][keys[0][i,0]+keys[1][0,j]] = str(np.round( vals[i][j] ,4  ))
                        tabYieldData.add_row( [ keys[0][i,0]+keys[1][0,j] , str(np.round( vals[i][j] ,4  )) ] )
                ## PLOTING THE FRACTIONS
                f,ax=plt.subplots(figsize=(7,7))
                c=ax.matshow(vals,cmap='Pastel1')
                for (i, j), z in np.ndenumerate(vals):
                    ax.text(j, i, '{:0.3f}'.format(z), ha='center', va='center',fontsize=10)
                ax.set_xticks(np.arange(nKeys)+0.5,np.arange(nKeys)+0.5,minor=True,color='w',fontsize=0.1)
                ax.set_yticks(np.arange(nKeys)+0.5,np.arange(nKeys)+0.5,minor=True,color='w',fontsize=0.1)
                ax.set_xticks(np.arange(nKeys),labels=keys[0][:,0],fontsize=12)
                ax.set_yticks(np.arange(nKeys),labels=keys[1][0,:],fontsize=12)
                ax.set_xlabel(dset)
                ax.grid(which='minor')
                ax2 = f.add_axes([ax.get_position().x1+0.005,ax.get_position().y0,0.02,ax.get_position().height])
                plt.colorbar(c,cax=ax2)
                ax2.set_yticks([])
                f.savefig(saveBase+'/btagSelectionEffciency_'+dset+'.png',bbox_inches='tight')
                plt.close(f)
               
                ## PLOTING THE FRACTIONS
                f,ax=plt.subplots(figsize=(7,7))
                c=ax.matshow(scoreMaps[dset],cmap='Pastel1')
                for (i, j), z in np.ndenumerate(scoreMaps[dset]):
                    ax.text(j, i, '{:0.3f}'.format(z), ha='center', va='center',fontsize=10)
                ax.set_xticks(np.arange(nKeys)+0.5,np.arange(nKeys)+0.5,minor=True,color='w',fontsize=0.1)
                ax.set_yticks(np.arange(nKeys)+0.5,np.arange(nKeys)+0.5,minor=True,color='w',fontsize=0.1)
                ax.set_xticks(np.arange(nKeys),labels=keys[0][:,0],fontsize=12)
                ax.set_yticks(np.arange(nKeys),labels=keys[1][0,:],fontsize=12)
                ax.set_xlabel(dset)
                ax.grid(which='minor')
                ax2 = f.add_axes([ax.get_position().x1+0.005,ax.get_position().y0,0.02,ax.get_position().height])
                plt.colorbar(c,cax=ax2)
                ax2.set_yticks([])
                f.savefig(saveBase+'/btagSelectionCounts_'+dset+'.png',bbox_inches='tight')
                
                plt.close(f)
            #print(tabYieldData)

        f,axSet=plt.subplots(1,2,figsize=(12,6))

        for i,dset in enumerate(['data','bkg']):
            ax=axSet[i]
            significance = scoreMaps['sig']*1e3/np.sqrt(scoreMaps['sig'] + scoreMaps[dset] )
            outDict[dset]['significance']={}
            tabYieldData=ptab.PrettyTable(['cut','significance'] ) 
            for i in range(vals.shape[0]):
                for j in range(vals.shape[0]):
                    outDict[dset]['significance'][keys[0][i,0]+keys[1][0,j]] = str(np.round( significance[i][j] , 4))
                    tabYieldData.add_row( [ keys[0][i,0]+keys[1][0,j] , str(np.round( significance[i][j] ,4  )) ] )
            c=ax.matshow(significance,cmap='Pastel1')
            
            for (i, j), z in np.ndenumerate(significance):
                ax.text(j, i, '{:0.3f}'.format(z), ha='center', va='center',fontsize=10)
            ax.set_xticks(np.arange(nKeys)+0.5,np.arange(nKeys)+0.5,minor=True,color='w',fontsize=0.1)
            ax.set_yticks(np.arange(nKeys)+0.5,np.arange(nKeys)+0.5,minor=True,color='w',fontsize=0.1)
            ax.set_xticks(np.arange(nKeys),labels=keys[0][:,0],fontsize=12)
            ax.set_yticks(np.arange(nKeys),labels=keys[1][0,:],fontsize=12)
            ax.set_xlabel(dset)
            ax2 = f.add_axes([ax.get_position().x1+0.005,ax.get_position().y0,0.02,ax.get_position().height])
            plt.colorbar(c,cax=ax2)
            ax2.set_yticks([])
            print(tabYieldData)
        f.savefig(saveBase+'/btagSignificance_'+dset+'.png',bbox_inches='tight')
        print("Output saved at : ",saveBase)
        foutname=saveBase+'/'+'pu_modelYields.json'
        with open(foutname, 'w') as f:
            json.dump(outDict,f,indent=4)
if __name__=='__main__':
    main( )


