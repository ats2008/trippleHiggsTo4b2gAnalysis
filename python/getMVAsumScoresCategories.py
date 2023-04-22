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
sigRegionCut='CMS_hgg_mass >= 115 && CMS_hgg_mass <= 135.0'

def main():
 
    parser = argparse.ArgumentParser()
    parser.add_argument('-v',"--version", help="Version of the specific derivation ",default='')
    parser.add_argument('-i',"--inputFile", help="Input File",default=None)
    parser.add_argument('-y',"--year", help="Year",default='run2')
    parser.add_argument("-o","--dest", help="destination To Use", default='workarea/results/plots/tmp/' )
    parser.add_argument("-c","--cuts", help="list of cuts to be applied", default=None )
    parser.add_argument("--cutstr", help="additional cutsrting to be applied", default=None )
    parser.add_argument("--doSR", help="Do Signal Region",default=False,action='store_true')
    args = parser.parse_args()
    version = args.version
    doIronTrnsformation=False
    inputFile  = args.inputFile
    saveOutput = True
    cutsToApply=[]
    cutsKey=[]
    cutIdx=0
    txt=[]
    if args.cuts:
        with open(args.cuts,'r') as f:
            txt_=f.readlines()
            for l in txt_:
                txt.append(l[:-1])
    cutStrs=[]
    if args.cutstr:
        items=args.cutstr.split("][")
        txt+=items
    for l in txt:
        if l[0]=='#':
            continue
        key='cut '+str(cutIdx)
        cut=l.split(',')[0]
        if ',' in l :
            key=l.split(',')[1]
        cutsToApply.append(cut)
        cutsKey.append(key)
        cutIdx+=1
    print("Cuts being applied : ")
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
    fullYields={}
    for yr  in yearsToProcess:
        fullYields[yr]={}
        rdataFrames[yr]={'sig':{},'bkg':{},'data':{}}
        try:
            fileName=fileDict[yr]['sig']['ggHHH']
        except:
            print(yr)
        treeName = "trees/ggHHH_125_13TeV"
        rdataFrames[yr]['sig']['ggHHH']= ROOT.RDataFrame(treeName, fileName).Filter(sigRegionCut)
        for cut,cutKey in zip(cutsToApply,cutsKey):
            rdataFrames[yr]['sig']['ggHHH'] = rdataFrames[yr]['sig']['ggHHH'].Filter(cut,cutKey)
        fullYields[yr]['sig']= np.sum( rdataFrames[yr]['sig']['ggHHH'].AsNumpy(['weight_bdt'])['weight_bdt'] )   
        print("Registering datset : ggHHH , ",yr," with tree",treeName)
        print("\t\t ",fileName)   
        allCutsReport = rdataFrames[yr]['sig']['ggHHH'].Report()
        allCutsReport.Print()          
        
        ky=list(fileDict[yr]['data'].keys())[0]
        fileName=fileDict[yr]['data'][ky]
        treeName = "trees/Data_13TeV_TrippleHTag_0"
        rdataFrames[yr]['data']['data']= ROOT.RDataFrame(treeName, fileName).Filter(data_blind).Filter(blind)
        for cut,cutKey in zip(cutsToApply,cutsKey):
            rdataFrames[yr]['data']['data'] = rdataFrames[yr]['data']['data'].Filter(cut,cutKey)
        fullYields[yr]['data']= np.sum( rdataFrames[yr]['data']['data'].AsNumpy(['weight_bdt'])['weight_bdt'] )   
            
        print("Registering datset : data , ",yr," with tree",treeName)
        print("\t\t ",fileName)   
        allCutsReport = rdataFrames[yr]['data']['data'].Report()
        allCutsReport.Print()          
        
        rdataFrames[yr]['bkg']={}
        treeName = "trees/bkg_13TeV_TrippleHTag_0"
        for bkg in bkgToProcess:
            if bkg not in fileDict[yr]['bkg']:
                print()
                print("FILE NOT FOUNG FOR : ",yr," background ",bkg)
                print()
                continue
            fileName = fileDict[yr]['bkg'][bkg]
            rdataFrames[yr]['bkg'][bkg]=ROOT.RDataFrame(treeName, fileName).Filter(sigRegionCut)
            for cut,cutKey in zip(cutsToApply,cutsKey):
                rdataFrames[yr]['bkg'][bkg] = rdataFrames[yr]['bkg'][bkg].Filter(cut,cutKey)
            
            fullYields[yr]['bkg']= np.sum( rdataFrames[yr]['bkg']['bkg'].AsNumpy(['weight_bdt'])['weight_bdt'] )   
            print("Registering datset : ",bkg," , ",yr," with tree",treeName)
            print("\t\t ",fileName)   
            
            allCutsReport = rdataFrames[yr]['bkg'][bkg].Report()
            allCutsReport.Print()          
    
    varsToGet=['weight_bdt','weight_binned','weight_v0','lumi','genRecoCategory']
    varsToGet+=[
            'hhhNR_mva_2023_v0',
            'sumScore_4j',
            'sumScore_3j'
        ]

    mvaCuts   =np.array([0.0,0.9 , 0.98,0.99,0.999 , 1.0  ])
    #sumTagCuts=np.array([0.0, 2.0, 2.8 , 4.0])
    sumTagCuts=np.array([0.0,4.0])
    
    #mvaCuts   =np.arange(0.8,1.0+0.00005,0.01)
    #sumTagCuts=np.arange(0.0,4.0+0.00005,0.05)
    
    print(f"{mvaCuts.shape=}")
    print(f"{sumTagCuts.shape=}")
    for yr in yearsToProcess:
        saveBase=prefixBase+'/'+yr+'/scategSumScoreVsMVA/'
        os.system('mkdir -p '+saveBase)
        outDict={}
        outDict['fullYields'] = {
                'data' : float(fullYields[yr]['data']),
                'sig' : float(fullYields[yr]['sig']),
                'bkg' : float(fullYields[yr]['bkg']),
            }
        scoreMaps={}

        for dset in ['bkg','data','sig']:
            outDict[dset]={}
            tabYieldData=ptab.PrettyTable(['cut','yield'] ) 
            for ky in rdataFrames[yr][dset]:
                print("Processing ",ky)
                varData=rdataFrames[yr][dset][ky].AsNumpy(varsToGet)
                
                vals=np.zeros( ( mvaCuts.shape[0]-1 , sumTagCuts.shape[0] -1) )    
                print(f"{vals.shape=}")

                for i in range(mvaCuts.shape[0]-1):
                    mvaMask=varData['hhhNR_mva_2023_v0'] >= mvaCuts[i]
                    mvaMask=np.logical_and(mvaMask, varData['hhhNR_mva_2023_v0'] < mvaCuts[i+1] )
                    for j in range(sumTagCuts.shape[0]-1):
                        btagMask = np.ones(varData['sumScore_4j'].shape ,dtype='bool')
                        btagMask = np.logical_and( btagMask ,  varData['sumScore_4j'] >= sumTagCuts[ j   ] )
                        btagMask = np.logical_and( btagMask ,  varData['sumScore_4j'] <  sumTagCuts[ j+1 ] )
                        mask=np.logical_and(mvaMask,btagMask)
                        vals[i][j] = np.sum(varData['weight_bdt'][mask])

                                
                #print("Selcted numbers : ",vals)
                norm=np.sum(varData['weight_bdt'])
                scoreMaps[dset]=vals
                #print(scoreMaps[dset])
                vals=scoreMaps[dset]/norm
                #print(scoreMaps[dset])
                outDict[dset]['evnts']={}
                outDict[dset]['fractions']={}
                print(f"{scoreMaps[dset].shape=}")
                for i in range(vals.shape[0]):
                    for j in range(vals.shape[1]):
                        tagName=f'MVA_{mvaCuts[i]:.2f}To_{mvaCuts[i+1]:.2f}_SS_{sumTagCuts[j]:.2f}To_{sumTagCuts[j+1]:.2f}'.replace('.','p')
                        #tagName=f'MVAmax_{mvaCuts[i]:.2f}_SSmin_{sumTagCuts[j]:.2f}'.replace('.','p')
                        outDict[dset]['evnts'][tagName] = str(np.round(scoreMaps[dset][i][j] , 4 ))
                        outDict[dset]['fractions'][tagName] = str(np.round( vals[i][j] ,4 ))
                        tabYieldData.add_row( [ tagName , str(np.round( vals[i][j] ,4  ))] )
                
                ## PLOTING THE FRACTIONS
                if dset=='sig':
                    f,ax=plt.subplots(1,1,figsize=(8,10))
                    c=ax.matshow(vals,aspect=1.0,cmap='cool')
                    for (i, j), z in np.ndenumerate(vals):
                        ax.text(j, i, '{:0.3f}'.format(z), ha='center', va='center',fontsize=10)
                    _=ax.set_yticks(ticks=[ i-0.5 for i in range(mvaCuts.shape[0]-1) ],
                                    labels=[np.round(mvaCuts[i],3) for i in range(mvaCuts.shape[0]-1)])
                    _=ax.set_xticks(ticks=[ i-0.5 for i in range(sumTagCuts.shape[0]) ],
                                    labels=[np.round(i,3) for i in sumTagCuts])
                    ax.grid(which='major',color='k',alpha=1)
                    ax.set_ylabel('MVA Cut')
                    ax.set_xlabel('SumScore Cut')
                    ax.set_xlim([-0.5 ,sumTagCuts.shape[0]-1-0.5])
                    ax.set_ylim([-0.5 ,mvaCuts.shape[0]-1-0.5])
                    ax2 = f.add_axes([ax.get_position().x1+0.005,ax.get_position().y0,0.02,ax.get_position().height])
                    plt.colorbar(c,cax=ax2)
                    ax2.set_yticks([])
                    f.savefig(saveBase+'/selectionEfffciency_'+dset+'.png',bbox_inches='tight')
                    plt.close(f)

                if True:
                    ## PLOTING THE COUNTS
                    f,ax=plt.subplots(1,1,figsize=(8,10))
                    c=ax.matshow(scoreMaps[dset],cmap='cool',aspect=1.0)
                    for (i, j), z in np.ndenumerate(scoreMaps[dset]):
                        ax.text(j, i, '{:0.3f}'.format(z), ha='center', va='center',fontsize=10)
                    _=ax.set_yticks(ticks=[ i-0.5 for i in range(mvaCuts.shape[0]-1) ],
                                    labels=[np.round(mvaCuts[i],3) for i in range(mvaCuts.shape[0]-1)])
                    _=ax.set_xticks(ticks=[ i-0.5 for i in range(sumTagCuts.shape[0]) ],
                                    labels=[np.round(i,3) for i in sumTagCuts])
                    ax.grid(which='major',color='k',alpha=1)
                    ax.set_ylabel('MVA Cut')
                    ax.set_xlabel('SumScore Cut')
                    ax.set_xlim([-0.5 ,sumTagCuts.shape[0]-1-0.5])
                    ax.set_ylim([-0.5 ,mvaCuts.shape[0]-1-0.5])
                    ax2 = f.add_axes([ax.get_position().x1+0.005,ax.get_position().y0,0.02,ax.get_position().height])
                    plt.colorbar(c,cax=ax2)
                    ax2.set_yticks([])
                    f.savefig(saveBase+'/selectionCounts_'+dset+'.png',bbox_inches='tight')
                    plt.close(f)


        #for i,dset in enumerate(['data','bkg']):
        for i,dset in enumerate(['bkg']):
            print(scoreMaps['sig'])
            print(scoreMaps[dset])
            significance = scoreMaps['sig']*1/np.sqrt( scoreMaps[dset] +1e-9)
            #significance[ scoreMaps[dset] == 0.0 ] = np.max(significance )
            print(significance)
            outDict[dset]['significance']={}

            f,ax=plt.subplots(1,1,figsize=(8,10))
            c=ax.matshow(significance,cmap='cool',aspect=1.0)
            tabYieldData=ptab.PrettyTable(['cut','significance'] ) 

            for (i, j), z in np.ndenumerate(significance):
                ax.text(j, i, '{:0.3f}'.format(z), ha='center', va='center',fontsize=10)
                tagName=f'MVA_{mvaCuts[i]:.2f}To_{mvaCuts[i+1]:.2f}_SS_{sumTagCuts[j]:.2f}To_{sumTagCuts[j+1]:.2f}'.replace('.','p')
                #tagName=f'MVAmax_{mvaCuts[i]:.2f}_SSmin_{sumTagCuts[j]:.2f}'.replace('.','p')
                outDict[dset]['significance'][tagName] = str(np.round( significance[i][j] , 4))

            _=ax.set_yticks(ticks=[ i-0.5 for i in range(mvaCuts.shape[0]-1) ],
                            labels=[np.round(mvaCuts[i],3) for i in range(mvaCuts.shape[0]-1)])
            _=ax.set_xticks(ticks=[ i-0.5 for i in range(sumTagCuts.shape[0]) ],
                            labels=[np.round(i,3) for i in sumTagCuts])

            ax.grid(which='major',color='k',alpha=1)
            ax.set_ylabel('MVA Cut')
            ax.set_xlabel('SumScore Cut')
            ax.set_xlim([-0.5 ,sumTagCuts.shape[0]-1-0.5])
            ax.set_ylim([-0.5 ,mvaCuts.shape[0]-1-0.5])
            ax2 = f.add_axes([ax.get_position().x1+0.005,ax.get_position().y0,0.02,ax.get_position().height])
            plt.colorbar(c,cax=ax2)
            ax2.set_yticks([])
            f.savefig(saveBase+'/cutSignificance'+dset+'.png',bbox_inches='tight')
            plt.close(f)
        print("Output saved at : ",saveBase)
        foutname=saveBase+'/'+'yields.json'
        with open(foutname, 'w') as f:
            json.dump(outDict,f,indent=4)

if __name__=='__main__':
    main( )

