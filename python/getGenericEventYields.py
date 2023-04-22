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
data_blind='CMS_hgg_mass < 115 || CMS_hgg_mass > 135.0'

def main():
 
    parser = argparse.ArgumentParser()
    parser.add_argument('-v',"--version", help="Version of the specific derivation ",default='')
    parser.add_argument('-i',"--inputFile", help="Input File",default=None)
    parser.add_argument('-y',"--year", help="Year",default='2018')
    parser.add_argument("-o","--dest", help="destination To Use", default='workarea/results/plots/tmp/' )
    parser.add_argument("-c","--cuts", help="list of cuts to be applied", default=None )
    parser.add_argument("--scaleFactor", help="additional scales to be multiplied  as json { year { bkg : { deset : <SF> } }}", default=None )
    parser.add_argument("--cutstr", help="additional cutsrting to be applied", default=None )
    parser.add_argument('-t',"--transformData", help="trnsoform data using IronTransformer ",default=False,action='store_true')
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
    miscScales={}
    if args.scaleFactor:
        with open(args.scaleFactor) as f:
            miscScales=json.load(f)
            print("Loaded the scales json : ",json.dumps(miscScales,indent=4))
    if args.cuts:
        with open(args.cuts,'r') as f:
            txt_=f.readlines()
            for l in txt_:
                txt.append(l[:-1])
        print("Cuts Being applied : ")
        for cut in cutsToApply:
            print("\t -> ",cut)
    
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
    print("Cuts Being applied : ")
    for cut in cutsToApply:
        print("\t -> ",cut)


    outDict={}
    blind=data_blind
    if args.doSR:
        blind='CMS_hgg_mass >= 115 && CMS_hgg_mass <= 135.0'
    print("Processing file list : ",args.inputFile)
    if inputFile:
        saveOutput=False
    fileDict={}
    with open(args.inputFile) as f:
        fileDict=json.load(f)
    
    prefixBase=args.dest
   
    yearsToProcess_all=['2018','2017','2016PreVFP','2016PostVFP','run2','extras']
    yearsToProcess=yearsToProcess_all
    
    if 'all' not in args.year:
        yearsToProcess=[]
        
        for yr in args.year.split(","):
            if yr not in yearsToProcess_all:
                print(yr," not in catalogue. Skipping !! ")
                continue
            yearsToProcess.append(yr)

    print()
    print("Outputs are stored in ",prefixBase)
    print()
    # We read the tree from the file and create a RDataFrame, a class that
    # allows us to interact with the data contained in the tree.
    rdataFrames={}
    for yr  in yearsToProcess:
        rdataFrames[yr]={'sig':{},'bkg':{},'data':{}}
        print( rdataFrames[yr] )
        if 'sig' in fileDict[yr]:
            try:
                fileName=fileDict[yr]['sig']['ggHHH']
            except:
                print(yr)
                print("Exiting")
                exit(1)
            treeName = "trees/ggHHH_125_13TeV"
            rdataFrames[yr]['sig']['ggHHH']= ROOT.RDataFrame(treeName, fileName).Filter(blind)
            for cut,cutKey in zip(cutsToApply,cutsKey):
                rdataFrames[yr]['sig']['ggHHH'] = rdataFrames[yr]['sig']['ggHHH'].Filter(cut,cutKey)
            
            print("Registering datset : ggHHH , ",yr," with tree",treeName)
            print("\t\t ",fileName)   
        else:
            rdataFrames[yr]['sig']={}

        if 'data' in fileDict[yr]:
            ky=list(fileDict[yr]['data'].keys())[0]
            fileName=fileDict[yr]['data'][ky]
            treeName = "trees/Data_13TeV_TrippleHTag_0"
            rdataFrames[yr]['data']['data']= ROOT.RDataFrame(treeName, fileName).Filter(data_blind).Filter(blind)
            for cut,cutKey in zip(cutsToApply,cutsKey):
                rdataFrames[yr]['data']['data'] = rdataFrames[yr]['data']['data'].Filter(cut,cutKey)
            print("Registering datset : data , ",yr," with tree",treeName)
            print("\t\t ",fileName)   
        else:
            rdataFrames[yr]['data']={}

        
        rdataFrames[yr]['bkg']={}
        treeName = "trees/bkg_13TeV_TrippleHTag_0"
        for bkg in fileDict[yr]['bkg']:
            if bkg not in fileDict[yr]['bkg']:
                print()
                print("FILE NOT FOUNG FOR : ",yr," background ",bkg)
                print()
                continue
            fileName = fileDict[yr]['bkg'][bkg]
            rdataFrames[yr]['bkg'][bkg]=ROOT.RDataFrame(treeName, fileName).Filter(blind)
            for cut,cutKey in zip(cutsToApply,cutsKey):
                rdataFrames[yr]['bkg'][bkg] = rdataFrames[yr]['bkg'][bkg].Filter(cut,cutKey)
            
            print("Registering datset : ",bkg," , ",yr," with tree",treeName)
            print("\t\t ",fileName)   
                       
    
    # reading the vars from file    
    dataStore={}
    varsToGet=['weight_v0','weight','weight_binned','weight_bdt','lumi']
    
    for yr in yearsToProcess:
        models={}
        print()
        print()
        print("Processing for year : ",yr)
        saveBase=prefixBase+'/'+yr+'/'
        os.system('mkdir -p '+saveBase)
        dataStore[yr]={}
        bkgStack=None
        
        outDict['yields']={}
        tabYieldData=ptab.PrettyTable(['Year','Process','Events','Raw Yield','Bin Reweighted','BDT Reweighted' ] ) 
        
        dataStore[yr]['bkg']={}
        for ky in rdataFrames[yr]['bkg']:
            dataStore[yr]['bkg'][ky]=rdataFrames[yr]['bkg'][ky].AsNumpy(varsToGet)

        dataStore[yr]['sig']={}
        if 'sig' in rdataFrames[yr]:
           for ky in rdataFrames[yr]['sig']:
                dataStore[yr]['sig'][ky]=rdataFrames[yr]['sig'][ky].AsNumpy(varsToGet)
        
           for ky in dataStore[yr]['sig']:
                wraw   =np.sum(dataStore[yr]['sig'][ky]['weight_v0']*dataStore[yr]['sig'][ky]['lumi'])
                wbinned=np.sum(dataStore[yr]['sig'][ky]['weight_binned'])
                wbdt   =np.sum(dataStore[yr]['sig'][ky]['weight_bdt'])
                outDict['yields'][ky]          =  {}
                outDict['yields'][ky]['Events']=str(len( dataStore[yr]['sig'][ky]['weight_bdt'] ))
                scf=1.0
                if yr in miscScales:
                    if 'sig' in miscScales[yr]:
                        if ky in miscScales[yr]['sig'] :
                            scf=miscScales[yr]['sig'][ky]
                outDict['yields'][ky]['raw']   =str(np.round(wraw*scf,5))
                outDict['yields'][ky]['binned']=str(np.round(wbinned*scf,5))
                outDict['yields'][ky]['bdt']   =str(np.round(wbdt*scf,5))
                tabYieldData.add_row([yr , ky ,
                                      outDict['yields'][ky]['Events'] ,
                                      outDict['yields'][ky]['raw'] ,
                                      outDict['yields'][ky]['binned'] ,
                                      outDict['yields'][ky]['bdt'] ])
        
        dataStore[yr]['data']={}
        if 'data' in rdataFrames[yr]['data']:
            dataStore[yr]['data']['data']=rdataFrames[yr]['data']['data'].AsNumpy(varsToGet)
            wraw   =np.sum(dataStore[yr]['data']['data']['weight_v0'])
            wbinned=np.sum(dataStore[yr]['data']['data']['weight_binned'])
            wbdt   =np.sum(dataStore[yr]['data']['data']['weight_bdt'])
            outDict['yields']['data']          =  {}
            outDict['yields']['data']['Events']= str(len(dataStore[yr]['data']['data']['weight_v0']))
            outDict['yields']['data']['raw']   =str(np.round(wraw,3))
            outDict['yields']['data']['binned']=str(wbinned)
            outDict['yields']['data']['bdt']=str(wbdt)
            tabYieldData.add_row([yr , 'data' ,
                                  outDict['yields']['data']['Events'] , 
                                  outDict['yields']['data']['raw'] , 
                                  outDict['yields']['data']['binned'] , 
                                  outDict['yields']['data']['bdt'] ])

        
        data={
                'data':{'x':None,'weight':None},
                'mc'  :{'x':None,'weight':None}
             }
        
        isFirst=True
        
        tabYield=ptab.PrettyTable(['Year','Process','Events','Raw Yield','Bin Reweighted','BDT Reweighted' ] ) 
        sumY={'raw':0.0 ,'binned':0.0 ,'bdt' : 0.0,'Events':0}
        for ky in dataStore[yr]['bkg']:
            wraw   =np.sum(dataStore[yr]['bkg'][ky]['weight_v0']*dataStore[yr]['bkg'][ky]['lumi'])
            wbinned=np.sum(dataStore[yr]['bkg'][ky]['weight_binned'])
            wbdt=np.sum(dataStore[yr]['bkg'][ky]['weight_bdt'])
            outDict['yields'][ky]          =  {}
            
            scf=1.0
            if yr in miscScales:
                if 'bkg' in miscScales[yr]:
                    if ky in miscScales[yr]['bkg'] :
                        scf=miscScales[yr]['bkg'][ky]

            outDict['yields'][ky]['Events']=str(len( dataStore[yr]['bkg'][ky]['weight_bdt'] ))
            outDict['yields'][ky]['raw']   =str(np.round(wraw*scf,3))
            outDict['yields'][ky]['binned']=str(np.round(wbinned*scf,3))
            outDict['yields'][ky]['bdt']   =str(np.round(wbdt*scf,3))
            if ky !='bkg':
                sumY['Events'] += len(dataStore[yr]['bkg'][ky]['weight_bdt'] )
                sumY['raw'] += wraw
                sumY['binned'] += wbinned
                sumY['bdt'] += wbdt
            tabYield.add_row([yr , ky ,
                              outDict['yields'][ky]['Events'] , 
                              outDict['yields'][ky]['raw'] , 
                              outDict['yields'][ky]['binned'],
                              outDict['yields'][ky]['bdt']  ])
        tabYield.add_row([yr , 'sum'  ,
                          np.round(sumY['Events'],3) ,
                          np.round(sumY['raw'],3) ,
                          np.round(sumY['binned'],3) ,
                          np.round(sumY['bdt'],3) ]
                          )


        print(tabYieldData)
        print()
        print(tabYield)
        foutname=saveBase+'/'+'modelYields.json'
        with open(foutname, 'w') as f:
            json.dump(outDict,f,indent=4)

if __name__=='__main__':
    main( )

