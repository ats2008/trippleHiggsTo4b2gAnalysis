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


hep.style.use("CMS")
data_blind='CMS_hgg_mass < 115 || CMS_hgg_mass > 135.0'
sigRegionCut='CMS_hgg_mass >= 115 && CMS_hgg_mass <= 135.0'
bkgRegionCut=data_blind

def main():
 
    parser = argparse.ArgumentParser()
    parser.add_argument('-v',"--version", help="Version of the specific derivation ",default='')
    parser.add_argument('-i',"--inputFile", help="Input File",default=None)
    parser.add_argument('-y',"--year", help="Year", default='run2')
    parser.add_argument('-t',"--tag" , help="Tag" , default='')
    parser.add_argument("-o","--dest", help="destination To Use", default='workarea/results/plots/tmp/' )
    parser.add_argument("-c","--cuts", help="list of cuts to be applied", default=None )
    parser.add_argument(     "--cats", help="Category Selction", default=None )
    parser.add_argument(     "--cutstr", help="additional cutsrting to be applied", default=None )
    parser.add_argument(     "--doSR", help="Do Signal Region",default=False,action='store_true')
    parser.add_argument(     "--cutFlow", help="Print cut flow",default=False,action='store_true')
    parser.add_argument(     "--scaleFactor", help="additional scales to be multiplied  as json { year { bkg : { deset : <SF> } }}", default=None )
    args = parser.parse_args()
    version = args.version
    doIronTrnsformation=False
    inputFile  = args.inputFile
    saveOutput = True
    cutsToApply=[]
    cutsKey=[]
    cutIdx=0
    txt=[]
    bkgRegionCut=data_blind
    if args.doSR:
        bkgRegionCut=sigRegionCut
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
    print("Cuts Being applied : ")
    for cut in cutsToApply:
        print("\t -> ",cut)
    miscScales={}
    if args.scaleFactor:
        with open(args.scaleFactor) as f:
            miscScales=json.load(f)
            print("Loaded the scales json : ",json.dumps(miscScales,indent=4))

    
    categorySelection={}
    catIdx=0
    if args.cats:
        with open(args.cats,'r') as f:
            txt=f.readlines()
            for l in txt:
                items=l[:-1].split(',')
                categorySelection[catIdx]={
                        'cut' :items[0],
                        'name':items[1].strip()
                    }
                catIdx+=1
    print()
    print("Categories being considered : ")
    for i in range(len(categorySelection)):
        print("\t",i," -> ",categorySelection[i]['name'],"   --> ", categorySelection[i]['cut']  )
    print()
    print()
    prefixBase=args.dest
    outDict={}
    blind=data_blind
    
    print("Processing file list : ",args.inputFile)
    if inputFile:
        saveOutput=False
    fileDict={}
    with open(args.inputFile) as f:
        fileDict=json.load(f)
    
   
    yearsToProcess_all=['2018','2017','2016PreVFP','2016PostVFP','run2','extras']
    yearsToProcess=yearsToProcess_all
    
    if 'all' not in args.year:
        yearsToProcess=[]
        for yr in args.year.split(","):
            if yr not in yearsToProcess_all:
                print(yr," not in catalogue. Skipping !! ")
                continue
            yearsToProcess.append(yr)

    bkgToProcess=[  'bkg',"*"
            #        'ggBox1Bjet',
            #        'ggBox2Bjet', 
            #        'ggBox', 
            #        'gJet20To40',
            #        'gJet40ToInf'
                   ]
    varToBinMap={}
    
    print()
    print("Outputs are stored in ",prefixBase)
    print()
    # We read the tree from the file and create a RDataFrame, a class that
    # allows us to interact with the data contained in the tree.
    rdataFrames={}
    for yr  in yearsToProcess:
        saveBase=prefixBase+'/'+yr
        os.system('mkdir -p '+saveBase)
        outDict={}
        rdataFrames[yr]={'sig':{},'bkg':{},'data':{}}
        try:
            fileName=fileDict[yr]['sig']['ggHHH']
        except:
            print(yr)
        treeName = "trees/ggHHH_125_13TeV"
        
        rdataFrames[yr]['sig']['ggHHH']= ROOT.RDataFrame(treeName, fileName)
        
        outDict['sig'] = {}
        NEvts=rdataFrames[yr]['sig']['ggHHH'].Count().GetValue()
        scf=1.0
        if yr in miscScales:
            if 'sig' in miscScales[yr]:
                if 'ggHHH' in miscScales[yr]['sig'] :
                    scf=miscScales[yr]['sig']['ggHHH']
        weightPrev   =rdataFrames[yr]['sig']['ggHHH'].Sum('weight_bdt').GetValue()*scf
        weightCurrent=rdataFrames[yr]['sig']['ggHHH'].Sum('weight_bdt').GetValue()*scf
        effStep=1.0
        effCumulative=1.0
        outDict['sig']['Orig']={'cut' : '' ,'nEvts':NEvts,'weight':weightCurrent , 'effStep' : effStep,"effCumulative":effCumulative}
        
        rdataFrames[yr]['sig']['ggHHH']= rdataFrames[yr]['sig']['ggHHH'].Filter(sigRegionCut,"SR cut")
        weightCurrent=rdataFrames[yr]['sig']['ggHHH'].Sum('weight_bdt').GetValue()*scf
        NEvts=rdataFrames[yr]['sig']['ggHHH'].Count().GetValue()
        effStep=weightCurrent/weightPrev
        effCumulative*=effStep
        outDict['sig']['SR']={'cut' : sigRegionCut  ,'nEvts':NEvts,'weight':weightCurrent , 'effStep' : effStep,"effCumulative":effCumulative}
        weightPrev=weightCurrent

        for cut,cutKey in zip(cutsToApply,cutsKey):
            rdataFrames[yr]['sig']['ggHHH'] = rdataFrames[yr]['sig']['ggHHH'].Filter(cut,cutKey)
            weightCurrent=rdataFrames[yr]['sig']['ggHHH'].Sum('weight_bdt').GetValue()*scf
            NEvts=rdataFrames[yr]['sig']['ggHHH'].Count().GetValue()
            effStep=weightCurrent/weightPrev
            effCumulative*=effStep
            outDict['sig'][cutKey]={'cut' : cut ,'nEvts':NEvts,'weight':weightCurrent , 'effStep' : effStep,"effCumulative":effCumulative}
            weightPrev=weightCurrent
        
        catDset=rdataFrames[yr]['sig']['ggHHH']
        outDict['sig']['category']={}
        for i in range(len(categorySelection)):
            cut=categorySelection[i]['cut']
            cutKey=categorySelection[i]['name']
            catEvts = catDset.Filter(cut,cutKey)
            weightCurrent=catEvts.Sum('weight_bdt').GetValue()*scf
            NEvts=catEvts.Count().GetValue()
            effStep=weightCurrent/weightPrev
            effCumulativeCat=effStep*effCumulative
            outDict['sig']['category'][cutKey]={'cat' : cut ,'nEvts':NEvts,'weight':weightCurrent , 'effStep' : effStep,"effCumulative":effCumulativeCat}
            
            #weightPrev=weightCurrent
            catDset = catDset.Filter("! ( "+cut+" ) ","NOT "+cutKey)

        print("Registering datset : ggHHH , ",yr," with tree",treeName)
        print("\t\t ",fileName)   
        allCutsReport = rdataFrames[yr]['sig']['ggHHH'].Report()
        allCutsReport.Print()          
        print()
        print("- -"*20)
        print()
        if 'data' in fileDict[yr]:
            ky=list(fileDict[yr]['data'].keys())[0]
            fileName=fileDict[yr]['data'][ky]
            treeName = "trees/Data_13TeV_TrippleHTag_0"
            
            rdataFrames[yr]['data']['data']= ROOT.RDataFrame(treeName, fileName).Filter(data_blind)
            outDict['data'] = {}
            NEvts=rdataFrames[yr]['data']['data'].Count().GetValue()
            weightPrev   =rdataFrames[yr]['data']['data'].Sum('weight_bdt').GetValue()
            weightCurrent=rdataFrames[yr]['data']['data'].Sum('weight_bdt').GetValue()
            effStep=1.0
            effCumulative=1.0
            outDict['data']['Orig']={'cut' : '' ,'nEvts':NEvts,'weight':weightCurrent , 'effStep' : effStep,"effCumulative":effCumulative}
  
            rdataFrames[yr]['data']['data']= rdataFrames[yr]['data']['data'].Filter(blind)
            weightCurrent=rdataFrames[yr]['data']['data'].Sum('weight_bdt').GetValue()
            NEvts=rdataFrames[yr]['data']['data'].Count().GetValue()
            effStep=weightCurrent/weightPrev
            effCumulative*=effStep
            outDict['data']['Sideband']={'cut' : blind ,'nEvts':NEvts,'weight':weightCurrent , 'effStep' : effStep,"effCumulative":effCumulative}
            weightPrev=weightCurrent
            
            for cut,cutKey in zip(cutsToApply,cutsKey):
                rdataFrames[yr]['data']['data'] = rdataFrames[yr]['data']['data'].Filter(cut,cutKey)
                weightCurrent=rdataFrames[yr]['data']['data'].Sum('weight_bdt').GetValue()
                NEvts=rdataFrames[yr]['data']['data'].Count().GetValue()
                effStep=weightCurrent/(1e-9 + weightPrev)
                effCumulative*=effStep
                outDict['data'][cutKey]={'cut' : cut ,'nEvts':NEvts,'weight':weightCurrent , 'effStep' : effStep,"effCumulative":effCumulative}
                #weightPrev=weightCurrent
            
            catDset=rdataFrames[yr]['data']['data']
            outDict['data']['category']={}
            for i in range(len(categorySelection)):
                cut=categorySelection[i]['cut']
                cutKey=categorySelection[i]['name']
                catEvts = catDset.Filter(cut,cutKey)
                weightCurrent=catEvts.Sum('weight_bdt').GetValue()
                NEvts=catEvts.Count().GetValue()
                effStep=weightCurrent/weightPrev
                effCumulativeCat=effCumulative*effStep
                outDict['data']['category'][cutKey]={'cat' : cut ,'nEvts':NEvts,'weight':weightCurrent , 'effStep' : effStep,"effCumulative":effCumulativeCat}
                
                weightPrev=weightCurrent
                catDset = catDset.Filter("! ( "+cut+" ) ","NOT "+cutKey)
            
            print("Registering datset : data , ",yr," with tree",treeName)
            print("\t\t ",fileName)   
            allCutsReport = rdataFrames[yr]['data']['data'].Report()
            allCutsReport.Print()          
            print()
            print("- -"*20)
        print()
        rdataFrames[yr]['bkg']={}
        treeName = "trees/bkg_13TeV_TrippleHTag_0"
        for bkg in fileDict[yr]['bkg']:
            if "*" not in bkgToProcess:
                if bkg not in bkgToProcess:
                    print()
                    print("FILE NOT FOUNG FOR : ",yr," background ",bkg)
                    print()
                    continue
            fileName = fileDict[yr]['bkg'][bkg]
            
            rdataFrames[yr]['bkg'][bkg]=ROOT.RDataFrame(treeName, fileName)
            outDict[bkg] = {}
            NEvts=rdataFrames[yr]['bkg'][bkg].Count().GetValue()

            scf=1.0
            if yr in miscScales:
                if 'bkg' in miscScales[yr]:
                    if bkg in miscScales[yr]['bkg'] :
                        scf=miscScales[yr]['bkg'][bkg]

            weightPrev   =rdataFrames[yr]['bkg'][bkg].Sum('weight_bdt').GetValue()*scf
            weightCurrent=rdataFrames[yr]['bkg'][bkg].Sum('weight_bdt').GetValue()*scf
            effStep=1.0
            effCumulative=1.0
            outDict[bkg]['Orig']={'cut' : '' ,'nEvts':NEvts,'weight':weightCurrent , 'effStep' : effStep,"effCumulative":effCumulative}
  
            rdataFrames[yr]['bkg'][bkg]=rdataFrames[yr]['bkg'][bkg].Filter(bkgRegionCut,"SR Cut")

            weightCurrent=rdataFrames[yr]['bkg'][bkg].Sum('weight_bdt').GetValue()*scf
            NEvts=rdataFrames[yr]['bkg'][bkg].Count().GetValue()
            effStep=weightCurrent/weightPrev
            effCumulative*=effStep
            outDict[bkg]['region']={'cut' :bkgRegionCut ,'nEvts':NEvts,'weight':weightCurrent , 'effStep' : effStep,"effCumulative":effCumulative}
            weightPrev=weightCurrent
        
            for cut,cutKey in zip(cutsToApply,cutsKey):
                rdataFrames[yr]['bkg'][bkg] = rdataFrames[yr]['bkg'][bkg].Filter(cut,cutKey)
                weightCurrent=rdataFrames[yr]['bkg'][bkg].Sum('weight_bdt').GetValue()*scf
                NEvts=rdataFrames[yr]['bkg'][bkg].Count().GetValue()
                effStep=weightCurrent/(1e-9 + weightPrev)
                effCumulative*=effStep
                outDict[bkg][cutKey]={'cut' : cut ,'nEvts':NEvts,'weight':weightCurrent , 'effStep' : effStep,"effCumulative":effCumulative}
                weightPrev=weightCurrent
 
            
            catDset=rdataFrames[yr]['bkg'][bkg]
            outDict[bkg]['category']={}
            for i in range(len(categorySelection)):
                cut=categorySelection[i]['cut']
                cutKey=categorySelection[i]['name']
                catEvts = catDset.Filter(cut,cutKey)
                weightCurrent=catEvts.Sum('weight_bdt').GetValue()*scf
                NEvts=catEvts.Count().GetValue()
                effStep=weightCurrent/(1e-9+weightPrev)
                effCumulativeCat=effCumulativeCat*effStep
                outDict[bkg]['category'][cutKey]={'cat' : cut ,'nEvts':NEvts,'weight':weightCurrent , 'effStep' : effStep,"effCumulative":effCumulativeCat}
                
                #weightPrev=weightCurrent
                catDset = catDset.Filter("! ( "+cut+" ) ","NOT "+cutKey)
    
    

            print("Registering datset : ",bkg," , ",yr," with tree",treeName)
            print("\t\t ",fileName)   
            
            allCutsReport = rdataFrames[yr]['bkg'][bkg].Report()
            allCutsReport.Print()          
            print()
            print("- -"*20)
            print()
        
        print("Output saved at : ",saveBase)
        foutname=saveBase+'/'+args.tag+'catYields.json'
        
        with open(foutname, 'w') as f:
            json.dump(outDict,f,indent=4)
        
        print("\n\n")

        for sample in outDict:
           
            print()
            print('= = ='*20)
            print()
            print("Sample  : ",sample)
            
            tabCuts    =ptab.PrettyTable(['Cut','Events','Yield','Cut Efficiency','Cumu. Efficiency'] ) 
            tabCutNames=ptab.PrettyTable(['Cut','Discription'] ) 
            
            tabCats    =ptab.PrettyTable(['Category','Events','Yield','Cut Efficiency','Cumu. Efficiency'] ) 
            tabCatNames=ptab.PrettyTable(['Category','Discription'] ) 
            
            for cut in outDict[sample]:
                if cut=='category':
                    
                    for cat in outDict[sample]['category']:
                        row=[
                                cat,
                                str(outDict[sample]['category'][cat]['nEvts']),
                                str(np.round(outDict[sample]['category'][cat]['weight'],3)),
                                str(np.round(outDict[sample]['category'][cat]['effStep'],3)),
                                str(np.round(outDict[sample]['category'][cat]['effCumulative'],3)),
                            ]
                        tabCats.add_row(row)

                        row=[ cat , outDict[sample]['category'][cat]['cat'] ]
                        tabCatNames.add_row(row)

                    continue

                row=[
                        cut,
                        str(outDict[sample][cut]['nEvts']),
                        str(np.round(outDict[sample][cut]['weight'],3)),
                        str(np.round(outDict[sample][cut]['effStep'],3)),
                        str(np.round(outDict[sample][cut]['effCumulative'],3)),
                    ]
                tabCuts.add_row(row)

                row=[ cut , outDict[sample][cut]['cut'] ]
                tabCutNames.add_row(row)
            
            print()
            print(tabCats)
            print()
            
            if args.cutFlow:
                print()
                print(tabCatNames)
                print()
           
                print()
                print(tabCuts)
                print()
           
                print()
                print(tabCutNames)
                print()

if __name__=='__main__':
    main( )

