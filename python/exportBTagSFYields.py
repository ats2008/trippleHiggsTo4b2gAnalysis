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

scf = 0.00017

hep.style.use("CMS")

def exportTreeForFit(final_dataset,treeNameToExport,export_filename): 
    export_file = urt.recreate(export_filename)
    print("\n\t  file exported  : ",export_filename)
    allBranchesToWrite={ ky : final_dataset[ky].dtype for ky in final_dataset}
    export_file.mktree( treeNameToExport ,allBranchesToWrite)
    export_file[treeNameToExport].extend(final_dataset)
    export_file.close()

def main():
    exportBranches=[
        'lumiID',
        'run',
        'event',
        'lumi',
        'year',
        'weight',
        'CMS_hgg_mass',
        'h1bb_mass',
        'h2bb_mass',
        'quadjet_0_deepJetScore',
        'quadjet_1_deepJetScore',
        'quadjet_2_deepJetScore',
        'quadjet_3_deepJetScore'
    ]




    data_blind='CMS_hgg_mass < 115 || CMS_hgg_mass > 135.0'
    sigRegionCut='CMS_hgg_mass >= 115 && CMS_hgg_mass <= 135.0'
    bkgRegionCut=data_blind
 
    parser = argparse.ArgumentParser()
    parser.add_argument("--cfg", help="cfg file" )
    parser.add_argument("--infoExport", help="Export selection output json  to a place !! ",default='./' )
    args = parser.parse_args()
    
    cfg={}
    with open(args.cfg,'r') as f:
        cfg=json.load(f)

    print(exportBranches)
    export=cfg['export']
    if export:
        sigRegionCut = ' 1==1 '
        print("Files will be exported")
    if 'extraVars' in cfg['extraVars'] :
        exportBranches+=cfg['extraVars']

    saveOutput = True
    cutsToApply=[]
    cutsKey=[]
    cutIdx=0
    txt=[]
    bkgRegionCut='CMS_hgg_mass < 115 || CMS_hgg_mass > 135.0'
    
    miscScales={}
    if 'scaleFactor' in cfg:
        with open(cfg['scaleFactor']) as f:
            miscScales=json.load(f)
            print("Loaded the scales json : ", cfg['scaleFactor'] ) #,json.dumps(miscScales,indent=4))
    if 'cuts' in cfg:
        with open(cfg['cuts'],'r') as f:
            txt_=f.readlines()
            for l in txt_:
                txt.append(l[:-1])
    
    for l in txt:
        if l[0]=='#':
            continue
        key='cut '+str(cutIdx)
        cut=','.join(l.split(',')[:-1])
        if ',' in l :
            key=l.split(',')[-1]
        cutsToApply.append(cut)
        cutsKey.append(key)
        cutIdx+=1
    print("Cuts Being applied : ")
    for cut in cutsToApply:
        print("\t -> ",cut)

    
    categorySelection={}
    catIdx=0
    if 'cats' in cfg:
        with open(cfg['cats'],'r') as f:
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
    prefixBase=cfg['destination']
    outDict={}
    blind=data_blind
    
    
    weight_var='weight_binScaledV1'  
    
    varToBinMap={}
    
    print()
    print("Outputs are stored in ",prefixBase)
    print()
    print("Weight Var set as : ",weight_var)
    
    # We read the tree from the file and create a RDataFrame, a class that
    # allows us to interact with the data contained in the tree.
    
    year=cfg['year']
    tag =cfg['tag']
    proc=cfg['proc']
    syst=cfg['syst']
    fileName=cfg['file']
    treeName=cfg['treeIn']
    treeNameOut=cfg['treeOut']
    print(f"Processing file     : {fileName}")
    print(f"        Input tree  : {treeName}")
    print(f"    Exporting tree  : {treeNameOut}")
    rdataFrames={
        year : {
            tag : {
              }
        }
    }
    try :
        rdataFrames[year][tag][proc]= ROOT.RDataFrame(treeName, fileName)    
        dd = rdataFrames[year][tag][proc].GetColumnNames()
    except:
        treeName = 'trees/ggHHH_125_13TeV'
        rdataFrames[year][tag][proc]= ROOT.RDataFrame(treeName, fileName)    
        dd = rdataFrames[year][tag][proc].GetColumnNames()

    for yr in rdataFrames:
        for tag in rdataFrames[yr]:
            for ky in rdataFrames[yr][tag]:
                print(f"Processing {yr} / {tag} / {ky}")
                allBranches=[ str(i) for i in rdataFrames[yr][tag][ky].GetColumnNames()]
                dFrame=rdataFrames[yr][tag][ky]
                if 'weight_binScaledV1' not in allBranches:
                    print("\tUpdating weight_binScaledV1 branch")
                    dFrame= dFrame.Define("weight_binScaledV1","weight_v0*lumi")
                if 'weight_cScaledV1' not in allBranches:
                    print("\tUpdating weight_cScaled branch")
                    dFrame= dFrame.Define("weight_cScaledV1","weight_v0*lumi")
                if 'weight_lumiScaled' not in allBranches:
                    print("\tUpdating weight_lumiScaled branch")
                    dFrame= dFrame.Define("weight_lumiScaled","weight_v0*lumi")
                if 'lumiID' not in allBranches:
                    exportBranches.remove('lumiID')
                rdataFrames[yr][tag][ky]=dFrame


    for yr  in rdataFrames:
        outDictGloabl={}
        for tag in rdataFrames[yr]:
            outDictGloabl[tag] = {}
            for proc in rdataFrames[yr][tag]:
                outDict={}
                #saveBase=prefixBase+'/'
                #cmd='mkdir -p '+saveBase ; print("[] $ ",cmd) ;  os.system(cmd)
                saveBase="./"
                print(f"Processing : {yr}/{tag}/{proc} | Systemetic : {syst} ")
                scf=1.0
                #print(json.dumps(miscScales,indent=4))
                if yr in miscScales:
                    print("GOT yr for scaling ",yr)
                    if tag in miscScales[yr]:
                        print("GOT Tag for scaling ",tag)
                        for prc in miscScales[yr][tag] :
                            if proc in prc:
                                print("GOT proc for scaling ",proc)
                                scf=miscScales[yr][tag][prc]
                                print(f"#INFO scaling weights for {yr}/{tag}/{proc} [ {prc} ]| Systemetic : {syst} by {scf:.5f}")
                NEvts=rdataFrames[yr][tag][proc].Count().GetValue()
                weightPrev    = rdataFrames[yr][tag][proc].Sum(weight_var).GetValue()*scf
                weightCurrent = rdataFrames[yr][tag][proc].Sum(weight_var).GetValue()*scf
                effStep=1.0
                effCumulative=1.0
                outDict['Orig']={'cut' : '' ,'nEvts':NEvts,'weight':weightCurrent , 'effStep' : effStep,"effCumulative":effCumulative}
                
                rdataFrames[yr][tag][proc]= rdataFrames[yr][tag][proc].Filter(sigRegionCut,"SR cut")
                weightCurrent=rdataFrames[yr][tag][proc].Sum(weight_var).GetValue()*scf
                NEvts=rdataFrames[yr][tag][proc].Count().GetValue()
                effStep=weightCurrent/(1e-14+weightPrev)
                effCumulative*=effStep
                outDict['SR']={'cut' : sigRegionCut  ,'nEvts':NEvts,'weight':weightCurrent , 'effStep' : effStep,"effCumulative":effCumulative}
                weightPrev=weightCurrent

                print("\tProcessing signal Cuts ")
                for cut,cutKey in zip(cutsToApply,cutsKey):
                    rdataFrames[yr][tag][proc] = rdataFrames[yr][tag][proc].Filter(cut,cutKey)
                    weightCurrent=rdataFrames[yr][tag][proc].Sum(weight_var).GetValue()*scf
                    NEvts=rdataFrames[yr][tag][proc].Count().GetValue()
                    effStep=weightCurrent/(1e-14+weightPrev)
                    effCumulative*=effStep
                    outDict[cutKey]={'cut' : cut ,'nEvts':NEvts,'weight':weightCurrent , 'effStep' : effStep,"effCumulative":effCumulative}
                    weightPrev=weightCurrent
                
                catDset=rdataFrames[yr][tag][proc]
                outDict['category']={}
                weightBase  = weightPrev
                cumuEffBase = effCumulative
                print("\tProcessing signal Categories ")
                for i in range(len(categorySelection)):
                    cut=categorySelection[i]['cut']
                    cutKey=categorySelection[i]['name']
                    catEvts = catDset.Filter(cut,cutKey)
                    weightCurrent=catEvts.Sum(weight_var).GetValue()*scf
                    NEvts=catEvts.Count().GetValue()
                    effStep=weightCurrent/(1e-14+weightBase)
                    effCumulative= cumuEffBase*effStep
                    outDict['category'][cutKey]={'cat' : cut ,'nEvts':NEvts,'weight':weightCurrent , 'effStep' : effStep,"effCumulative":effCumulative}
                    print("\t\t Making category ",cutKey)
                    
                    # ## Exporting the augmented dataset
                    if cfg['export']:       
                        final_dataset = catEvts.AsNumpy(columns=exportBranches +["dZ"]+[weight_var]) 
                        final_dataset["weight"]= scf * final_dataset[weight_var] /final_dataset["lumi"]
                        export_filename=saveBase+f'/{tag}_{proc}_{yr}_{syst}_{cutKey}_M125_13TeV.root' 
                        treeNameToExport = treeNameOut.replace( '13TeV',f'13TeV_{cutKey}'  )
                        exportTreeForFit(final_dataset,treeNameToExport,export_filename)

                    catDset = catDset.Filter("! ( "+cut+" ) ","NOT "+cutKey)
       
                cutKey='CATX'
                if cfg['exportX']:       
                    final_dataset = catDset.AsNumpy(columns=exportBranches +["dZ"]) 
                    final_dataset["weight"]= scf * final_dataset["weight"] /final_dataset["lumi"]
                    export_filename=saveBase+f'/{tag}_{proc}_{yr}_{syst}_{cutKey}_M125_13TeV.root' 
                    treeNameToExport = treeNameOut.replace( '13TeV',f'13TeV_{cutKey}'  )
                    exportTreeForFit(final_dataset,treeNameToExport,export_filename)

 
                allCutsReport = rdataFrames[yr][tag][proc].Report()
                allCutsReport.Print()          
                print()
                print("- -"*20)
                print()
                
                outDictGloabl[tag][proc]=outDict      

        yieldPerCat={}           
        for ky in outDictGloabl:
            outDict=outDictGloabl[ky]

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

                            if cat not in yieldPerCat:
                                yieldPerCat[cat]= ptab.PrettyTable(['Process','Events','Yield','Cumu. Efficiency'] )
                            row=[ sample                                        , str(outDict[sample]['category'][cat]['nEvts']),
                                  str(np.round(outDict[sample]['category'][cat]['weight'] ,3 )), str(np.round(outDict[sample]['category'][cat]['effCumulative'],5)) ]
                            
                            yieldPerCat[cat].add_row(row)
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
                
                if cfg['cutFlow']:
                    print()
                    print(tabCuts)
                    print()
               
                    print()
                    print(tabCutNames)
                    print()

                print()
                print(tabCats)
                print()
        
        for cat in yieldPerCat:
            print()
            print("CAT : ",cat)
            print()
            print( yieldPerCat[cat]  )
            print()
           
           #print()
           #print(tabCatNames)
           #print()
        foutname=f'{cfg["destination"]}/Global_Yields_{yr}.json'
        with open(foutname,'w') as f:
            print(f"Exporting resuts to : {foutname}")
            exportDict={'selections' : outDictGloabl , 'cfg' : cfg}
            json.dump(exportDict,f,indent=4)
        
       
if __name__=='__main__':
    main( )

