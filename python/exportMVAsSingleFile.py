#import matplotlib.pyplot as plt
#import mplhep as hep
import uproot as urt
import ROOT
import json,sys,os,argparse
import numpy as np
#import Util as utl
import pickle,itertools
#import prettytable as ptab
#import scaleFactorUtil as scl

def getScoreForArray(model,features,allVarDict,classIdx=[]):
    x=np.stack([allVarDict[var] for var in features ]).T
    scores = model.predict_proba(x)
    
    if classIdx==[]:
        return scores
    else:
        return [ scores[i] for i in classIdx]




def exportTreeForFit(final_dataset,treeNameToExport,export_filename): 
    export_file = urt.recreate(export_filename)
    print("\n\t  file exported  : ",export_filename)
    allBranchesToWrite={ ky : final_dataset[ky].dtype for ky in final_dataset}
    export_file.mktree( treeNameToExport ,allBranchesToWrite)
    export_file[treeNameToExport].extend(final_dataset)
    export_file.close()

btagBranches=["btag_scale","btag_scale_down_hf","btag_scale_down_hfstats1","btag_scale_down_hfstats2","btag_scale_down_cferr1","btag_scale_down_cferr2",
                           "btag_scale_down_jes","btag_scale_down_lf","btag_scale_down_lfstats1","btag_scale_down_lfstats2",
                           "btag_scale_up_hf","btag_scale_up_hfstats1","btag_scale_up_hfstats2","btag_scale_up_jes","btag_scale_up_cferr1","btag_scale_up_cferr2",
                         "btag_scale_up_lf","btag_scale_up_lfstats1","btag_scale_up_lfstats2"]

def main():
    exportBranches=[
        'lumiID',
        'run',
        'event',
        'lumi',
        'year',
        'weight',
        "weight_v0",
        'weight_cScaledV1',
        'CMS_hgg_mass',
        'h1bb_mass',
        'h2bb_mass',
        'leadingPhoton_pt',
        'subleadingPhoton_pt',
        'leadingPhoton_eta',
        'subleadingPhoton_eta',
        'leadingPhoton_phi',
        'subleadingPhoton_phi',
        'quadjet_0_deepJetScore',
        'quadjet_1_deepJetScore',
        'quadjet_2_deepJetScore',
        'quadjet_3_deepJetScore',
        'quadjet_0_flavour',
        'quadjet_1_flavour',
        'quadjet_2_flavour',
        'quadjet_3_flavour',
	    'eta_h1leadJ',
	    'eta_h1subleadJ',
	    'eta_h2leadJ',
	    'eta_h2subleadJ',
	    'phi_h1leadJ',
	    'phi_h1subleadJ',
	    'phi_h2leadJ',
	    'phi_h2subleadJ',
	    'pT_h1leadJ',
	    'pT_h1subleadJ',
	    'pT_h2leadJ',
	    'pT_h2subleadJ',
        'quadjet_0_mass',
        'quadjet_1_mass',
        'quadjet_2_mass',
        'quadjet_3_mass',
    ]
    preselBranches=['quadjet_0_isTight','quadjet_1_isTight','quadjet_2_isTight','quadjet_3_isTight',
                    'quadjet_0_isPUTight','quadjet_1_isPULoose','quadjet_2_isPUTight','quadjet_3_isPULoose',
                    'customLeadingPhotonIDMVA','customSubLeadingPhotonIDMVA','mva_supressionBDT_v16','mva_vsPeakingBDT_vsTTH_v27p1_0']
    
    exportBranches=exportBranches+preselBranches

    data_blind='CMS_hgg_mass < 115 || CMS_hgg_mass > 135.0'
    sigRegionCut='CMS_hgg_mass >= 115 && CMS_hgg_mass <= 135.0'
    bkgRegionCut=data_blind
 
    parser = argparse.ArgumentParser()
    parser.add_argument("--cfg", help="cfg file" )
    parser.add_argument("--ignoreBTag", help="do not check btag",action='store_true' )
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
    
    prefixBase=cfg['destination']
    outDict={}
    
    weight_var='weight_binScaledV1'  
    varToBinMap={}
    
    print()
    #aprint("Outputs are stored in ",prefixBase)
    print()
    print("Weight Var set as : ",weight_var)
    
    # We read the tree from the file and create a RDataFrame, a class that
    # allows us to interact with the data contained in the tree.
    
    year=cfg['year']
    tag =cfg['tag']
    proc=cfg['proc']
    syst=cfg['syst']
    if 'morphedData' == cfg['proc']:
        args.ignoreBTag=True
    if  not args.ignoreBTag:
        if syst=='nominal':
            print("Adding btag branches")
            exportBranches+=btagBranches
        elif 'JEC' in syst:
            print("Adding btag branches")
            exportBranches+=btagBranches

        
    fileName=cfg['file']
    treeName=cfg['treeIn']
    treeNameOut=cfg['treeIn']
    print(f"Processing file     : {fileName}")
    print(f"        Input tree  : {treeName}")
    print(f"    Exporting tree  : {treeNameOut}")
    
    mvaTags  =cfg['mvaTags']
    mvaFiles =cfg['mvaFiles']
    supressionModel={}   
    varsForMVA=[]
    for modelTag,fName in zip(mvaTags,mvaFiles):  
        print(f"Loading BDT {modelTag} from {fName}")
        with open(fName, 'rb') as f:
            modelDict = pickle.load(f)
            supressionModel[modelTag]=modelDict
            varsForMVA+=modelDict['features']

    varsForMVA=np.unique(varsForMVA)
    for v in varsForMVA:
        if v not in exportBranches:
            exportBranches.append(v)
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
        try :
            treeName = 'trees/bkg_13TeV_TrippleHTag_0'
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
                dFrame=rdataFrames[yr][tag][ky]
                dset=dFrame.AsNumpy(columns=[i  for i in ["weight_v0","lumi"]])
                print("  Dset  : ",len(dset['weight_v0'])," events  | Sum weights : ",np.sum(dset['weight_v0']) )
                allBranches=[ str(i) for i in rdataFrames[yr][tag][ky].GetColumnNames()]
                dFrame= dFrame.Define("weight_raw","weight_v0")
                #if False and (not args.ignoreBTag):
                if not args.ignoreBTag:
                    dFrame=dFrame.Redefine("weight_v0","weight_v0*btag_scale")
                if 'log_leadingPhotonSigOverE' not in allBranches:
                    print("\tUpdating log_leadingPhotonSigOverE branch")
                    dFrame= dFrame.Define("log_leadingPhotonSigOverE","log(leadingPhotonSigOverE)")
                if 'minPhotonID' not in allBranches:
                    print("\tUpdating triHigs_ptOverMass branch")
                    dFrame= dFrame.Define("triHigs_ptOverMass","trihiggs_pt/trihiggs_mass")
                if 'minPhotonID' not in allBranches:
                    print("\tUpdating minPhotonID branch")
                    dFrame= dFrame.Define("minPhotonID","min(customSubLeadingPhotonIDMVA,customLeadingPhotonIDMVA)")
                if 'maxPhotonID' not in allBranches:
                    print("\tUpdating maxPhotonID branch")
                    dFrame= dFrame.Define("maxPhotonID","max(customSubLeadingPhotonIDMVA,customLeadingPhotonIDMVA)")
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
                for br in btagBranches:
                    if (br not in allBranches) and (br in exportBranches):
                        exportBranches.remove(br)
                for br in allBranches:
                    if ('btag' in br) and (br not in exportBranches):
                        exportBranches.append(br)
                rdataFrames[yr][tag][ky]=dFrame


    for yr  in rdataFrames:
        outDictGloabl={}
        for tag in rdataFrames[yr]:
            outDictGloabl[tag] = {}
            for proc in rdataFrames[yr][tag]:
                outDict={}
                saveBase="./"
                print(f"Processing : {yr}/{tag}/{proc} | Systemetic : {syst} ")
                final_dataset = dFrame.AsNumpy(columns=exportBranches +["dZ"]+[weight_var])
                dset=dFrame.AsNumpy(columns=[i  for i in varsForMVA])
                mvaExportNames=[]
                for model_tag in supressionModel:
                    model=supressionModel[model_tag]
                    varsForClassifier=model['features']
                    clf=model['classifier']
                    scores=getScoreForArray(clf,varsForClassifier,dset)
                    for ii in range(len(scores[0])):
                        name=f'{model_tag}_{ii}'
                        final_dataset[name]=scores[:,ii]
                        mvaExportNames+=[name]
                        print(f"  > {tag}/{ky}/{name} : {np.average(scores[:,ii])}")
                export_filename=saveBase+f'/{tag}_{proc}_{yr}_{syst}_M125_13TeV.root' 
                treeNameToExport = treeNameOut
                exportTreeForFit(final_dataset,treeNameToExport,export_filename)

                
                print()
                print("- -"*20)
                print()
                
                outDictGloabl[tag][proc]=outDict      
       
       
if __name__=='__main__':
    main( )

