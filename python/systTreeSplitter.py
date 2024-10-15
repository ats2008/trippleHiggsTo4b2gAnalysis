#!/usr/bin/env python

import ROOT
import json
import pickle as pkl
import os,argparse
import Util as utl


parser = argparse.ArgumentParser()
parser.add_argument("--cfg", help="cfg file" )

args=parser.parse_args()
cfg={}
with open(args.cfg,'r') as f:
    cfg=json.load(f)

yr=cfg['yr']
tag=cfg['tag']
hashStr=cfg['hash']
dest=cfg['destination']

doNotRedo=False
doNotRedo=True
if 'doNotRedo' in cfg:
    doNotRedo=cfg['doNotRedo']
doOnlyNominal=cfg['doOnlyNominal']
rdataFrames={yr : { tag : {}} }
treeFolderName=cfg['path']
scale=1.0
varsToGet=['weight']
if 'scale' in cfg:
    scale=cfg['scale']
print("Scale : ",scale)
idx=0
for fname in cfg['files']:
    print("Processing ",fname)   
    fname=fname.replace("root://se01.indiacms.res.in:1094/","/tmp/")
    print(fname)
    fileIn=ROOT.TFile.Open(fname)
    systDict={}
    baseTag=None
    treeDir=fileIn.Get(treeFolderName)
    treeDir.Print()
    for ky in treeDir.GetListOfKeys():
        name=str(ky.GetName())
        base     ='_'.join(name.split('_')[:-1])
        syst_tag =name.split('_')[-1]
        ofname=f"{dest}/{tag}_{yr}_{hashStr}_{idx}_{syst_tag}_scaled_split.root"
        if doNotRedo:
            if os.path.exists(ofname):
                continue
        if syst_tag=='0':
            syst_tag='nominal'
            baseTag=name
        if doOnlyNominal:
            if 'nominal' not in syst_tag:
                continue
        treeName=f'{treeFolderName}/{name}'
        print("Getting ",treeName)
        systDict[syst_tag]=  ROOT.RDataFrame(treeName, fname) 
        if systDict[syst_tag]:
            print("\tGot Sucessfully !",systDict[syst_tag].Count().GetValue())
    
    for sysTag in systDict:
        
        rdfIn     =  systDict[sysTag] 
        dataStore =  rdfIn.AsNumpy(varsToGet)
        scaledWeight = dataStore["weight"]*scale
        print(f"  Sum of input weights  : {rdfIn.Sum('weight').GetValue()} " )
        defnString   = 'auto to_eval = std::string("scaledWeight[") + std::to_string(rdfentry_) + "]"; '
        defnString  += 'return float(TPython::Eval(to_eval.c_str()));'
        rdfOut       =  rdfIn.Redefine("weight",defnString)
        print(f"  Sum of output weights : {rdfOut.Sum('weight').GetValue()} " )
        ofname=f"{dest}/{tag}_{yr}_{hashStr}_{idx}_{sysTag}_scaled_split.root"
        print(f"Writing out syst {sysTag} to the file {ofname}")
        treeOName='tagsDumper/bkg_13TeV_TrippleHTag_0'
        rdfOut.Snapshot(treeOName,ofname)

        #f=ROOT.TFile(ofname,"RECREATE")
        #f.mkdir(treeFolderName)
        #f.cd(treeFolderName)
        #treeOut=systDict[sysTag].CloneTree()
        #treeOut.SetName(baseTag)
        #treeOut.Write()
        #f.Purge()
        #f.Close()
        
    idx+=1


#if __name__=='__main__':
#    main()    
