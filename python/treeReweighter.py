#!/usr/bin/env python

import ROOT
import json
import pickle as pkl
import os,argparse
import Util as utl

import numpy as np
import awkward as ak
import vector
vector.register_awkward()

def getMomentumArrayFromObjects(raw_data,obj_tagList=["Photon"],extraColums=['mass']):
    allColums=['pt','eta','phi']+extraColums
    objects = ak.zip({
        obj_tag: ak.zip({
             ky : raw_data[f"{obj_tag}_{ky}"] for ky in allColums
        }, with_name="Momentum4D")  for obj_tag in obj_tagList
    }, depth_limit=1)
    return objects

def getObjectAsArray(raw_data,object_tag="Photon",colums=["pt","eta","phi","mass"],nameTag="Momentum4D"):
#     print(object_tag)
#     print(colums)
    objectArr=[ 
                ak.zip({var : raw_data[f"{object_tag}{i}_{var}"]  for var in colums},with_name=nameTag)
                 for i in range(1,6+1)
            ]
    objectArr=[ ak.unflatten(objArr,counts=1,axis=-1) for objArr in objectArr]
    objectArr= ak.concatenate(objectArr,axis=1)
    return objectArr    


    
parser = argparse.ArgumentParser()
parser.add_argument("--cfg", help="cfg file" )


args=parser.parse_args()
cfg={}
with open(args.cfg,'r') as f:
    cfg=json.load(f)

dest=cfg['destination']
treeName=cfg['tree']
scale_value=cfg['scale']
varsToGet=['weight']
print(" Scale : ", scale_value)

idx=0
for fname in cfg['files']:
    print(f"Processing {treeName} in {fname}")
    rdfIn= None
    rdfIn =   ROOT.RDataFrame(treeName, fname)
    if not rdfIn:
        print(f"{treeName} failed defaulting to  : tagsDumper/trees/bkg_13TeV_TrippleHTag_0 ")
        rdfIn =   ROOT.RDataFrame("tagsDumper/trees/bkg_13TeV_TrippleHTag_0", fname)
    
    dataStore= rdfIn.AsNumpy(varsToGet)
    scaledWeight = dataStore["weight"]*scale_value
    print(f"  Sum of input weights  : {rdfIn.Sum('weight').GetValue()} " )
    defnString   = 'auto to_eval = std::string("scaledWeight[") + std::to_string(rdfentry_) + "]"; '
    defnString  += 'return float(TPython::Eval(to_eval.c_str()));'
    rdfOut       =  rdfIn.Redefine("weight",defnString)
    print(f"  Sum of output weights : {rdfOut.Sum('weight').GetValue()} " )

    outFname = dest +'/'+ fname.split('/')[-1].replace(".root",'_scaled.root')
    f=ROOT.TFile(outFname,"RECREATE")
    f.mkdir('tagsDumper/trees')
    print(f"Writing Outputs to {treeName} in {outFname} [ {rdfOut.Count().GetValue()} events]")
    rdfOut.Snapshot(treeName,outFname)
    print()
    f.Close()

#if __name__=='__main__':
#    main()    
