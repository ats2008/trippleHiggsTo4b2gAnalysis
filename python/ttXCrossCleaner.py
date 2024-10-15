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


if True:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cfg", help="cfg file" )
    parser.add_argument("--nMax", help="Max number of allowed gen-photons ",default=0 )

    
    args=parser.parse_args()
    cfg={}
    with open(args.cfg,'r') as f:
        cfg=json.load(f)

    dest=cfg['destination']
    if os.path.exists(dest):
        print("Making the destination directory ",dest)
        os.system('mkdir -p '+dest)

    treeName=cfg['tree']
    nGPhotonMax = 0
    if 'nMax' in cfg:
        nGPhotonMax = cfg['nMax']
    if args.nMax > 0:
        nGPhotonMax = args.nMax

    varsToGet=[]
    tmp=[
         'gen_promptG@@IDX_e','gen_promptG@@IDX_eta','gen_promptG@@IDX_phi','gen_promptG@@IDX_pt',
         'gen_promptG@@IDX_isFHPFS','gen_promptG@@IDX_isHard','gen_promptG@@IDX_grandmother',
         'gen_promptG@@IDX_mother', 'gen_promptG@@IDX_numberOfDaughters','gen_promptG@@IDX_pdgId'
        ]
    for i in range(1,7):
        for ky in tmp:
            varsToGet+=[ky.replace('@@IDX',f'{i}')]

    idx=0
    for fname in cfg['files']:
        print(f"Processing {treeName} in {fname}")
        rdfIn =   ROOT.RDataFrame(treeName, fname)
       
        dataStore= rdfIn.AsNumpy(varsToGet)
        for i in range(1,7):
            dataStore[f'gen_promptG{i}_mass']=np.zeros_like(dataStore[varsToGet[0]])
        eCols=['isFHPFS', 'isHard', 'mother','grandmother', 'numberOfDaughters', 'pdgId']
        
        genPhotons=getObjectAsArray(dataStore, object_tag='gen_promptG',colums=eCols+["pt","eta","phi","mass","e"])
        genPhotons=ak.drop_none(ak.mask(genPhotons, genPhotons.pt > -1 ))
        
        # Making the mask for the photons our interest , ie the ones that can overlap from tt+nG samples
        genPhoMask =  np.abs(genPhotons.mother)==6 
        # Light Quarks not from top Deacay
        l_mask =  np.abs(genPhotons.mother)==1 
        l_mask = l_mask | ( np.abs(genPhotons.mother)==2 )
        l_mask = l_mask | ( np.abs(genPhotons.mother)==3 )
        l_mask = l_mask | ( np.abs(genPhotons.mother)==4 )
        l_mask = l_mask | ( np.abs(genPhotons.mother)==5 )
        l_mask = l_mask & ( np.abs(genPhotons.grandmother)!=24 )
        
        genPhoMask = genPhoMask | l_mask
        genPhotonSelected = genPhotons[genPhoMask]
        
        
        numPhotons = ak.num(genPhotonSelected).to_numpy()

        print(f"Number of events with  1 selected g-photon : {np.sum( numPhotons==1 )} / {len(numPhotons)} " )
        print(f"Number of events with  2 selected g-photon : {np.sum( numPhotons==2 )} / {len(numPhotons)} " )
        print(f"Number of events with >2 selected g-photon : {np.sum( numPhotons >2 )} / {len(numPhotons)} " )
        
        defnString   = 'auto to_eval = std::string("numPhotons[") + std::to_string(rdfentry_) + "]"; '
        defnString  += 'return float(TPython::Eval(to_eval.c_str()));'

     #   defnString = "Numba::getNPhotons(rdfentry_)"
        rdfOut = rdfIn.Define("nGPhotons",defnString)
        nOut = rdfOut.AsNumpy(['nGPhotons'])
        print( ' in : ',np.unique(numPhotons) )
        print( ' out : ',np.unique(nOut['nGPhotons']))
        print("Imposing the filter ",f"nGPhotons < {nGPhotonMax+1}")
        rdfOut = rdfOut.Filter( f"nGPhotons < {nGPhotonMax+1}" )
        
        print(f"Number of events  {rdfIn.Count().GetValue()} -> {rdfOut.Count().GetValue()}")
        outFname = dest +'/'+ fname.split('/')[-1].replace(".root",'_cleaned.root')
        print(f"Writing Outputs to {treeName} in {outFname} [ {rdfOut.Count().GetValue()} events]")
        rdfOut.Snapshot(treeName,outFname)
        print()
