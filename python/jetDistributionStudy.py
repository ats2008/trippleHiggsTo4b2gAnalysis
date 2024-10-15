#!/usr/bin/env python

import ROOT
import json
import pickle as pkl
import numpy as np
import os,argparse
import Util as utl

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


def getObjectAsArray(raw_data,object_tag="Photon",colums=["pt","eta","phi","mass"],nMax=6,nameTag="Momentum4D"):
    objectArr=[ 
                ak.zip({var : raw_data[f"{object_tag}_{i}_{var}"]  for var in colums},with_name=nameTag)
                 for i in range(0,nMax)
            ]
    objectArr=[ ak.unflatten(objArr,counts=1,axis=-1) for objArr in objectArr]
    objectArr= ak.concatenate(objectArr,axis=1)
    return objectArr

def getActualJets(jets):
    jpt =jets.pt
    jeta=jets.eta

    padding=ak.unflatten(ak.zeros_like(jpt[:,0]),counts=1)
    jpt=np.concatenate([padding,jpt],axis=1)

    padding=ak.unflatten(ak.zeros_like(jeta[:,0]),counts=1)
    jeta=np.concatenate([padding,jeta],axis=1)

    actualJets_mask = ( jpt[:,:-1]!= jpt[:,1:]) &  ( jeta[:,:-1]!= jeta[:,1:])
    actualJets_mask = actualJets_mask & (jpt[:,1:] > 0)
    aJets = ak.drop_none(ak.mask(jets,mask=actualJets_mask))
    return aJets


def main():
    
    fileDict={
        'run2' :{
            'sig' :{},
            'bkg' :{ 
                'gjet' : '/home/aravind/Documents/dumpX/output_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8_6.root',
            },
            'data':{}
        }
    }

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
    rdataFrames={yr : { tag : {}} }
    treeName=cfg['treeName']
    
    idx=0
    for fls in cfg['files']:
        rdataFrames[yr][tag][idx] = ROOT.RDataFrame(treeName, fls)
        idx+=1   
    # In[8]:
    varsToGet=[]
    baseVars=[
            'weight',
            'CMS_hgg_mass',
            'corrMET',
            'EGMLeadingPhotonIDMVA',
            'EGMSubLeadingPhotonIDMVA',
            'customLeadingPhotonIDMVA',
            'customSubLeadingPhotonIDMVA',
            'subleadingPhoton_pt',
            'diphoton_eta',
            'diphoton_phi',
            'diphoton_pt',
            'subleadingPhoton_eta',
            'subleadingPhoton_phi',
            'leadingPhoton_pt',
            'leadingPhoton_eta',
            'leadingPhoton_phi',
            'lumi',
            'weight',
            "year"
        ]
    varsToGet+=baseVars
    tmp=[
        "jet_@@IDX_bJetRegCorr",
        "jet_@@IDX_bJetRegRes",
        "jet_@@IDX_csvScore",
        "jet_@@IDX_deepCSVScore",
        "jet_@@IDX_deepJetScore",
        "jet_@@IDX_deepJetScore_b",
        "jet_@@IDX_deepJetScore_bb",
        "jet_@@IDX_deepJetScore_lepb",
        "jet_@@IDX_eta",
        "jet_@@IDX_flavour",
        "jet_@@IDX_isLoose",
        "jet_@@IDX_isTight",
        "jet_@@IDX_isTight2017",
        "jet_@@IDX_isTight2018",
        "jet_@@IDX_isValid",
        "jet_@@IDX_mass",
        "jet_@@IDX_pFlavour",
        "jet_@@IDX_particleNetAK4_B",
        "jet_@@IDX_particleNetAK4_CvsB",
        "jet_@@IDX_particleNetAK4_CvsL",
        "jet_@@IDX_particleNetAK4_QvsG",
        "jet_@@IDX_particleNetAK4_puIdDisc",
        "jet_@@IDX_phi",
        "jet_@@IDX_pt",
        "jet_@@IDX_puJetIdMVA"
    ]
    for i in range(0,12):
        for ky in tmp:
            varsToGet+=[ky.replace('@@IDX',f'{i}')]
    
    # In[11]:
    
    
    dataStore={}
    for yr in rdataFrames:
        dataStore[yr]={}
        for ky in rdataFrames[yr][tag]:
            print(f"Loading the year {yr} / {ky}  ")
            dataStore[yr][ky]=rdataFrames[yr][tag][ky].AsNumpy(varsToGet)
            dataStore[yr][ky]['leadingPhoton_mass']=np.zeros_like(dataStore[yr][ky][varsToGet[0]])
            dataStore[yr][ky]['subleadingPhoton_mass']=np.zeros_like(dataStore[yr][ky][varsToGet[0]])
    
    
       
       
    
    # In[16]:
    
    
    eCols=["bJetRegCorr",
    "bJetRegRes",
    "csvScore",
    "deepCSVScore",
    "deepJetScore",
    "deepJetScore_b",
    "deepJetScore_bb",
    "deepJetScore_lepb",
    "eta",
    "flavour",
    "isLoose",
    "isTight",
    "isTight2017",
    "isTight2018",
    "isValid",
    "mass",
    "pFlavour",
    "particleNetAK4_B",
    "particleNetAK4_CvsB",
    "particleNetAK4_CvsL",
    "particleNetAK4_QvsG",
    "particleNetAK4_puIdDisc",
    "puJetIdMVA",
    "phi",
    "pt"]
    
    
    data_awk={}
    
    for yr in dataStore:
        data_awk[yr]={}
        dataToWrite=None
        isFirst=True
        for ky in dataStore[yr]:
            print(f"Processing {yr} / {ky}")
            dset=getMomentumArrayFromObjects(dataStore[yr][ky],obj_tagList=['leadingPhoton','subleadingPhoton'])
            jets=getObjectAsArray(dataStore[yr][ky],
                                         object_tag='jet',
                                         colums=eCols,nMax=12)
            corr_jets = getActualJets(jets)
            dset=ak.with_field(dset,what=corr_jets,where='jet')
            for var in baseVars:
                dset=ak.with_field(dset,what=ak.Array(dataStore[yr][ky][var]),where=var)

            if isFirst:
                dataToWrite=dset
                isFirst=False
            else:
                dataToWrite=ak.concatenate([dataToWrite,dset])

        print(ak.num( dataToWrite.jet )) 
        dOut= { 'data' : dataToWrite }
        oFname=f'{dest}/out_{yr}_{tag}_{hashStr}.pkl'
        
        with open(oFname,'wb') as f:
            print("\t Wring out to ",oFname)
            pkl.dump(dOut,f)    

if __name__=='__main__':
    main()    
