#!/usr/bin/env python
# coding: utf-8

"""
    
    for i in */* ; do echo $i ; hadd $i.root $i/*/*/*  ; done

"""
import matplotlib.pyplot as plt
import uproot as urt
import ROOT
import json
import numpy as np
import os,argparse
import pickle as pkl
import Util as utl
import sys,os

#sys.path.append('/home/aravind/cernbox/work/trippleHiggs/hhhTo4b2gamma/genAnalysis/python/analysis/python/')

import glob


prefix ='./'
cutFname=None #'cuts/PreselectionRaw.cuts'

# 1.3146618961745309 scale_msbr=0.9455667969895265
scaleFactors={
    'run2' : {
        'bkg' : {
            'morphedData'    : 0.946,  #0.939 ,#  0.92 , # 0.934   , #0.946  , 1.0633791300781563 scale_msbr=0.9390015334810019
            'ggBox1Bjet'     : 1.311,  #1.063 ,#  0.9788 , # 1.247 , #1.311  ,
            'ggBox2Bjet'     : 1.311,  #1.063 ,#  0.9788 , # 1.247 , #1.311  ,
            'ggBox'          : 1.311,  #1.063 ,#  0.9788  # 1.247  , #1.311  ,
        }
    }
}

fileDict={
  'run2' : {
    'bkg' : {
    'morphedData' :'data/morphedData.root',
    'ggBox'       : 'bkg/mc_ggBox.root',
    'ggBox1Bjet'  : 'bkg/mc_ggBox1Bjet.root',
    'ggBox2Bjet'  : 'bkg/mc_ggBox2Bjet.root'
   }
 }
}

# In[4]:


def getTH1FromNumpHist(nph,nph_w2=None):
    hist = ROOT.TH1D("","",len(nph[1])-1,np.array(nph[1],dtype=float))
    hist.Sumw2()
    for i,freq in enumerate(nph[0]):
        hist.SetBinContent(i+1,freq)
    if nph_w2!=None:
        for i,freq_w2 in enumerate(nph_w2[0]):
            hist.SetBinError(i+1,np.sqrt(freq_w2))
    return hist


# In[5]:


mcDict={
    'QCD':['qcd30To40','qcd40ToInf'],
    'gJet':['gJet20To40','gJet40ToInf'],
    'ttX' :['ttjj','ttgj','ttgg'],
    'ggX' :['ggBox','ggBox1Bjet','ggBox2Bjet'],
    'morphedSBR' : ['morphedData']
}

color_scheme={
    'QCD' :'c',
    'gJet':'g',
    'ggX' :'r',
    'ttX' :'b',
    'morphedSBR' : 'y'
    
}

label_scheme={
    'QCD' :'QCD',
    'gJet':r'$\gamma$ + Jets',
    'ttX' :r't$\bar{t}$ + X',
    'ggX' :r'$\gamma\gamma$ + Jets',
    'morphedSBR' : 'Morphed SBR'
    
}


# In[7]:



yearsToProcess_all=['2018','2017','2016PreVFP','2016PostVFP','run2','2016']
yearsToProcess=['run2'] #yearsToProcess_all
bkgToProcess=['all', 'ggBox1Bjet','ggBox2Bjet', 'ggBox','gJet20To40','gJet40ToInf','qcd30To40', 'qcd40ToInf']

# if 'all' not in args.year:
#     yearsToProcess=[]
#     for yr in args.year.split(","):
#         if yr not in yearsToProcess_all:
#             print(yr," not in catalogue. Skipping !! ")
#             continue
#         yearsToProcess.append(yr)


# In[9]:


unblind=True
cutsToApply=[]
cutsKey=[]
cutIdx=0
txt=[]
if cutFname:
    with open( cutFname ,'r') as f:

        txt_=f.readlines()
        for l in txt_:
            if l[0]=='#':
                continue
            if len(l) <2: 
                continue
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



# In[10]:


rdataFrames={}
for yr  in fileDict:
    rdataFrames[yr]={}
    if 'sig' in fileDict[yr]:
        rdataFrames[yr]['sig']={}
        fileName=fileDict[yr]['sig']['ggHHH']
        treeName = "trees/ggHHH_125_13TeV"
        print("Registering datset : ggHHH , ",yr," withh tree",treeName)
        print("\t\tFilename :", fileName)
        rdataFrames[yr]['sig']['ggHHH']= ROOT.RDataFrame(treeName, fileName)
        n =rdataFrames[yr]['sig']['ggHHH'].Count().GetValue()
        for cut,cutKey in zip(cutsToApply,cutsKey):
            nDiff =rdataFrames[yr]['sig']['ggHHH'].Count().GetValue()
            rdataFrames[yr]['sig']['ggHHH'] = rdataFrames[yr]['sig']['ggHHH'].Filter(cut,cutKey)
            nDiff-=rdataFrames[yr]['sig']['ggHHH'].Count().GetValue()
            print(f"\t\t cut {cutKey} applied | killed {nDiff} Events [  cumulative eff : {1.0  - nDiff/n}]")

    if 'data' in fileDict[yr]:
        ky=list(fileDict[yr]['data'].keys())[0]
        fileName=fileDict[yr]['data'][ky]
        treeName = "trees/Data_13TeV_TrippleHTag_0"
        if unblind:
            rdataFrames[yr]['data']['data']= ROOT.RDataFrame(treeName, fileName)
            n =rdataFrames[yr]['data']['data'].Count().GetValue()
            for cut,cutKey in zip(cutsToApply,cutsKey):
                nDiff =rdataFrames[yr]['data']['data'].Count().GetValue()
                rdataFrames[yr]['data']['data'] = rdataFrames[yr]['data']['data'].Filter(cut,cutKey)
                nDiff-=rdataFrames[yr]['data']['data'].Count().GetValue()
                print(f"\t\t cut {cutKey} applied | killed {nDiff} Events [  cumulative eff : {1.0  - nDiff/n}]")
        else:
            rdataFrames[yr]['data']['data']= ROOT.RDataFrame(treeName, fileName).Filter(' CMS_hgg_mass < 115.0  || CMS_hgg_mass >135.0')
            n =rdataFrames[yr]['data']['data'].Count().GetValue()
            for cut,cutKey in zip(cutsToApply,cutsKey):
                nDiff =rdataFrames[yr]['data']['data'].Count().GetValue()
                rdataFrames[yr]['data']['data'] = rdataFrames[yr]['data']['data'].Filter(cut,cutKey)
                nDiff-=rdataFrames[yr]['data']['data'].Count().GetValue()
                print(f"\t\t cut {cutKey} applied | killed {nDiff} Events [  cumulative eff : {1.0  - nDiff/n}]")
        print("Registering datset : data , ",yr," withh tree",treeName)
        print("\t\tFilename :", fileName)
     
    if 'bkg' in  fileDict[yr]: 
        rdataFrames[yr]['bkg']={}
        for bkg in  fileDict[yr]['bkg']:
            treeName = "trees/bkg_13TeV_TrippleHTag_0"
            if bkg=='morphedData':
                treeName='trees/Data_13TeV_TrippleHTag_0'
            print("Registering datset : ",bkg," , ",yr," withh tree",treeName)
            fileName = fileDict[yr]['bkg'][bkg]
            rdataFrames[yr]['bkg'][bkg]=ROOT.RDataFrame(treeName, fileName)
            #rdataFrames[yr]['bkg'][bkg]=ROOT.RDataFrame(treeName, fileName).Filter(' CMS_hgg_mass < 115.0  || CMS_hgg_mass >135.0')
            print("\t\tFilename :", fileName)
            n=rdataFrames[yr]['bkg'][bkg].Count().GetValue()
            for cut,cutKey in zip(cutsToApply,cutsKey):
                nDiff=rdataFrames[yr]['bkg'][bkg].Count().GetValue()
                rdataFrames[yr]['bkg'][bkg] = rdataFrames[yr]['bkg'][bkg].Filter(cut,cutKey)
                nDiff-=rdataFrames[yr]['bkg'][bkg].Count().GetValue()
                print(f"\t\t cut {cutKey} applied | killed {nDiff} Events [  cumulative eff : {1.0  - nDiff/(n+1e-9)}]")


# In[14]:
if not os.path.exists(prefix):
    os.system('mkdir -p '+prefix)
    print("Making directory ",prefix)

rdataFrames_updated={}

for yr in rdataFrames:
    rdataFrames_updated[yr]={}
    for tag in rdataFrames[yr]:
        rdataFrames_updated[yr][tag]={}
        for ky in rdataFrames[yr][tag]:
            print(f" Processing {yr}/{tag}/{ky}")
            outFname=fileDict['run2'][tag][ky].split('/')[-1]
            outFname=outFname.replace('.root','_combReWeighted.root')
            outFname=prefix+'/'+outFname
            if  True:
                allBranches=[ str(i) for i in rdataFrames[yr][tag][ky].GetColumnNames()]
                dFrame = rdataFrames[yr][tag][ky]
                if 'minPhoIDMVA' not in allBranches:
                    print("## defining the minPhoIDMVA")
                    dFrame = rdataFrames[yr][tag][ky].Define("minPhoIDMVA","min(customLeadingPhotonIDMVA,customSubLeadingPhotonIDMVA)")
                dset=dFrame.AsNumpy(['weight_v0','lumi','rdfentry_','minPhoIDMVA'])
                sf=np.ones_like(dset['weight_v0'])
                print(f"\nProcessing {yr} / {tag} / {ky}")
                if yr in scaleFactors:
                    if tag in scaleFactors[yr]:
                        if ky in scaleFactors[yr][tag]:
                            sf=scaleFactors[yr][tag][ky]*sf
                            print(f"\t Setting SF = {sf}")
                mask= dset['minPhoIDMVA'] < -0.7
                print(f" Number of events in LSR : {np.sum(mask)}")
                sf[mask]=1.0
                rdfentryDef=dset['rdfentry_']
                rdfEntry=np.arange(len(rdfentryDef))
                rdfMap = { i:j for i,j in zip( rdfentryDef ,rdfEntry ) }
                defnString  = 'auto to_eval = std::string("rdfMap[") + std::to_string(rdfentry_) + "]"; '
                defnString += 'return int(TPython::Eval(to_eval.c_str()));'
                dFrame      = dFrame.Redefine("rdfentry_",defnString)

                weight_v1=dset['weight_v0']*dset['lumi']*sf
                defnString  = 'auto to_eval = std::string("weight_v1[") + std::to_string(rdfentry_) + "]"; '
                defnString += 'return float(TPython::Eval(to_eval.c_str()));'
                
                try:
                    dFrame= dFrame.Define("weight_cScaledV1",defnString)
                except:
                    print("REDEFINING CSCALE")
                    dFrame= dFrame.Redefine("weight_cScaledV1",defnString)

                #dFrame= dFrame.Define("weight_lumiScaled","weight_v0*lumi")
                
                treeName='trees/bkg_13TeV_TrippleHTag_0'
                if ky=='ggHHH':
                    treeName='trees/ggHHH_125_13TeV'
                if ky=='data':
                    treeName='trees/Data_13TeV_TrippleHTag_0'
                print(f"\tMaking {outFname}")
                dFrame.Snapshot(treeName,outFname)
            fileDict[yr][tag][ky]=outFname

with open(prefix+'/filelist.json','w') as f:
    json.dump(fileDict,f,indent=4)
