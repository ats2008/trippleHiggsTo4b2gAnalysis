import os, sys
import ROOT
import numpy as np


def saveTheDictionary(aCollection,fname=None,folder=None,fileExists=False):
    outfile=folder
    if fname!=None and outfile==None:
        outfile=ROOT.TFile(fname,"RECREATE")
    outfile.cd()
    
    for key in aCollection:
        if type(aCollection[key])==type({}):
            inFolder=outfile.mkdir(key)
            saveTheDictionary(aCollection[key],None,inFolder)
        else:
            aCollection[key].Write()
    if fname!=None:
        outfile.Close()

def dPhi(phi1,phi2):
    return np.arccos(np.cos(phi2-phi1))

def dEta(eta1,eta2):
    return abs(eta1-eta2)

def deltaR(eta1,phi1,eta2,phi2):
    dphi=np.arccos(np.cos(phi2-phi1))
    return np.sqrt((eta1-eta2)*(eta1-eta2) + dphi*dphi)

def getListOfStringsFromConfigs(cfgTxt,beg,end):
    isFlist=False
    allFnames=[]
    for l in cfgTxt:
        if beg in l:
            isFlist=True
            continue
        if end in l:
            isFlist=False
            continue
        if isFlist:
            allFnames.append(l[:-1])
    return allFnames

def getValueFromConfigs(cfgTxt,tag,default=""):
    isParams=False
    val=default
    for l in cfgTxt:
        if "#PARAMS_BEG" in l:
            isParams=True
            continue
        if '#PARAMS_END' in l:
            isParams=False
            continue
        if isParams:
            ll=l[:-1].split('=')
            if tag==ll[0]:
                val=ll[1]
    return val

def getBoolFromConfigs(cfgTxt,tag,default):
    isParams=False
    val=''
    for l in cfgTxt:
        if "#PARAMS_BEG" in l:
            isParams=True
            continue
        if '#PARAMS_END' in l:
            isParams=False
            continue
        if isParams:
            ll=l[:-1].split('=')
            if tag==ll[0]:
                val=ll[1]
    if val=="True" or val=="1" or val=='true':
        return True
    if val=="False" or val=="0" or val=='false':
        return False
    else:
        return default


def getTheObjectsFromFile(aFile):
    histStore={}
    for key in aFile.GetListOfKeys():
        kk=key.GetName()
        curObj=aFile.Get(kk)
        if curObj.Class_Name()=='TDirectoryFile':
            histStore[kk]=getTheObjectsFromFile(curObj)
        else:
            histStore[kk]=curObj       
    return histStore

def getSumHistDicts(histList):
    histOut={}
    for kk in histList[0]:
        curObj=histList[0][kk]
        if type(curObj) == type({}):
            histOut[kk]=getSumHistDicts([ hDict[kk] for hDict in histList ])
        else:
            histOut[kk]=None
#             print(kk)
            for hDict in histList:
                if histOut[kk]==None:
                    histOut[kk]=hDict[kk].Clone()
                else:
                    histOut[kk].Add(hDict[kk])
    return histOut

def rebinTheHistogram(hist,n,binEdges=None):

    if binEdges==None:
        bEdges=0    
    else:
        bEdges=np.array(binEdges)

    hout=hist.Rebin(n,hist.GetName()+"_rbinned",bEdges)
    return hout

def getAClonedTree(eTree,name=""):
    if name=="":
        name=eTree.GetName()
    oTree= ROOT.TTree(name)
    branches=OrderedDict(eTree.__dict__.keys())
    tofill = OrderedDict(zip(branches, [np.nan]*len(branches)))
    
    return oTree

