import uproot as urt
import ROOT
import json,sys,os,argparse
import numpy as np
import prettytable as ptab

import scaleFactorUtil as scl

fileIn="/grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashggFinalFits/CMSSW_10_2_13/src/flashggFinalFit/Trees2WS/oct16/singleH.root"
folderIn='trees'

print("Opening {fileIn}")
f=urt.open(fileIn)

print("\tOpening {folderIn}")
print(f.keys())
treeCollection=f[folderIn]

for ky in treeCollection:
    iTree=treeCollection[ky]
    dataStore=iTree.arrays(['weight'])
    print(ky,"\t\t",len(dataStore["weight"]),f"\t\t{np.sum(dataStore['weight'])}")
f.close()
exit()

