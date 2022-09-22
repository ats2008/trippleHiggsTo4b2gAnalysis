#!/usr/bin/env python3
import os,sys
import ROOT
import numpy as np

rBinCount=1
fname=""
hname=""
if len(sys.argv) > 2:
    fname=sys.argv[1]
    hnameROC=sys.argv[2]
else :
    print("Exiting  ! ")
    exit(0)
if len(sys.argv) > 4:
    rBinCount=int(sys.argv[4])

print("Getting File : ",fname)
print("Getting Hist : ",hname)

f=ROOT.TFile.Open(fname)
hist=f.Get(hnameROC)

print("ROC = ",hist.Integral("width"), "[ ",hist.Integral(), " ]")

