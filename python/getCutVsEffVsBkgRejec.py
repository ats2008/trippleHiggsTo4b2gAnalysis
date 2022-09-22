#!/usr/bin/env python3
import os,sys
import ROOT
import numpy as np

rBinCount=1
fname=""
hname=""
if len(sys.argv) > 2:
    fname=sys.argv[1]
    hnameS=sys.argv[2]
    hnameB=sys.argv[3]
else :
    print("Exiting  ! ")
    exit(0)
if len(sys.argv) > 4:
    rBinCount=int(sys.argv[4])

print("Getting File : ",fname)
print("Getting Hist : ",hname)

f=ROOT.TFile.Open(fname)
histSigEff=f.Get(hnameS)
histBkgEff=f.Get(hnameB)

nBins=histSigEff.GetNbinsX()

print(-1 ," ","dT","\t:\t", 
       "CUT" ,"\t:\t",
        "SigEff" ,'\t:\t',
       "BkgEff" ,"\t:\t",
        "significane")
vals={}
delta=0.001
vals['maxSig']=()
maxSig=-1.0

effsToSee=[0.95,0.9,0.8,0.65,0.5]
valsToSee=[0.1,0.05,0.01,0.005,0.001,0.0001]
cutTosee=[0.9,0.99]
for i in range(nBins):
    binXsig=histSigEff.GetBinCenter(i)
    binXbkg=histSigEff.GetBinCenter(i)
    sigEff=histSigEff.GetBinContent(i)
    bkgEff=histBkgEff.GetBinContent(i)
    significance = sigEff/(sigEff+bkgEff + 1e-4)**0.5
    binX=binXsig
    if significance > maxSig:
        vals['maxSig']=(binX,sigEff,bkgEff,significance)
        maxSig=significance

#    print(i," ",binXbkg-binXsig,"\t:\t", 
#            "{0:0.3f}".format(binX),"\t:\t",
#            "{0:0.3f}".format(sigEff),'\t:\t',
#            "{0:0.3f}".format(bkgEff),"\t:\t",
#            "{0:0.3f}".format(significance))
    bkgEff=np.round(bkgEff,5)
    
    inTol=False
    for eff in effsToSee:
        if ( abs(sigEff-eff) < delta ):
            inTol=True
    for eff in valsToSee:
        if ( abs(bkgEff-eff) < delta ):
            inTol=True
    if inTol:
        if bkgEff in vals:
            continue
        vals[bkgEff]=(binX,sigEff,bkgEff,significance)
i=0
for v in vals: 
    i+=1
    print(i,"\t:\t", 
               "cut {0:0.3f}".format(vals[v][0]),"\t:\t",
            "sigEff {0:0.3f}".format(vals[v][1]),'\t:\t',
            "bkgEff {0:0.3f}".format(vals[v][2]),"\t:\t",
            "signif {0:0.3f}".format(vals[v][3]),end=" ")
    if v == 'maxSig':
        print("  : ",v)
    else:
        print()
