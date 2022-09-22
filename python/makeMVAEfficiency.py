
import os,sys
import ROOT
import numpy as np

rBinCount=1
fname=""
hname=""
if len(sys.argv) > 2:
    fname=sys.argv[1]
    hname=sys.argv[2]
else :
    print("Exiting  ! ")
    exit(0)
if len(sys.argv) > 3:
    rBinCount=int(sys.argv[3])

print("Getting File : ",fname)
print("Getting Hist : ",hname)

f=ROOT.TFile.Open(fname)
histSigEff=f.Get(hname)

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

effsToSee=[0.9]
valsToSee=[0.1,0.05,0.01,0.005,0.001,0.0001]
cutTosee=[0.9,0.99]

cumu=histSigEff.GetCumulative()
#cumu.Draw()
#c=3
#input(c)
cumu.Scale( 1.0/histSigEff.Integral() )
for i in range(nBins):
    cut     = cumu.GetBinCenter(i)
    sigEff  = cumu.GetBinContent(i)
    print(cut,sigEff)

#    for ct in cutTosee:
#        if abs(cut - ct ) < 0.005:
#            print("\t for ct ",ct," | ",abs(cut - ct ) < 0.05)
#            print("cut place for cut " ,ct ," |  acutual cut ",cut," rejection eff ",sigEff)
#    
#    for eff in effsToSee:
#        if (abs((1.0-sigEff) - eff ) < 0.01):
#            print("\t for eff ",eff," | ",abs(  (1.0-sigEff) - eff ) < 0.05)
#            print("eff place for eff " ,eff ," |  acutual cut ",cut," acceptance  eff ",1.0-sigEff, " | ",abs(cut - ct ) < 0.05)
    
