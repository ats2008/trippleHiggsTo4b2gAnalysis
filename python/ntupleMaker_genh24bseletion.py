from __future__ import print_function
import ROOT 
import numpy as np
from trippleHiggsGenAnlayzer import *
from Util  import *

import os,sys

cfgFileName=''
if len(sys.argv) <2:
    print("Usage\n\t ~$python recoAnalyzer.py <configFile>\n")
    exit(1)
else:
    cfgFileName=sys.argv[1]
maxEvtsSuperSeeder=-1
if len(sys.argv) >2:
    maxEvtsSuperSeeder=int(sys.argv[2])
print("Loading cfg file ",cfgFileName)
f=open(cfgFileName,'r')
cfgTxt=f.readlines()
f.close()

allFnames=getListOfStringsFromConfigs(cfgTxt,"#FNAMES_BEG","#FNAMES_END")
foutName=getValueFromConfigs(cfgTxt,"OutpuFileName")
treeName=getValueFromConfigs(cfgTxt,"treeName",default="tagsDumper/trees/Data_13TeV_TrippleHTag_0")

maxEvents=-1
tmp_=getValueFromConfigs(cfgTxt,"MaxEvents")
if tmp_!='':
    maxEvents=int(tmp_)

for i in allFnames:
    print(" file : ",i)
if(maxEvtsSuperSeeder > 0):
    maxEvents=maxEvtsSuperSeeder
print("maxevents : ",maxEvents)

nNoHtoGammaGamma=0
nNoHHto4B=0
totalEvents=0

isMC = True

sumEntries=ROOT.TH1F("sumEvts","sumEvts",1,0.0,1.0)
sumEntries.SetCanExtend(ROOT.TH1.kAllAxes)
sumWeights=ROOT.TH1F("sumWeighs","sumWeighs",1,0.0,1.0)
sumWeights.SetCanExtend(ROOT.TH1.kAllAxes)

th1Store={}
th1Store['acceptance2d']=ROOT.TH1F("acceptance2d","acceptance2d",1,0.0,1.0)
th1Store['acceptance2d'].SetCanExtend(ROOT.TH1.kAllAxes)

bLV=[]
for i in range(4):
    bLV.append(ROOT.TLorentzVector())
gammaLV=[]
for i in range(4):
    gammaLV.append(ROOT.TLorentzVector())

for fname in allFnames:
    
    simFile = ROOT.TFile(fname,'READ')
    eTree=simFile.Get(treeName)
    if not eTree:
        eTree=simFile.Get('Ntuple/tree')
        
    maxEvents_ = eTree.GetEntries()
    if(maxEvents >0  and (totalEvents+maxEvents_) > maxEvents):
        maxEvents_= (maxEvents - totalEvents)
    totalEvents+=maxEvents_

    for i in range(maxEvents_):
        eTree.GetEntry(i)
        sumEntries.Fill('allEvents', 1)
        bIdx=0
        gammaIdx=0        
        
        if(i%100==0):
            print("Doing i = ",i," / ",maxEvents_)
        

        bDaus=getHiggsDauP4s(eTree,5)       
        bIdx=len(bDaus)
        gammaDaus=getHiggsDauP4s(eTree,22)       
        gammaIdx=len(gammaDaus)
        
        bLV=bDaus
        gammaLV=gammaDaus
        if( (gammaLV[0]+gammaLV[1]).Pt() < 100 ): continue;

        if(gammaIdx!=2 or bIdx!=4):
            print("Problem !!! ",gammaIdx,bIdx)
        
        nHiggsAbove300=0
        
        if ( bLV[0] +bLV[1] ).Pt() > 300: nHiggsAbove300+=1
        if ( bLV[2] +bLV[3] ).Pt() > 300: nHiggsAbove300+=1
        if ( gammaLV[0] +gammaLV[1] ).Pt() > 300: nHiggsAbove300+=1
        if ( gammaLV[0].Pt() < 30.0 ): continue
        if ( gammaLV[1].Pt() < 20.0 ): continue
        if ( abs(gammaLV[0].Eta()) > 2.5 ): continue
        if ( abs(gammaLV[1].Eta()) > 2.5 ): continue

        sumEntries.Fill('allEventsAfterHgg', 1)
        
        isB0InAcceptance=False
        isB1InAcceptance=False
        isB2InAcceptance=False
        isB3InAcceptance=False

        if ( abs( bLV[0].Eta() ) < 4.7  and abs( bLV[0].Pt() ) > 25.0 ) :
            isB0InAcceptance=True
            sumEntries.Fill('isB0InAcceptance', 1)
        if ( abs( bLV[1].Eta() ) < 4.7  and abs( bLV[1].Pt() ) > 25.0 ) :
            isB1InAcceptance=True
            sumEntries.Fill('isB1InAcceptance', 1)
        if ( abs( bLV[2].Eta() ) < 4.7  and abs( bLV[2].Pt() ) > 25.0 ) :
            isB2InAcceptance=True
            sumEntries.Fill('isB2InAcceptance', 1)
        if ( abs( bLV[3].Eta() ) < 4.7  and abs( bLV[3].Pt() ) > 25.0 ) :
            isB3InAcceptance=True
            sumEntries.Fill('isB3InAcceptance', 1)
        nout=0
        if not isB0InAcceptance:            nout+=1;    
        if not isB1InAcceptance:            nout+=1;    
        if not isB2InAcceptance:            nout+=1;    
        if not isB3InAcceptance:            nout+=1;
        
        nleadOut=0
        if not isB0InAcceptance:            nleadOut+=1;    
        if not isB2InAcceptance:            nleadOut+=1;    
        
        nSubLeadOut=0
        if not isB1InAcceptance:            nSubLeadOut+=1;    
        if not isB3InAcceptance:            nSubLeadOut+=1;    
        
        h1Out= not ( isB0InAcceptance and isB1InAcceptance)
        h2Out= not ( isB2InAcceptance and isB3InAcceptance)
        
        sumEntries.Fill(str(nout)+'_JetsNotInAcceptance',1)
        sumEntries.Fill(str(nleadOut)+'_LeadNotInAcceptance',1)
        sumEntries.Fill(str(nSubLeadOut)+'_SubLeadNotInAcceptance',1)
        if h1Out:
            sumEntries.Fill('H1NotInAcceptance', 1)
        if h2Out:
            sumEntries.Fill('H2NotInAcceptance', 1)
            
        if isB0InAcceptance and isB1InAcceptance and isB2InAcceptance and isB3InAcceptance :
            sumEntries.Fill('allBJetsInAcceptance', 1)
        else:
            pass
            #print(  isB0InAcceptance  ,isB1InAcceptance , isB2InAcceptance ,isB3InAcceptance ,
            #        " --> ", isB0InAcceptance and isB1InAcceptance and isB2InAcceptance and isB3InAcceptance )     
    simFile.Close()           

ofile=ROOT.TFile( foutName  ,"RECREATE")
filemode="RECREATE"
fout = ROOT.TFile(foutName, filemode)
dir_ = fout.mkdir('tagsDumper')
dir_ = dir_.mkdir('trees')
dir_.cd()

sumEntries.Write()
sumWeights.Write()
print(" File written out  : ",foutName)
