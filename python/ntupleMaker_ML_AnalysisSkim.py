##################################
from __future__ import print_function
import warnings
warnings.filterwarnings("ignore")
##################################

from collections import OrderedDict
import ROOT 
import numpy as np
from trippleHiggsUtils import *
import Util as utl
from branches import *
from array import array
from trippleHiggsSelector import *
from TMVA_Model import *

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

headers=getListOfStringsFromConfigs(cfgTxt,"#HEADER_BEG","#HEADER_END")
for header in headers:
    print("Loading cfg file ",header)
    f=open(header,'r')
    tmp=f.readlines()
    f.close()
    for l in tmp:
        cfgTxt.append(l)

def insideTheEllipse( x,y,x1,y1,x2,y2,a):
    return np.sqrt( (x1-x)*(x1-x)+(y1-y)*(y1-y) ) + np.sqrt( (x2-x)*(x2-x)+(y2-y)*(y2-y) ) < a           

allFnames=getListOfStringsFromConfigs(cfgTxt,"#FNAMES_BEG","#FNAMES_END")
foutName=getValueFromConfigs(cfgTxt,"OutpuFileName","fggHists.root")
processID=getValueFromConfigs(cfgTxt,"processID",default="DATA").lower()
treeName=getValueFromConfigs(cfgTxt,"treeName",default="tagsDumper/trees/Data_13TeV_TrippleHTag_0")
outTreeName=getValueFromConfigs(cfgTxt,"outTreeName",default="Data_13TeV_TrippleHTag_0")
maskSignalMgg=getBoolFromConfigs(cfgTxt,"maskSignalMgg",default=True)
doBjetCounting=getBoolFromConfigs(cfgTxt,"doBjetCounting",default=False)
minNBjetsFromMC=int(getValueFromConfigs(cfgTxt,"minNBjetsFromMC",default="0"))
maxNBjetsFromMC=int(getValueFromConfigs(cfgTxt,"maxNBjetsFromMC",default="1000"))
weightScale=float(getValueFromConfigs(cfgTxt,"WeightScale",default="1.0"))
doPtReWeighting=getBoolFromConfigs(cfgTxt,"doPtReWeighting",default=False)
pTReweitingFile=getValueFromConfigs(cfgTxt,"pTReweitingFile",default="")
pTReweitingHistName=getValueFromConfigs(cfgTxt,"pTReweitingHistName",default="")
pTReweitingValFile=getValueFromConfigs(cfgTxt,"pTReweitingValFile",default="")
pTReweitingHistValName=getValueFromConfigs(cfgTxt,"pTReweitingHistValName",default="")
resetWeight=float(getValueFromConfigs(cfgTxt,"resetWeight",default=-1e5))
mlScoreTag=getValueFromConfigs(cfgTxt,"mlScoreTag","")
etaMax=float(getValueFromConfigs(cfgTxt,"etaMax",default=2.5))
pTMin=float(getValueFromConfigs(cfgTxt,"pTMin",default=25.0))
treesToWrite=getValueFromConfigs(cfgTxt,"TreesToStore","*").split(',')
print("allFnames   :  ",              allFnames)
print("foutName   :  ",               foutName)
print("processID   :  ",              processID)
print("treeName   :  ",               treeName)
print("maskSignalMgg   :  ",          maskSignalMgg)
print("doBjetCounting   :  ",         doBjetCounting)
print("minNBjetsFromMC   :  ",        minNBjetsFromMC)
print("maxNBjetsFromMC   :  ",        maxNBjetsFromMC)
print("resetWeight   :  ",            resetWeight)
print("weightScale   :  ",            weightScale)
print("doPtReWeighting   :  ",        doPtReWeighting)
print("pTReweitingFile   :  ",        pTReweitingFile)
print("pTReweitingHistName   :  ",    pTReweitingHistName)
print("pTReweitingValFile   :  ",     pTReweitingValFile)
print("pTReweitingHistValName   :  ", pTReweitingHistValName)


scaler=mcScaler()
scalerVal=mcScaler()
if doPtReWeighting:
    scaler.setSFHistFromFile(pTReweitingFile,pTReweitingHistName)
    scalerVal.setSFHistFromFile(pTReweitingValFile,pTReweitingHistValName)
print("")

maxEvents=-1
tmp_=getValueFromConfigs(cfgTxt,"MaxEvents")
if tmp_!='':
    maxEvents=int(tmp_)

for i in allFnames:
    print(" file : ",i)
if(maxEvtsSuperSeeder > 0):
    maxEvents=maxEvtsSuperSeeder
print("maxevents : ",maxEvents)


th1Store={}

nNoHtoGammaGamma=0
nNoHHto4B=0
totalEvents=0

filemode="RECREATE"
fout = ROOT.TFile(foutName, filemode)
dir_ = fout.mkdir('tagsDumper')
dir_ = dir_.mkdir('trees')

branches_skimmed=['M1jj','M2jj','CMS_hgg_mass','dZ','weight',
                  'topScored_0','topScored_1','topScored_2','topScored_3',
                  'topScored_0_deepCSVScore','topScored_1_deepCSVScore','topScored_2_deepCSVScore','topScored_3_deepCSVScore',
                  'sumScoresTop4','sumScoresTop3','nJetVld']

branches_skimmed=np.unique(branches_skimmed)
ntuple={}
ntuple['all'] = ROOT.TNtuple(outTreeName+'_NOTAG', outTreeName, ':'.join(branches_skimmed))
tagList=['']
if processID=='data':
    tagList.append('fitRange')

catList=['cat0','cat1','cat2']
cutList = [ 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9  ]
catListToIgnore=['cat2','all']

for tag in tagList:
    if filemode=='update':
        for cat in catList:
            ntuple[cat+tag] = fout.Get(outTreeName+'_'+cat)
    else:
        for cat in catList:
            ntuple[cat+tag] = ROOT.TNtuple(outTreeName+tag+'_'+cat, outTreeName+tag+'_'+cat, ':'.join(branches_skimmed))

histStore={}
for tag in tagList:
    histStore[tag]={}
    for cat in catList:
        histStore[tag][cat]={}
        histStore[tag][cat]['jets_score']   = ROOT.TH1F( "jets_score" , "jets_score",22,-0.05,1.05 )
        histStore[tag][cat]['sortedIndexVsScore'] = ROOT.TH2F("sortedIndexVsScore" , "sortedIndexVsScore" ,20,0.0,1.0, 9, -0.5,8.5)
        histStore[tag][cat]['sortedIndexVsCumuScore'] = ROOT.TH2F("sortedIndexVsCumuScore" , "sortedIndexVsCumuScore" ,80,0.0,4.0, 9, -0.5,8.5)
        
        for k in range(8):
            histStore[tag][cat]['topScored_'+str(k)+'_ptIndex']  = ROOT.TH1F('topScored_'+str(k)+'_ptIndex','topScored_'+str(k)+'_ptIndex',9,-0.5,8.5) 
            histStore[tag][cat]['topScored_'+str(k)+'_btagIndex']= ROOT.TH1F('topScored_'+str(k)+'_btagIndex','topScored_'+str(k)+'_btagIndex',9,-0.5,8.5)
            histStore[tag][cat]['topScored_'+str(k)+'_deepCSVScore']= ROOT.TH1F('topScored_'+str(k)+'_deepCSVScore','topScored_'+str(k)+'_deepCSVScore',50,-2.5,2.5)
            histStore[tag][cat]['topScored_'+str(k)+'_score']    = ROOT.TH1F('topScored_'+str(k)+'_score','topScored_'+str(k)+'_score',20,0.0,1.0)
            histStore[tag][cat]['topScored_'+str(k)+'_scoreVsMgg'] = ROOT.TH2F( 'topScored_'+str(k)+'_scoreVsMgg', 'topScored_'+str(k)+'_scoreVsMgg',40,100,180,10 ,0.0,1.0)
            histStore[tag][cat]['topScored_'+str(k)+'_cumuSum']   = ROOT.TH1F('topScored_'+str(k)+'_cumuSum','topScored_'+str(k)+'_cumuSum',80,0.0,4.0) 
        
        for cut in cutList:
            ct=str(cut).replace('.','p')
            histStore[tag][cat]['min5Jets_'+ct+'_mgg'] = ROOT.TH1F('min5Jets_'+ct+'_mgg','min5Jets_'+ct+'_mgg',400,0.0,200.0) 
            histStore[tag][cat]['n4Jets_'+ct+'_mgg'  ] = ROOT.TH1F('n4Jets_'+ct+'_mgg','n4Jets_'+ct+'_mgg',400,0.0,200.0)
            histStore[tag][cat]['n3Jets_'+ct+'_mgg'  ] = ROOT.TH1F('n3Jets_'+ct+'_mgg','n3Jets_'+ct+'_mgg',400,0.0,200.0)
            histStore[tag][cat]['max3Jets_'+ct+'_mgg'] = ROOT.TH1F('max3Jets_'+ct+'_mgg','max3Jets_'+ct+'_mgg',400,0.0,200.0)
histStore['BiasStudy']={}
for cat in catList:
    histStore['BiasStudy'][cat]={}
    for k in range(8):
        histStore['BiasStudy'][cat]['topScored_'+str(k)+'_scoreVsMgg'] = ROOT.TH2F( 'topScored_'+str(k)+'_scoreVsMgg', 'topScored_'+str(k)+'_scoreVsMgg',40,100,180,10 ,0.0,1.0)

skimmedDataDict = OrderedDict(zip(branches_skimmed, [np.nan]*len(branches_skimmed)))

print("len(branches) : " , len(branches))

m1m2 = ROOT.TH2F("m1jj_m2jj","H1bb , H2bb mass",300,0.0,300.0,300,0.0,300. )
mvaAll=ROOT.TH1F("mvaPk_v1_all","",24,-0.1,1.1)
mvaPass=ROOT.TH1F("mvaPk_v1_pass","",24,-0.1,1.1)


sumWeights=ROOT.TH1F("sumEvts","sumEvts",1,0.0,1.0)
sumWeights.SetCanExtend(ROOT.TH1.kAllAxes)

for fname in allFnames:
    print("Opening file : ",fname)
    simFile = ROOT.TFile(fname,'READ')
    eTree=simFile.Get(treeName)
    print(" NEntries = ", eTree.GetEntries() )
    if not eTree:
        eTree=simFile.Get('tagsDumper/trees/Data_13TeV_TrippleHTag_0')
    maxEvents_ = eTree.GetEntries()
    if(maxEvents >0  and (totalEvents+maxEvents_) > maxEvents):
        maxEvents_= (maxEvents - totalEvents)
    totalEvents+=maxEvents_
    allBranches=[]
    for ky in eTree.GetListOfBranches():
        allBranches.append(ky.GetName())

    for i in range(maxEvents_):
        eTree.GetEntry(i)
        tofill={}
        
        for ii in range(4):
            tofill['topScored_'+str(ii)]=-0.05
            tofill['topScored_'+str(ii)+'_deepCSVScore']=-3.0

        tofill['sumScoresTop4']=0.0
        tofill['sumScoresTop3']=0.0
        tofill['nJetVld']=0.0
        
        if(i%500==0):
            print("   Doing i = ",i," / ",maxEvents_,
                  " n No HiggsCand : ", nNoHtoGammaGamma,
                  " n No 4BCands  : ", nNoHHto4B)
            print(" gg mass , h1MassVsh2Mass  : ",eTree.CMS_hgg_mass ,"( ",eTree.M1jj , eTree.M2jj,")" )
            
        isMasked =  processID=="data" and ( eTree.CMS_hgg_mass > 115.0) and (eTree.CMS_hgg_mass < 135.0)
        
        tag=''
        if isMasked:
            tag='fitRange'

        wei=eTree.weight
        if resetWeight > -1e4:
            wei=resetWeight
        if weightScale > -1e4:
            wei*=weightScale
        if doPtReWeighting:
            pT=eTree.diphoton_pt
            scaleFactor=scaler.getSFForX(pT)
            wei*=scaleFactor

        if processID=="MC":
            if doBjetCounting:
                nBs=getNBsFromCand(eTree)
                if nBs > maxNBjetsFromMC:
                    continue
                if nBs < minNBjetsFromMC:
                    continue
        ## Control and signal region
        sumWeights.Fill('total', 1)
        sumWeights.Fill('total_wei', wei)
        m1m2.Fill(eTree.M1jj,eTree.M2jj)
        cat='cat0'
        
        ##########   Categorization  ###################

        jetScores=[]
        jetBTagScores=[]
        jetIdx=[]
        nJetVld=0
        for i in range(8):
            tofill['jet_'+str(i)+'_isValid']=getattr(eTree , 'jet_'+str(i)+'_isValid'  )
            if abs(getattr(eTree,'jet_'+str(i)+'_eta') ) > etaMax:
                tofill['jet_'+str(i)+'_isValid']=0
            if abs(getattr(eTree,'jet_'+str(i)+'_pt') )  < pTMin:
                tofill['jet_'+str(i)+'_isValid']=0
            if tofill['jet_'+str(i)+'_isValid'] <0.25:
                continue
            nJetVld+=1
            jetIdx.append(i)
            print(eTree.event ," |  score[",i,"] -->  ",  'jet_'+str(i)+mlScoreTag+'_score' )  )
            jetScores.append( getattr(eTree , 'jet_'+str(i)+mlScoreTag+'_score' ) )
            jetBTagScores.append( getattr(eTree , 'jet_'+str(i)+'_deepCSVScore' ) )
            histStore[tag][cat]['jets_score'].Fill(jetScores[-1],wei)
        
        tofill['nJetVld']=float(nJetVld    )
        jetScores=np.array(jetScores)
        srtIdx=np.argsort(jetScores*-1)
        bTagSrtIdx=np.argsort(jetBTagScores)
        
        if len(jetScores)>=4:
            tofill['sumScoresTop4']=sum(jetScores[srtIdx[:4]])
        elif len(jetScores)>=3:
            tofill['sumScoresTop3']=sum(jetScores[srtIdx[:3]])
        
        if nJetVld < 3:
            continue
        if tofill['sumScoresTop4'] > 3.2:
            cat='cat0'
        elif tofill['sumScoresTop4'] > 0.6 :
            cat='cat1'
        else :
            cat='cat2'

        ########### Histograms  ##################
        #print()
        #print('jetScores : ',jetScores)
        #print('srtIdx : ',srtIdx)
        #print('bTagSrtIdx : ',bTagSrtIdx)
        cumuScoreSum=0.0
        for k,i in enumerate(srtIdx):
            histStore[tag][cat]['topScored_'+str(k)+'_ptIndex'].Fill(i,wei)
            histStore[tag][cat]['topScored_'+str(k)+'_btagIndex'].Fill(bTagSrtIdx[k],wei)
            histStore[tag][cat]['topScored_'+str(k)+'_score'].Fill( jetScores[i]  , wei )
            histStore[tag][cat]['topScored_'+str(k)+'_deepCSVScore'].Fill( jetBTagScores[i]  , wei )
            histStore[tag][cat]['topScored_'+str(k)+'_scoreVsMgg'].Fill(eTree.CMS_hgg_mass , jetScores[i]  ,wei)
            
            cumuScoreSum+=jetScores[i]
            tofill['topScored_'+str(k)]=jetScores[i]
            tofill['topScored_'+str(k)+'_deepCSVScore']=jetBTagScores[i]

            histStore[tag][cat]['topScored_'+str(k)+'_cumuSum'].Fill( cumuScoreSum  , wei )
            
            histStore['BiasStudy'][cat]['topScored_'+str(k)+'_scoreVsMgg'].Fill(eTree.CMS_hgg_mass , jetScores[i]  ,wei)
            histStore[tag][cat]['sortedIndexVsScore'].Fill( jetScores[i] , k , wei )
            histStore[tag][cat]['sortedIndexVsCumuScore'].Fill( cumuScoreSum , k , wei )
        
        for cut in cutList:
            n= sum(jetScores>cut)
            ct=str(cut).replace('.','p')
            if   n  > 4  :
                histStore[tag][cat]['min5Jets_'+ct+'_mgg'].Fill( eTree.CMS_hgg_mass , wei  )
            elif n == 4  :
                histStore[tag][cat]['n4Jets_'+ct+'_mgg'  ].Fill( eTree.CMS_hgg_mass , wei  )
            elif n==3  :
                histStore[tag][cat]['n3Jets_'+ct+'_mgg'  ].Fill( eTree.CMS_hgg_mass , wei  )
            else       :
                histStore[tag][cat]['max3Jets_'+ct+'_mgg'].Fill( eTree.CMS_hgg_mass , wei  )
 
        ################################################

        for ky in skimmedDataDict:
            if ky in allBranches:
                skimmedDataDict[ky]=getattr(eTree,ky)
            elif ky=='dZ' and processID=='data':
                continue
            elif ky in tofill:
                skimmedDataDict[ky]=tofill[ky]
            else:
                skimmedDataDict[ky]=-1e9
                print("Requesting  branch not available in the tree ! ",ky)
                exit(1)
        
        skimmedDataDict['weight']=wei

        sumWeights.Fill('total_'+cat, 1)
        sumWeights.Fill('total_wei'+cat, wei)
        ntuple[cat].Fill(array('f', skimmedDataDict.values()))
        if (tag!='') and ( cat not in catListToIgnore ):
            ntuple[cat+tag].Fill(array('f', skimmedDataDict.values()))
        if 'all' not in catListToIgnore:
            ntuple['all'].Fill(array('f', skimmedDataDict.values()))
        

        ########## ------------- ##################

    simFile.Close()           
    print("Closing file : ",fname)

dir_.cd()    
sumWeights.Write()
m1m2.Write()
mvaAll.Write()
mvaPass.Write()

for cat in ntuple:
    if "*" in treesToWrite :
        ntuple[cat].Write()
    else:
        if cat in catListToIgnore:
            continue
        ntuple[cat].Write()

histDir=fout.mkdir('Histograms')

for tag in histStore:
    if tag=='':
        tagDir=histDir.mkdir('all')
    else:
        tagDir=histDir.mkdir(tag)
    #tagDir.cd()
    for cat in histStore[tag]:
        catDir=tagDir.mkdir(cat)
        catDir.cd()
        keys=sorted(histStore[tag][cat].keys())
        for ky in keys:
            #histStore[tag][cat][ky].SetName(tag+cat+'_'+ky)
            histStore[tag][cat][ky].Write()

fout.Purge()
fout.Close()

print(" File written out  : ",foutName)

