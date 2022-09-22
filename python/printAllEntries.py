from PlotUtil import *
from Util import *
import os
from collections import OrderedDict

ext=None
if len(sys.argv) >1:
    ext=sys.argv[1]

def getBinVsContent(hist):
    xAxis=hist.GetXaxis()
    nBins=hist.GetNbinsX()
    cutFlow={}
    
    for i in range(nBins+1):
        binName=xAxis.GetBinLabel(i)
        binContent=hist.GetBinContent(i)
        if binName=="":
            continue
        cutFlow[binName]=binContent#np.round(float(binContent),5)
    return cutFlow


if __name__ == "__main__":
    
    for preReweight in [False]:
        fileDict=OrderedDict()
        fileDict[ "ggHHH"          ] =   "workarea/gHHHto4b2gamma_UL17_13TeV_accptance.root"
        fileDict[ "ggHHH"          ] =   "workarea/gHHHto4b2gamma_UL17_13TeV_accptance4e7.root"
        fileDict[ "ggHHH"          ] =   "workarea/gHHHto4b2gamma_UL17_13TeV_accptance4e7_100GeVHggPtMin.root"
        fileDict[ "ggHHH"          ] =   "workarea/batch/jetIDStudy/ggHHHto4b2gamma_UL17_13TeV_reco4e7.root"
        fileDict[ "ggHHH"          ] =   "workarea/gHHHto4b2gamma_UL17_13TeV_pvtNtuple.root"
        fileDict[ "ggHHH"          ] =   "workarea/batch/jetIDStudy/ggHHHto4b2gamma_UL17_13TeV_reco4e7_pTHggGr100.root"
        fileDict[ "ggHHH"          ] =   "workarea/batch/jetIDStudy/ggHHHto4b2gamma_UL17_13TeV_reco2p5_noHFlavmatching_0p9.root"
        fileDict[ "ggHHH"          ] =   "workarea/batch/jetIDStudy/ggHHHto4b2gamma_UL17_13TeV_reco2p5_withHFlavmatching_0p9.root"
        fileDict[ "ggHHH"          ] =   "workarea/batch/jetIDStudy/ggHHHto4b2gamma_UL17_13TeV_reco4p7_withHFlavmatching_0p9.root"
        fileDict[ "ggHHH"          ] =   "workarea/batch/jetIDStudy/ggHHHto4b2gamma_UL17_13TeV_reco2p5.root"
        fileDict[ "ggHHH"          ] =   "workarea/batch/jetStudy_accEff/ggHHHto4b2gamma_UL17_13TeV_v2eta2p5_dr0p4.root"
           
        if ext !=  None:
            fileDict[ "ggHHH"          ]= ext
        #binsToPrint=['total', 'total_wei', 'total_cat0', 'total_weicat0' , 'total_cat1', 'total_weicat1']
        binsToPrint=['total','isROI','isRecoed_0','isRecoed_1','isRecoed_2','isRecoed_3',
                     'allBJetsIsRecoed','0_JetsNotIsRecoed','1_JetsNotIsRecoed','2_JetsNotIsRecoed','3_JetsNotIsRecoed','4_JetsNotIsRecoed',
                     '0_LeadNotIsRecoed','1_LeadNotIsRecoed','2_LeadNotIsRecoed',
                     '0_SubLeadNotIsRecoed','1_SubLeadNotIsRecoed','2_SubLeadNotIsRecoed',
                     'H1NotIsRecoed','H2NotIsRecoed',
                     'hasQuad',
                     'hasAllGenMatchedBJetsInQuad','hasH1GenMatchedInQuad','hasH2GenMatchedInQuad','hasAllJetsMatchedProperly',
                     'hasAllJetsMatchedOneToOneProperly','hasEachHMatchedProperly'

                     ]
 #       binsToPrint=["*"]
        histStore={}
        fileStore={}
        for tag in fileDict:
            print("="*40)
            print(tag,fileDict[tag].split("/")[-1])
            print("File : ",fileDict[tag])
            fileStore=ROOT.TFile(fileDict[tag],'READ')
            histStore=getTheObjectsFromFile(fileStore)
            nHist=histStore['tagsDumper']['trees']['sumEvts'].Clone()
            counts=getBinVsContent(nHist)
            #print(counts)
            tw=0
            t=0
            print("====== Entries    =====")
            if '*' in binsToPrint:
                for head in counts:    
                    print("\t",head, " : ",counts[head])
            else:
                for head in binsToPrint:
                    if head in counts:
                        print("\t",head, " : ",counts[head])
            nHist=histStore['tagsDumper']['trees']['sumWeighs'].Clone()
            counts=getBinVsContent(nHist)
            print("====== Weights    =====")
            if '*' in binsToPrint:
                for head in counts:    
                    print("\t",head, " : ",counts[head])
            else:
                for head in binsToPrint:
                    if head in counts:
                        print("\t",head, " : ",counts[head])
            nHist=histStore['tagsDumper']['trees']['sumWeighs'].Clone()
            counts=getBinVsContent(nHist)
            
            
            fileStore.Close()
            print()
            print()

