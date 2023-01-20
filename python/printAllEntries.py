from __future__ import print_function
from Util import *
import os
from collections import OrderedDict
import argparse

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
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',"--file", default='out_data_data_2018_2018D_3p0_0.root',help="file to read" )
    args = parser.parse_args()
    
    for preReweight in [False]:
        fileDict=OrderedDict()
        fileDict[ "ggHHH"          ] =   args.file
        print(" -- > file : ",fileDict[ "ggHHH"          ])          
        binsToPrint=['total','isROI','isRecoed_0','isRecoed_1','isRecoed_2','isRecoed_3',
                     'allBJetsIsRecoed','0_JetsNotIsRecoed','1_JetsNotIsRecoed','2_JetsNotIsRecoed','3_JetsNotIsRecoed','4_JetsNotIsRecoed',
                     '0_LeadNotIsRecoed','1_LeadNotIsRecoed','2_LeadNotIsRecoed',
                     '0_SubLeadNotIsRecoed','1_SubLeadNotIsRecoed','2_SubLeadNotIsRecoed',
                     'H1NotIsRecoed','H2NotIsRecoed',
                     'hasQuad', 'hasMatchInQuad_0', 'hasMatchInQuad_1','hasMatchInQuad_2','hasMatchInQuad_3','hasMatchInQuad_4',
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
            #nHist=histStore['tagsDumper']['trees']['sumWeighs'].Clone()
            #counts=getBinVsContent(nHist)
            # print("====== Weights    =====")
            # if '*' in binsToPrint:
            #     for head in counts:    
            #         print("\t",head, " : ",counts[head])
            # else:
            #     for head in binsToPrint:
            #         if head in counts:
            #             print("\t",head, " : ",counts[head])
            # nHist=histStore['tagsDumper']['trees']['sumWeighs'].Clone()
            # counts=getBinVsContent(nHist)
            
            
            fileStore.Close()
            print()
            print()

