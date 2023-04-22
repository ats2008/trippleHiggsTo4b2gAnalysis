import matplotlib.pyplot as plt
import mplhep as hep
import uproot as urt
import ROOT
import json,sys,os,argparse
import numpy as np
import Util as utl
import pickle,flatdict

import scaleFactorUtil as scl

import hep_ml as hepml


hep.style.use("CMS")
blind='CMS_hgg_mass < 115 || CMS_hgg_mass > 135.0'

def main():
 
    parser = argparse.ArgumentParser(prog='testTrainSplitter.py',description='For all the years in filelist , and all the tags inside , the files are split into test and train based on eventID%STRIDE==0 or not [ STRIDE = nEVTS*_frac]')
    parser.add_argument('-v',"--version", help="Version of the specific derivation ",default='')
    parser.add_argument('-i',"--inputFile", help="Input File", required=True)
    parser.add_argument('-o',"--outputDestination", help="output Destination", required=True)
    parser.add_argument('-t',"--testFraction", help="fraction of events for test", type=float ,default=0.25)
    parser.add_argument("--slim", help="slim the files to have less events", type=int ,default=-1)
    args = parser.parse_args()
    version = args.version
    inputFile  = args.inputFile
    outputDestination  = args.outputDestination
    saveOutput = True
    testFraction=args.testFraction
    nTestStride= int(1/testFraction  )

    fileDict={}
    with open(inputFile) as f:
        fileDict=json.load(f)
        fileDict_={}
        for year in fileDict:
            fileDict_[year]=flatdict.FlatDict(fileDict[year])
        fileDict=fileDict_
    yearsToProcess=['2018','2017','2016PreVFP','2016PostVFP','run2','2016']
    yearsToProcess=list(fileDict.keys())
    
    varToBinMap={}

    # We read the tree from the file and create a RDataFrame, a class that
    # allows us to interact with the data contained in the tree.
    trainBase=outputDestination+'/train'
    testBase =outputDestination+'/test'
    if args.slim > 0:
        trainBase=outputDestination+'/slim_train'
        testBase =outputDestination+'/slim_test'
        
    os.makedirs(outputDestination,exist_ok=True) 
    os.makedirs(trainBase,exist_ok=True) 
    os.makedirs(testBase,exist_ok=True) 
    print("Test Fraction : " , testFraction)
    print("Test Stride   : " , nTestStride)
    print("Train files to be stored in : ",trainBase)
    print("Test files to be stored in : " ,testBase)

    rdataFrames={}
    nMax=0
    for yr  in yearsToProcess:
        for fTag in fileDict[yr]:
         nMax+=1
    n=0
    for yr  in yearsToProcess:
        for fTag in fileDict[yr]:
            fileName=fileDict[yr][fTag]
            treeName = "trees/bkg_13TeV_TrippleHTag_0"
            if 'sig' in fileName.split('/')[-1]:
                treeName = "trees/ggHHH_125_13TeV"
            elif 'data' in fileName.split('/')[-1] :
                treeName = "trees/Data_13TeV_TrippleHTag_0"
            dset=fTag.split(":")[-1]
            n+=1
            print(" [",n,"/",nMax,"] Processing dataset : ",dset," from ",yr,"with tree",treeName)
            rdataFrame= ROOT.RDataFrame(treeName, fileName)

            train_dset=rdataFrame.Filter( 'rdfentry_%'+str(nTestStride)+'!=0' ) 
            test_dset =rdataFrame.Filter( 'rdfentry_%'+str(nTestStride)+'==0' )    
            if args.slim > 0:
                train_dset=train_dset.Filter( 'rdfentry_ < '+str(args.slim) ) 
                test_dset =test_dset.Filter(  'rdfentry_ < '+str(args.slim) )    
                

            trainOFileName=trainBase+'/'+fileName.split('/')[-1]
            testOFileName=testBase+'/'+fileName.split('/')[-1]
            print("\t Total Number of events : ",rdataFrame.Count().GetValue())
            print("\t       Train Events     : ",train_dset.Count().GetValue())
            print("\t       Test Events      : ",test_dset.Count().GetValue())
            if os.path.exists(trainOFileName):
                trainOFileName=trainOFileName.replace('.root','_train.root')
            print("\t writing train file : ",trainOFileName," [ tree : ", treeName ," ]")
            train_dset.Snapshot(treeName,trainOFileName)
            if os.path.exists(testOFileName):
                testOFileName=testOFileName.replace('.root','_test.root')
            print("\t writing test  file : ",testOFileName)
            test_dset.Snapshot(treeName ,testOFileName)
            print("\t DONE ! \n")
if __name__=='__main__':
    main()
