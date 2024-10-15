usage="""

# SET THE YEARLIST and FILELIST  
python3 ttHCorrection.py

"""
lumiMap={'2018':58,'2017':41.5,'2016':36.3 , '2016PreVFP':19.5,'2016PostVFP':16.8,'run2':137.61}


import uproot as urt
import ROOT
import json,sys,os,argparse
import numpy as np
import prettytable as ptab


def exportTreeForFit(final_dataset,treeNameToExport,export_filename): 
    export_file = urt.recreate(export_filename)
    #print("\n\t  file exported  : ",export_filename)
    allBranchesToWrite={ ky : final_dataset[ky].dtype for ky in final_dataset}
    export_file.mktree( treeNameToExport ,allBranchesToWrite)
    export_file[treeNameToExport].extend(final_dataset)
    export_file.close()

def main():

    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    BASE='/tmp/atsFilesForSkimm/'
    if os.path.exists(BASE):
        cmd=f'rm -rf {BASE}'
        print(cmd);os.system(cmd)

    tthFileList=["merged_ttH_UL16Post_13TeV.root","merged_ttH_UL16Pre_13TeV.root","merged_ttH_UL17_13TeV.root","merged_ttH_UL18_13TeV.root"]
    tthFileList=["mergedFiles/"+ i for i in tthFileList]
    yrMap      =["2016PostVFP","2016PreVFP","2017","2018"]
    print(f"Opening file {tthFileList} ")
    
    fileList=[]
    
    normalYield={}
    for fl in tthFileList:
        fileIn=urt.open(fl)
        treeList=fileIn['trees'].keys()
        for i,ky in enumerate(treeList):
            if 'sigma' in ky:
                continue
            cat=ky.split(";")[0].split('_')[-1]
            if cat not in normalYield:
                normalYield[cat]=0.0
            dataIn=fileIn['trees'][ky].arrays(library='np')
            n_orig=len(dataIn['CMS_hgg_mass'])
            w=np.sum(dataIn['weight'])
            print(" processed ",ky," from ",fl," with yield ",w)
            normalYield[cat]+=w
        fileIn.close()
    for ky in normalYield:
        print(f" yield in  {ky} : ",normalYield[ky])
    for fl,yr in zip(tthFileList,yrMap):
        cmd='mkdir -p '+BASE
        print(cmd);    os.system(cmd)
        fileIn=urt.open(fl)
        treeList=fileIn['trees'].keys()
        W_NORM={}
        for i,ky in enumerate(treeList):
            if 'sigma' not in ky:
                cat=ky.split(";")[0].split('_')[-1]
                dataIn=fileIn['trees'][ky].arrays(library='np')
                W_NORM[cat]=np.sum(dataIn['weight'])
                print("Setting W_NORM as ",W_NORM)
        if not W_NORM:
            print("DID NOT GET NORM !! for ",fl, " in cat ",cat )
            exit()
        print(" W_NORM  = ",W_NORM)
        for i,ky in enumerate(treeList):
            kyV=ky.split(';')[0]
            dataIn=fileIn['trees'][ky].arrays(library='np')
            dataOut=dataIn
            if 'sigma' not in ky:
                cat=ky.split(";")[0].split('_')[-1]
                w=normalYield[cat]*lumiMap[yr]/lumiMap['run2']
                dataOut={ky : dataIn[ky] for ky in dataIn}
                dataOut['weight']=dataOut['weight']*0.0 + w/len(dataOut['weight'])
                w_out=np.sum(dataOut['weight'])
                print(f"\tNominal : Processing {kyV} to a corrected weight {normalYield[cat]}*{lumiMap[yr]}/{lumiMap['run2']} = {w_out=} ",end="\n")
            else:
                cat=ky.split(";")[0].split('_')[-2]
                w=normalYield[cat]*lumiMap[yr]/lumiMap['run2']
                dataOut={ky : dataIn[ky] for ky in dataIn}
                w_orig=np.sum(dataOut['weight'])
                dataOut['weight']=dataOut['weight']*0.0 + w*(w_orig/W_NORM[cat])/ len(dataOut['weight']) 
                w_out=np.sum(dataOut['weight'])
                print(f"\t   syst :Processing {kyV} to a corrected weight {w}*{w_orig}/{W_NORM[cat]} = {w_out=} ",end="\n")

            export_filename=BASE+f'/{kyV}_{i}.root' 
            treeNameToExport=f"trees/{kyV}"
            #print("w = ",w,f" {n_orig= } {n=} {w_out} ")
            #print("Exporting to tmp file name : ",export_filename)
            exportTreeForFit(dataOut,treeNameToExport,export_filename)
        finName=fl.split('/')[-1]
        fout=finName.replace('.root','_normUpdated.root')
        cmd=f'hadd -f -v 1 {fout} {BASE}/*.root'
        print(cmd);os.system(cmd)
        cmd=f'rm -rf {BASE}'
        print(cmd);os.system(cmd)

if __name__=='__main__':
    main( )

