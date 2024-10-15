import uproot as urt
import ROOT
import json,sys,os,argparse
import numpy as np
import prettytable as ptab
usage="""

python3 displayFile.py -i merged_c3_m1p5_d4_m0p5_UL2017_13TeV.root 

"""

scf = 0.00017

def exportTreeForFit(final_dataset,treeNameToExport,export_filename): 
    export_file = urt.recreate(export_filename)
    #print("\n\t  file exported  : ",export_filename)
    allBranchesToWrite={ ky : final_dataset[ky].dtype for ky in final_dataset}
    export_file.mktree( treeNameToExport ,allBranchesToWrite)
    export_file[treeNameToExport].extend(final_dataset)
    export_file.close()

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--inputFile", help="Input file" )
    parser.add_argument("-n","--onlyNominal", help="Print only nominal",action='store_true'  )
    args = parser.parse_args()
    print(f"Opening file {args.inputFile} ")
    
    fileIn=urt.open(args.inputFile)
    treeList=fileIn['trees'].keys()
    for i,ky in enumerate(treeList):
        if args.onlyNominal:
            if 'sigma' in ky:
                continue
        kyV=ky.split(';')[0]
        dataIn=fileIn['trees'][ky].arrays(library='np')
        n_orig=len(dataIn['CMS_hgg_mass'])
        w=np.sum(dataIn['weight'])
        print(f"\tProcessing {kyV} with {n_orig} events w : {w} ",end="\n")
    
    print()

if __name__=='__main__':
    main( )

