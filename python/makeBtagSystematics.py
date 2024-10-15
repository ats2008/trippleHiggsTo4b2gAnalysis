import uproot as urt
import json,sys,os,argparse
import numpy as np
usage="""

python3 slimExportedFile.py  -i merged_c3_m1p5_d4_m0p5_UL2017_13TeV.root --proc c3_m1p5_d4_m0p5 --year UL2017

"""
BTAG_DEF_SCALE="btag_scale"
allBtagSFSyst=[
        "btag_scale_down_hf","btag_scale_down_hfstats1",
        "btag_scale_down_hfstats2","btag_scale_down_jes","btag_scale_down_lf",
        "btag_scale_down_lfstats1","btag_scale_down_lfstats2","btag_scale_up_hf",
        "btag_scale_up_hfstats1","btag_scale_up_hfstats2","btag_scale_up_jes",
        "btag_scale_up_lf","btag_scale_up_lfstats1","btag_scale_up_lfstats2"
     ]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--inputFile", help="Input file" )
    parser.add_argument("-p","--proc", help="Process",default='proc' )
    parser.add_argument("-y","--year", help="Process",default='run2' )
    #parser.add_argument("-n","--numberToexport" , help="Number Of events to be exported !! ",default=None ,type=int  )
    args = parser.parse_args()
    
    print(f"Opening file {args.inputFile} ")
    systErrDict={}   
    fileIn=urt.open(args.inputFile)
    treeList=fileIn['trees'].keys()
    systErrDict["METADATA"]={
        "fileName" : args.inputFile
    }
    for i,ky in enumerate(treeList):
        if '13TeV' not in ky:
            continue
        
        dataIn=fileIn['trees'][ky].arrays(library='np')
        n_orig=len(dataIn['CMS_hgg_mass'])
        w=np.sum(dataIn['weight'])
        def_w=dataIn['weight']/dataIn['btag_scale']
        print(f" > Processing {ky}")
        systErrDict[ky]={'WEIGHT' : float(w) }
        for syst in allBtagSFSyst:
            w_sys=np.sum(def_w*dataIn[syst])
            systErrDict[ky][syst]=float(w_sys)
            print(f"\r   > {syst} : {w_sys/w:.3f}",end="                               ")
        print()
    foutname=f'btag_{args.proc}_{args.year}.json'
    print("Exported to : ",foutname)
    with open(foutname,'w') as     f:  
        json.dump(systErrDict,f,indent=4)
   

if __name__=='__main__':
    main( )

