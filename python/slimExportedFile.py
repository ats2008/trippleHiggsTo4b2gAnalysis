import uproot as urt
import ROOT
import json,sys,os,argparse,uuid,time
import numpy as np
import prettytable as ptab
usage="""

 python3 slimExportedFile.py  -i  mergedFiles/merged_sig_UL16Post_13TeV*.root -f 0.022
 python3 slimExportedFile.py  -i  mergedFiles/merged_sig_UL17_13TeV*.root -f 0.015
 python3 slimExportedFile.py  -i  mergedFiles/merged_ggHHH_UL18_13TeV*.root -f 0.015
 

 python3 slimExportedFile.py  -i  mergedFiles/merged_sig_UL16Post_13TeV.root -f 0.022
 python3 slimExportedFile.py  -i  mergedFiles/merged_sig_UL17_13TeV.root -f 0.015
 python3 slimExportedFile.py  -i  mergedFiles/merged_ggHHH_UL18_13TeV.root -f 0.015
 python3 slimExportedFile.py  -i merged_c3_m1p5_d4_m0p5_UL2017_13TeV.root -f 0.5

"""

branches_to_export=['run','event','lumi','year','weight','CMS_hgg_mass',
#                    'h1bb_mass','h2bb_mass','quadjet_0_deepJetScore',
#                    'quadjet_1_deepJetScore','quadjet_2_deepJetScore','quadjet_3_deepJetScore',
                     'mva_vsNonResonantBDT_v24p0_1','mva_vsPeakingBDT_vsTTH_v34p0_1','dZ']

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
    parser.add_argument("-d","--prefix", help="Output destination" , default='./' )
    parser.add_argument("-p","--printStats", help="print statistics per cat", action='store_true' )
    parser.add_argument("-b","--branch_only", help="do only branches slim", action='store_true' )
    parser.add_argument("-f","--fraction" , help="Fraction of inputr to be exported !! ",default='1.0'   )
    parser.add_argument("-n","--minNumberToexport" , help="minimum number Of events to be exported !! ",default=7000 ,type=int  )
    #parser.add_argument("-n","--numberToexport" , help="Number Of events to be exported !! ",default=None ,type=int  )
    args = parser.parse_args()
    specific=uuid.uuid4().hex.upper()[0:6]
    BASE=f'/tmp/atsFilesForSkimm/{specific}/'
    if os.path.exists(BASE):
        cmd=f'rm -rf {BASE}'
        print(cmd);os.system(cmd)
    
    cmd='mkdir -p '+BASE
    print(cmd)
    os.system(cmd)
    time.sleep(0.2)
    print(BASE)
    fracs={}
    if ',' not in args.fraction:
        fracs['all']=float(args.fraction)
    else:
        vals=args.fraction.split(",")
        for v in vals:
            ky,fr=v.split(":")
            fracs[ky]=max(float(fr),1.0)
    print("--> fraction dict : ",fracs)
    print(f"Opening file {args.inputFile} ")
    print(f"Sampling a fraction :  {args.fraction} ")
    fileIn=urt.open(args.inputFile)
    treeList=fileIn['trees'].keys()
    for i,ky in enumerate(treeList):
        kyV=ky.split(';')[0]
        if 'all' in fracs:
            fraction=fracs['all']
        else:
            for ct in fracs:
                if ct in ky:
                    fraction=fracs[ct]
                    break
            if fraction<0:
                print(f"No fraction found exiting ! [ {ky} ]")
                exit(1)

        if args.printStats:
            if 'sigma' in ky:
                continue
        dataIn=fileIn['trees'][ky].arrays(branches_to_export,library='np')
        n_orig=len(dataIn['CMS_hgg_mass'])
        w=np.sum(dataIn['weight'])
        n=int(fraction*n_orig+1)
        n=max(n,args.minNumberToexport)
        print(f"\tProcessing {kyV} with {n_orig} events (  -> {n} events .. frc : {fraction:.4f}) ",end="\n")
        if args.printStats:
            continue
        if not args.branch_only:
            rnd_placer = np.random.permutation(n_orig) 
            rnd_placer=rnd_placer[:n]
            dataOut={ky : dataIn[ky][rnd_placer] for ky in dataIn}
            dataOut['weight']=dataOut['weight']*w/np.sum(dataOut['weight'])
            w_out=np.sum(dataOut['weight'])
        else:
            dataOut=dataIn  
        export_filename=BASE+f'/{kyV}_{i}.root' 
        treeNameToExport=f"trees/{kyV}"
        #print("w = ",w,f" {n_orig= } {n=} {w_out} ")
        #print("Exporting to tmp file name : ",export_filename)
        exportTreeForFit(dataOut,treeNameToExport,export_filename)
    print()
    finName=args.inputFile.split('/')[-1]
    fout=args.prefix+'/'+finName.replace('.root','_skimmed.root')
    cmd=f'hadd -f {fout} {BASE}/*.root'
    print(cmd);os.system(cmd)
    cmd=f'rm -rf {BASE}'
    #print(cmd);os.system(cmd)


if __name__=='__main__':
    main( )

