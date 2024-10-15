import glob
import sys,os,argparse


"""
python3 makeHaddFiles.py --force -d ./mergedFiles/  -b sig/         --execute
python3 makeHaddFiles.py --force -d ./mergedFiles/  -b singleH/     --execute
python3 makeHaddFiles.py --force -d ./mergedFiles/  -b doubleH/     --execute
python3 makeHaddFiles.py --force -d ./mergedFiles/  -b kappaScan/   --execute

python3 makeHaddFiles.py --doOnlyNominal --force -d ./mergedFilesNominal/ --execute -b sig/
python3 makeHaddFiles.py --doOnlyNominal --force -d ./mergedFilesNominal/ --execute -b singleH/
python3 makeHaddFiles.py --doOnlyNominal --force -d ./mergedFilesNominal/ --execute -b doubleH/
python3 makeHaddFiles.py --doOnlyNominal --force -d ./mergedFilesNominal/ --execute -b kappaScan/

"""
parser = argparse.ArgumentParser()
parser.add_argument('-b',"--base" , help="base path" , default='./')
parser.add_argument('-t',"--tag" , help="Tag" , default='merged')
parser.add_argument('-d',"--dest" , help="destination path" , default='./')
parser.add_argument("--doOnlyNominal", help="had only nominal files",default=False,action='store_true')
parser.add_argument("--execute", help="execute hadd",default=False,action='store_true')
parser.add_argument("--force", help="execute hadd with forced rewrite",default=False,action='store_true')
args = parser.parse_args()

base=args.base
searchPath=f"{base}/*/*/*/*.root"
print(searchPath)
allFiles= glob.glob(searchPath)
#allFiles= glob.glob(f"kappaLambda_scan/c3_19_d4_19/*/*/*.root")
tag=args.tag

if not os.path.exists(args.dest):
    cmd=f"mkdir -p {args.dest}"
    print(f"[]$ {cmd} ")
    os.system(cmd)

if tag[:-1]!='_':
    tag=tag+'_'

fileMap={}

for fls in allFiles:
    items=fls.split('/')
    proc,syst,year=items[-4],items[-3],items[-2]
    if proc not in fileMap:
        print(f"Adding Process : {proc}")
        fileMap[proc]={}
    if year not in fileMap[proc]:
        print(f"Adding year : {proc}/{year}")
        fileMap[proc][year]={}
    if syst not in fileMap[proc][year]:
        fileMap[proc][year][syst]=[]
    fileMap[proc][year][syst].append(fls)
i=0
for proc in fileMap:
    for year in fileMap[proc]:
        print(f"Processing {proc}/{year}")
        foutName=f"{args.dest}/{tag}{proc}_{year}_13TeV.root"
        cmd=f"hadd {foutName} "
        if args.force:
            cmd=cmd.replace("hadd ","hadd -f ")
        nFiles=0
        for syst in fileMap[proc][year]:
            if args.doOnlyNominal and ('nominal' not in syst):
                continue
            #print(f"\r processing {syst}      ",end="    \t\t")
            for fls in fileMap[proc][year][syst]:
                cmd+=f" {fls}"
                nFiles+=1
        print(f"Hadding {nFiles} files into {foutName}")
        if args.execute:
            os.system(cmd)
