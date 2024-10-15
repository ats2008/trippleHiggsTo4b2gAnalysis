import uproot,argparse,os
import awkward as awk
import pickle as pkl

iData=['run',
     'luminosityBlock',
     'PV_npvs',
     'L1PreFiringWeight_Nom',
     'L1PreFiringWeight_Up',
     'L1PreFiringWeight_Dn',
     'event',
     'LHEWeight_originalXWGTUP',
     'LHEScaleWeight', # 9 Nos  QCDScale
     'LHEPdfWeight', # 102
     'PSWeight' # 4 , Parton Scale weight
]

parser = argparse.ArgumentParser()
parser.add_argument('-i',"--inputFile", help="Input File",default=None)
parser.add_argument('-t',"--tag", help="output tag",default=None)
args = parser.parse_args()
ipf=args.inputFile
if ipf.startswith('/store'):
    ipf="root://cms-xrd-global.cern.ch/"+ipf
print("Getting file ",ipf)
#os.system(f'xrdcp {ipf} .')
#exit()
iFile=args.inputFile.split('/')[-1]
iFile=ipf
file=uproot.open(iFile,xrootd_handler=uproot.source.xrootd.XRootDSource)
print(f"File {iFile} opened sucessfully !")
Events=file['Events']
data=Events.arrays(iData)
print("Extracted ",len(data)," events ")
ofname=f"output_{args.tag}.root"
ofile = uproot.recreate(ofname)
dDict={ky :data[ky].to_numpy() for ky in data.fields}
for ky in dDict:
    print(ky,dDict[ky].shape,data[ky].ndim)
ofile["Events"]=dDict
ofile.close()
#os.system('rm -r '+iFile)

