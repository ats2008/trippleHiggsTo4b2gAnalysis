import os,argparse
"""
  []$  ls singleH/* doubleH/* sig/* -d > fls
       python3 makeSystematicsSplit.py  -i fls -t v5 -j 1
"""
parser = argparse.ArgumentParser()
parser.add_argument('-i',"--inputfile" , help="input file" , default='fls')
parser.add_argument('-t',"--tag" , help="tag" , default='')
parser.add_argument('-j',"--filesPerJob" , help="number hadd cmds per job" , default=5, type=int)
args = parser.parse_args()

condor_base=f'_Condor_{args.tag}'
fileInName=args.inputfile

if not os.path.exists(condor_base):
    cmd= f'mkdir -p {condor_base}'; print(" [] $ ",cmd)
    os.system(cmd)

flsList=[]
with open(fileInName) as f:
    flsList=f.readlines()
i=1
deleteFileCMDs=[]
nFLSPerJob=args.filesPerJob
flsCount=0
subScript="""
executable = $(filename)
request_cpus = 3
output = $Fp(filename)run.$(Cluster).$(ProcId).stdout
error = $Fp(filename)run.$(Cluster).$(ProcId).stderr
log = $Fp(filename)run.$(Cluster).$(ProcId).log
queue filename matching (@@QUEUECMD)

"""
header="""#!/bin/bash
source /grid_mnt/t3home/athachay/.bashrc
source /cvmfs/cms.cern.ch/cmsset_default.sh 
set -x
export HOME=/home/athachay
export X509_USER_PROXY=/home/athachay/.proxy/x509up_u56621
"""
condor_base=os.path.abspath(condor_base)
subScript=subScript.replace("@@QUEUECMD",f"{condor_base}/*/sub_*.sh")
f=None
summary=[]
for fl in flsList:  
    fl=fl[:-1]
    print("Processing ",f"{i} / {len(flsList)}  | " ,fl)
    i+=1
    flsCount+=1
    if flsCount> nFLSPerJob or (not f):
        if f:
            f.close()

        cdir=f'{condor_base}/Job_{i}' ; cdir=os.path.abspath(cdir)
        if not os.path.exists(cdir):
            cmd= f'mkdir -p {cdir}'; print(" [] $ ",cmd)
            os.system(cmd)
        f=open(f"{condor_base}/Job_{i}/sub_{i}.sh","w")
        f.write(header)
        flsCount=1
    #sig_syst_vH_UL17_MaterialForwardDown01sigma_v26p1
    tag=fl.split("/")[-2] #items[0]
    items=fl.replace("doubleH_private_","").replace("singleH_private_","").split('/')[-1].replace("doubleH_","").split('_')
    if 'private' in items:
        items.remove('private')
    print(fl,items)
    era=items[-2].replace("private","")
    proc='_'.join(items[0:1])
    if 'ggHH_kl' in fl:
        proc='_'.join(items[0:2])
    # kappaScan/kappaScan_UL2018_c3_m1p5_d4_m0p5_nominal
    if ('c3' in fl) and ( 'd4' in fl):   
        proc='_'.join(items[2:6])
        era=items[1]
    syst=items[-1]
    ver='v33p0'
    detail=f"{tag=} {proc=} {era=} {syst=} {ver=}"
    print(detail)
    summary.append(detail)
    #exit(1)
    # doubleH_ggZTo2BHHTo2B2G_privateUL17_FNUFEBDown01sigma_v26p1
    # items=fl.split('/')[-1].split('_')
    # tag=items[1]
    # era=items[2]
    # syst=items[3]
    # ver=items[4]
    
    pth=f'{ver}/{tag}/{proc}/{syst}/{era}/' ; pth=os.path.abspath(pth)
    if not os.path.exists(pth):
        cmd= f'mkdir -p {pth}'; print(" [] $ ",cmd)
        os.system(cmd)
    absFl=os.path.abspath(fl)
    cmd= f'hadd -f  {pth}/{tag}_{era}_{syst}.root {absFl}/*.root'; #print(" [] $ ",cmd)
    f.write(cmd+'\n')
#     os.system(cmd)
    cmd=f"rm -r {fl}" ; deleteFileCMDs.append(cmd)
        
#    if i > 5:
#        break
        
with open('cleanup.sh','w') as f:
    i=0
    for l in deleteFileCMDs:
        i+=1
        f.write(f"echo {i}/{len(deleteFileCMDs)} , {l}\n")
        f.write(l+"\n")
os.system('chmod +x cleanup.sh')  
os.system(f'chmod +x {condor_base}/*/*.sh')    
cfname=f'{condor_base}/job_sub.sub'
with open(cfname,'w') as f:
    f.write(subScript)

print()    
print()    
for i,su in enumerate(summary):
    print(summary[i])    
print(f"csub {cfname}")
