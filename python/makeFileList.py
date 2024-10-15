import glob,json,os

search_path="*/*.root"
foutname='_flist_out.json'

flist=glob.glob(search_path)

dstore={
  'run2':{
           'singleH' : {},
           'doubleH' : {},
           'data' : {},
           'kappaScan' : {},
           'bkg' : {},
           'sig' : {}
        }
  }

for fl in flist:
    pth=os.path.abspath(fl)
    ky=fl.split("/")[-1].replace(".root","").replace("mc_","")
    tag='bkg'
    if 'HHH' in ky:
        tag='sig'
    elif 'HH' in ky:
        tag='doubleH'
    elif 'H' in ky:
        tag='singleH'
    elif 'data' in ky:
        tag='data'
    elif 'd4' in ky:
        tag='kappaScan'

    print(tag,"\n  ->",ky," : \t> ",pth)
    dstore['run2'][tag][ky]=pth

with open(foutname,'w') as f:
    json.dump(dstore,f,indent=4)
    

