import glob,json
import sys,os,argparse
import uproot as urt
import numpy as np
import tabulate
usage="""
DONOM=" -n "
DOSB=" --doSB "
DOSR=" --doSR "
python3 makeYields.py -l -b singleH/ -t singleH -d summary/            $DONOM    $DOSB  $DOSR 
python3 makeYields.py -l -b doubleH/ -t doubleH -d summary/            $DONOM    $DOSB  $DOSR 
python3 makeYields.py -l -b sig/ -t ggHHH -d summary/ -s 100           $DONOM    $DOSB  $DOSR 
python3 makeYields.py -l -b kappaScan/ -t kappaScan -d summary/ -s 100 $DONOM    $DOSB  $DOSR 
python3 makeYields.py -l -b bkg/ -t bkg -d summary/                    $DONOM    $DOSB  $DOSR
python3 makeYields.py -l -b ttX/ -t ttX -d summary/                    $DONOM    $DOSB  $DOSR
python3 makeYields.py -l -b data/ -t data -d summary/ --isData         $DONOM    $DOSB  $DOSR
"""
parser = argparse.ArgumentParser()
parser.add_argument('-b',"--base" , help="base path" , default='./')
parser.add_argument('-t',"--tag" , help="Tag" , default='merged')
parser.add_argument('-d',"--dest" , help="destination path" , default='./')
parser.add_argument("-l","--doLumiWeights", help="Do luminosity weighting for yields",default=False,action='store_true')
parser.add_argument("-n","--nominalOnly", help="get only nimnal yield",default=False,action='store_true')
parser.add_argument(    "--doSR", help="get yield from SR",default=False,action='store_true')
parser.add_argument(    "--doSB", help="get yield from SB",default=False,action='store_true')
parser.add_argument(    "--isData", help="This data",default=False,action='store_true')
parser.add_argument("-s","--scale", help="scale the yields by a number",default=1.0,type=int)
args = parser.parse_args()

lumiMap={'2018':58,'UL18':58.0,'2017':41.5,'2016':36.3 ,'RR17':41.5,'RR18':48.0,'RR16':19.5+16.8, '2016PreVFP':19.5,'UL16Post':16.8,
        'UL2016Pre':19.5,'UL16Pre':19.5,'UL2016Post':16.8,'2016PostVFP':16.8,'run2':137.61,'UL17':41.5,'UL2017':41.5,'UL2018':58.0}

base=args.base
allFiles= glob.glob(f"{base}/*/*/*/*.root")
tag=args.tag
scl=args.scale
if not os.path.exists(args.dest):
    cmd=f"mkdir -p {args.dest}"
    print(f"[]$ {cmd} ")
    os.system(cmd)

if tag[:-1]!='_':
    tag=tag+'_'
if args.doSR:
    tag=tag+'SR_'
if args.doSB:
    tag=tag+'SB_'


fileMap={}

fileMapForYield={}
for fls in allFiles:
    items=fls.split('/')
    proc,syst,year=items[-4],items[-3],items[-2]
    if proc=='sig':
        proc='ggHHH'
    cat=[ i for  i in fls.split('/')[-1].split('_') if 'CAT'  in i]
    if len(cat)==1:
        cat=cat[0]
    else:
        print("CAT not resolved for  : ",fls)
        exit(1)
        continue
    if cat not in fileMapForYield:
        print(f"Adding Cat : {cat}")
        fileMapForYield[cat]={}
    if proc not in fileMapForYield[cat]:
#         print(f"Adding Process : {proc}")
        fileMapForYield[cat][proc]={}
    if syst not in fileMapForYield[cat][proc]:
        fileMapForYield[cat][proc][syst]={}
    if year not in fileMapForYield[cat][proc][syst]:
#         print(f"Adding year : {proc}/{year}")
        fileMapForYield[cat][proc][syst][year]={}
    fileMapForYield[cat][proc][syst][year]=fls
allTables={}
allData={}
summaryData={}
summaryTable={}
allYears=['RR16','RR17','RR18','U1L6']
i0=1
for cat in fileMapForYield:
    print(f"{i0}/{len(fileMapForYield)} Doing cat {cat} ") ;i0+=1
    allTables[cat]={}
    summaryTable[cat]={'data':[],'year':[]}
    i1=1
    for proc in fileMapForYield[cat]:
        print(f"\t {i1}/{len(fileMapForYield[cat])} Doing proc {proc} ") ;i1+=1
        years=[]
        allTables[cat][proc]={'data':[],'year':[]}
        N=len(fileMapForYield[cat][proc])
        i=0
        for syst in fileMapForYield[cat][proc]:
            if proc not in allData:
                allData[proc]={}
            if cat not in allData[proc]:
                allData[proc][cat]={}
            if syst not in allData[proc][cat]:
                allData[proc][cat][syst]={}
            if syst=='nominal':
                print(f"\t\t {i}/{N} Doing syst {syst} ",end="\t\t\t") ;i+=1
                w={}
                wl={}
                wl2={}
                years = list( fileMapForYield[cat][proc][syst].keys() )
                for yr in years:
                    w[yr]=0.0
                    wl[yr]=0.0
                    wl2[yr]=0.0
                    if yr in fileMapForYield[cat][proc][syst]:
                        with urt.open(fileMapForYield[cat][proc][syst][yr]) as f:
                            kk=list(f['trees'].keys())[0]
                            dt_=f['trees'][kk].arrays(['weight','lumi','CMS_hgg_mass'],library='np')
                            if len(dt_['weight']) < 1:
                                dt_={'weight':np.array([0.0]), 'lumi' : np.array([0.0]),'CMS_hgg_mass' :np.array([0.0])}
                            if args.doSR:
                                mask=np.logical_and( dt_['CMS_hgg_mass'] >=115.0 ,  dt_['CMS_hgg_mass'] <= 135.0 )
                                for ky in dt_:
                                    dt_[ky]=dt_[ky][mask]
                            elif args.doSB:
                                mask=np.logical_or(  dt_['CMS_hgg_mass'] < 115.0 ,  dt_['CMS_hgg_mass'] > 135.0 )
                                for ky in dt_:
                                    dt_[ky]=dt_[ky][mask]

                            if len(dt_['weight']) < 1:
                                dt_={'weight':[0.0], 'lumi' : [0.0] }
                            lmi    = np.unique(dt_['lumi'])[0]
                            lmi    = lumiMap[yr]
                            if args.isData:
                                lmi=1.0
                            #print("\n",yr,lmi,lumiMap[yr])
                            w[yr]  = np.sum(dt_['weight'])
                            wl[yr] = np.sum(dt_['weight'])*lmi
                            wl2[yr] = np.sqrt(np.sum(np.array(dt_['weight'])*np.array(dt_['weight'])))*lmi
                
                if not args.doLumiWeights:
                    allTables[cat][proc]['data'].append( [syst] +[ f"{w[i]*1e3*scl:.6f}" for i in years ] )
                    summaryTable[cat]['data'].append( [proc] +[ " | ".join([ f"{i} : {w[i]*1e3:.6f}" for i in years ]) ] )
                else:
                    allTables[cat][proc]['data'].append( [syst] +[ f"{wl[i]*scl:.6f}" for i in years ] )
                    summaryTable[cat]['data'].append( [proc] +[ " | ".join([ f"{i} : {wl[i]*scl:.6f}" for i in years ]) ] )
                wNominal=w
                wNominalL=wl
                for yr in years:
                    allData[proc][cat][syst][yr]=[float(wl[yr]),float(wl2[yr])]
                break
        #for syst in fileMapForYield[cat][proc]:
        #    if args.nominalOnly:
        #        break
            else:
                print(f"\r\t\t {i}/{N} Doing syst {syst} ",end="                          ") ;i+=1
                if syst=='nominal':
                    continue
                w={}
                for yr in years:
                    w[yr]=0.0
                    wl[yr]=0.0
                    wl2[yr]=0.0
                    if yr in fileMapForYield[cat][proc][syst]:
                        with urt.open(fileMapForYield[cat][proc][syst][yr]) as f:
                            kk=list(f['trees'].keys())[0]
                            dt_=f['trees'][kk].arrays(['weight','lumi','CMS_hgg_mass'],library='np')
                            if len(dt_['weight']) < 1:
                                dt_={'weight':[0.0], 'lumi' : [0.0] }
                            w[yr]=np.sum(dt_['weight'])
                            lmi=np.unique(dt_['lumi'])[0]
                            wl[yr]=np.sum(dt_['weight'])*lmi*scl
                            wl2[yr] = np.sqrt(np.sum(np.array(dt_['weight'])*np.array(dt_['weight'])))*lmi*scl

                if not args.doLumiWeights:
                    allTables[cat][proc]['data'].append( [syst] +[ f"{w[i]*scl:.6f}" for i in years ] )
                else:
                    allTables[cat][proc]['data'].append( [syst] +[ f"{wl[i]*scl:.6f}" for i in years ] )
                for yr in years:
                    allData[proc][cat][syst][yr]=[float(wl[yr]),float(wl2[yr])]
        print()
        allTables[cat][proc]['year']=years

if not args.nominalOnly:
    for cat in allTables:
        for proc in allTables[cat]:
            p=tabulate.tabulate(allTables[cat][proc]['data'],headers=["SYST"] + allTables[cat][proc]['year'],tablefmt='mixed_grid')
            foutName=f"{args.dest}/summary_{tag}{proc}_{cat}.txt"
            with open(foutName,'w') as f:
                f.write(p)

foutName=f"{args.dest}/{tag}summary.txt"
if args.doLumiWeights:
    foutName=f"{args.dest}/{tag}summary_lumiScaled.txt"
print(f"Summary File saved as {foutName}")
with open(foutName,'w') as f:
    f.write(f"\n\n-----------    >  Yields here are scaled by {scl}  <    -----------\n")
    for cat in summaryTable:
        p=tabulate.tabulate(summaryTable[cat]['data'],headers=["PROC",'Yields'] ,tablefmt='mixed_grid')
        f.write(f"\n ========  CATEGORY : {cat} ========\n")
        f.write(p)
        f.write("\n\n")
        
foutName=f"{args.dest}/{tag}yields.json"
with open(foutName,'w') as f:
    print("Exporting yields to : ",foutName)
    json.dump(allData,f,indent=4)



