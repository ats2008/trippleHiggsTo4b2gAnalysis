import json,argparse
import tabulate as tab
import numpy as np
import glob

sr_files=glob.glob("summary/*SR*.json")
sb_files=glob.glob("summary/*SB*.json")

parser = argparse.ArgumentParser()
parser.add_argument('-c',"--cat" , help="Tag" , default='CAT0')
parser.add_argument("--doCAT", help="Print cat wise table",default=False,action='store_true')
parser.add_argument("-l","--doLatex", help="Print Latex table",default=False,action='store_true')
args = parser.parse_args()

dataStore_sr={}
for fl in sr_files:
    with open(fl) as f:
        dta    =json.load(f)
        for tag in dta:
            dataStore_sr[tag]=dta[tag]
        
dataStore_sb={}
for fl in sb_files:
    with open(fl) as f:
        dta    =json.load(f)
        for tag in dta:
            dataStore_sb[tag]=dta[tag]


for proc in dataStore_sb.keys():
    for cat in dataStore_sb[proc].keys():
        for yr in dataStore_sb[proc][cat]['nominal']:
            dataStore_sb[proc][cat]['nominal'][yr]=np.array([dataStore_sb[proc][cat]['nominal'][yr][0],
                                                             dataStore_sb[proc][cat]['nominal'][yr][1]**2])  
for proc in dataStore_sr.keys():
    for cat in dataStore_sr[proc].keys():
        for yr in dataStore_sr[proc][cat]['nominal']:
            dataStore_sr[proc][cat]['nominal'][yr]=np.array([dataStore_sr[proc][cat]['nominal'][yr][0],
                                                             dataStore_sr[proc][cat]['nominal'][yr][1]**2]
                                                           )           

merger_maps =    {
        'TTX' : ['mc_ttGG', 'mc_ttGJ', 'mc_ttJJ'],
        'GGJ' : ['mc_ggBox', 'mc_ggBox1Bjet', 'mc_ggBox2Bjet'],
        'H'   : ['ggH', 'bbh', 'ttH', 'vbfH', 'vH'],
        'vHH' : ['ZToBBHHTo2B2G','WToQQHHTo2B2G'],
        'HH'  : ['ggHH', 'ttHH','ZToBBHHTo2B2G','WToQQHHTo2B2G'],
        'Expected' : ['morphedData','mc_ttGG', 
                      'mc_ttGJ', 'mc_ttJJ',
                      'mc_ggBox', 'mc_ggBox1Bjet', 'mc_ggBox2Bjet'
                     ],
        'ExpectedAll' : ['morphedData','mc_ttGG', 
                      'mc_ttGJ', 'mc_ttJJ',
                      'mc_ggBox', 'mc_ggBox1Bjet', 'mc_ggBox2Bjet','H','HH'
                     ],


    }




yieldStore={}
for tag,dstore in zip(["SR","SB"],[dataStore_sr,dataStore_sb]) :
    yields={}
    for i in range(5):
        cat=f"CAT{i}"
        yields[cat]={}
        for proc in dstore.keys():
            dta=dstore[proc][cat]['nominal']
#             print(dta)
            if len(dta) < 4:
        #         raise ValueError
                if proc in ["WToQQHHTo2B2G","ZToBBHHTo2B2G"]:
                    dta={ky:dta[ky]*3.315 for ky in dta}
                elif proc in ["ggHH"]:
                    if 'UL16Post' not in dta:
                        dta['UL16Post']=dta['UL16Pre']*0.8615
                elif proc in ['morphedData']:
                    pass
                else:
#                     print(proc,dta.keys())
                    continue
        #     print(proc,dta)
            yields[cat][proc]=sum([dta[ky] for ky in dta])
#             break
        for ky in merger_maps:
            yields[cat][ky]=0.0
            for pr in merger_maps[ky]:
                yields[cat][ky]+=yields[cat][pr]
        yieldStore[tag]=yields

procs_to_print=[
    'mc_ttJJ','mc_ttGJ','mc_ttGG','morphedData',
    'GGJ' , 
    'ggH' ,'ttH','vbfH','bbh','vH',  'H',
    'ggHH','ttHH', 'vHH','HH',
    'ExpectedAll','data', 'ggHHH', 
]

if not args.doCAT:
    yield_table=[["-","ggHHH", "NR Backgroud Exp. SR","H","HH","Data SB"]]
    for i in range(5):
        cat=f'CAT{i}'
        yield_table.append([cat])    
        
        for proc in ["ggHHH","Expected","H","HH"]:
            yld=yieldStore['SR'][cat][proc]
            if yld[0] > 1.0:
                yield_table[-1]+=[f"{yld[0]:>.2f} +/- {np.sqrt(yld[1]):>.2f}"]
            elif yld[0] > 1e-3:
                yield_table[-1]+=[f"{yld[0]:>.4f} +/- {np.sqrt(yld[1]):>.4f}"]
            else:
                yield_table[-1]+=[f"-"]
        yield_table[-1]+=[f"{yieldStore['SB'][cat]['data'][0]:.0f} +/- {float(np.sqrt(yieldStore['SB'][cat]['data'][1])):.2f}"]

    if args.doLatex:
        tab_str=tab.tabulate(yield_table,tablefmt='latex').replace("\n"," \hline\n")
        tab_str=tab_str.replace('+/-','$\pm$')
        print(tab_str)
    else:
        print(tab.tabulate(yield_table))


else:
    cat=args.cat
    yield_table=[["Process","Sideband Yield","Signal Region Yield"]]
    for proc in procs_to_print:
        yld_sr=yieldStore['SR'][cat][proc]
        yld_sb=yieldStore['SB'][cat][proc]
        c=np.log10(yld_sr)
        yield_table.append([f"{proc:>12}"])
        if yld_sb[0] > 1.0:
            yield_table[-1]+=[f"{yld_sb[0]:>.2f} +/- {np.sqrt(yld_sb[1]):>.2f}"]
        elif yld_sb[0] > 1e-3:
            yield_table[-1]+=[f"{yld_sb[0]:>.4f} +/- {np.sqrt(yld_sb[1]):>.4f}"]
        else:
            yield_table[-1]+=[f"-"]
            
        
        if yld_sr[0] > 1.0:
            yield_table[-1]+=[f"{yld_sr[0]:>.2f} +/- {np.sqrt(yld_sr[1]):>.2f}"]
        elif yld_sr[0] > 1e-3:
            yield_table[-1]+=[f"{yld_sr[0]:>.4f} +/- {np.sqrt(yld_sr[1]):>.4f}"]
        else:
            yield_table[-1]+=[f"-"]
            
        
    if args.doLatex:
        tab_str=tab.tabulate(yield_table,tablefmt='latex').replace("\n"," \hline\n")
        tab_str=tab_str.replace('+/-','$\pm$')
        print(tab_str)
    else:
        print(tab.tabulate(yield_table))



# tag_map={
#     "mc\_ttJJ" : r"\PQt\PAQt+ jets",
#     "mc\_ttGJ" : r"\PQt\PAQt+\PGg+ jets",
#     "mc\_ttGG" : r"\PQt\PAQt+\PGg\PGg",  
#     "morphedta" :"Morphed LSR",        
#     "GGJ" :      "\PGg\PGg + jets",    
#     "ggH" :      "ggH",        
#     "ttH" :      "ttH",        
#     "vbfH" :     "VBFH",       
#     "bbh" :      "bbH",        
#     "vH" :       "VH",                
#     "H" :        "{\bf Total H }",    
#     "ggHH" :     "ggHH",       
#     "ttHH" :     "ttHH",    
#     "vHH"  :     "VHH",               
#     "HH" :       "{\bf Total HH}",    
#     "Expected" : "Total MC Bkg.",    
#     "data" :     "Data",             
#     "ggHHH" :    "\textbf{ggHHH}"    
# }


# for i in range(len(yield_table)):
#     raw=yield_table[i]
#     if raw[0].strip() in ['H','HH','GGJ']:
#         yield_table[i][2]+=' \hline'
#     if raw[0].strip() in tag_map:
#         yield_table[i][0]=tag_map[raw[0].strip()]
        
