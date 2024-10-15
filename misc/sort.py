import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i',"--input", help="Input file ",default=None)
parser.add_argument('-n',"--numberToPrint", help="Number of inputs to print",default=20,type=int)
args = parser.parse_args()

with open(args.input) as f:
    txt=f.readlines()
fldr=[]
size=[]

for l in txt:
    kv=l[:-1].strip().split()
    fldr.append(kv[1])
    size.append(eval(kv[0].replace("M","*1000").replace("G","*1e6").replace("T","1e9").replace("K","*1e0"))/1e3)

srt_idxs=np.argsort(np.array(size)*-1)
for j in range(20):
    i=srt_idxs[j]
    print(fldr[i], "  :  " ,size[i]/1e3," GB")
