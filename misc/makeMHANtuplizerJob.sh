version=1p1

./misc/condorJobMakerGeneric.py \
       /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntuple_forMHA_Signal.py \
       fileList/ntuples/ggHHH_UL2018.fls \
       misc/runPython.tpl.sh \
       misc/ntuplizerSignal.tpl.cfg \
       results/mc/ntuples/UL18/ml_ggHHH_$version/ \
       10000 \
       1 \
       -1 \
       ggHHH_$version \
       250

