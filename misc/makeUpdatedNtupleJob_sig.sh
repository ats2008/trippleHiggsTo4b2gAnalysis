version=2p3
job=updater
./misc/condorJobMakerGeneric.py \
       /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntupleUpdater.py \
       fileList/ntuples/ggHHH_UL2018.fls \
       misc/runPython.tpl.sh \
       misc/ntuplizerSignal.tpl.cfg \
       results/mc/updated_fggNtuples/UL18/ggHHH_$version/ \
       10000 \
       1 \
       -1 \
       ${job}_ggHHH_$version \
       250 \
       10001

