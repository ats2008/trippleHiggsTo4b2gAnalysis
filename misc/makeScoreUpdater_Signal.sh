version=6p0
./misc/condorJobMakerGeneric.py \
       --exe       /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntupleMaker_mlScoreAdder.py \
       --fsrc      fileList/ntuples/ggHHH_UL2018_mlScoreUpdated.fls \
       --runScript misc/runPython.tpl.sh \
       --cfg       misc/cfg/ntupleMLScoreUpdater_Signal.tpl.cfg \
       --dest      results/mc/ntuples/UL18/ntupleUpdatedWithMLScore_ggHHH_$version/ \
       --jn        1000 \
       --fn        1 \
       --maxEvt    -1 \
       --tag       ggHHH_MLScoreUpdatedNtuples_$version \
       --maxMeterialize    250 

