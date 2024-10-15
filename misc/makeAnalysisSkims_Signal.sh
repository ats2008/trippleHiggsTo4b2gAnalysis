version=3p1
./misc/condorJobMakerGeneric.py \
       --exe       /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntupleMaker_analysisSkimer_v1.py \
       --fsrc      fileList/ntuples/ggHHH_UL2018_mlScoreUpdated.fls \
       --runScript misc/runPython.tpl.sh \
       --cfg       misc/cfg/analyzer_Signal.tpl.cfg \
       --dest      results/analysis/analysisSkims/
       --jn        10000 \
       --fn        1 \
       --maxEvt    -1 \
       --tag       ggHHH_AnalysisV1_$version \
       --maxMeterialize    250 

