
version=0p0

./misc/condorJobMakerGeneric.py \
    --exe                /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntupleMaker_analysisSkim_legacy.py \
    --fsrc               fileList/ntuples/ggHHH_UL2018_mlScoreUpdated.fls \
    --runScript          misc/runPython.tpl.sh \
    --cfg                misc/cfg/analsis_ggHHH_2018.tpl.cfg \
    --dest               results/analysis/analysisSkims/analysis_ggHHH_$version/ \
    --jn                 1000 \
    --fn                 1 \
    --maxEvt             -1 \
    --tag                analysis_ggHHH_$version \
    --maxMeterialize     250

