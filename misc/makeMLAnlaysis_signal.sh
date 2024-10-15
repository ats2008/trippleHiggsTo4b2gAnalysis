
version=6p0a1

./misc/condorJobMakerGeneric.py \
    --exe                /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntupleMaker_ML_AnalysisSkim.py \
    --fsrc               fileList/ntuples/ggHHH_UL2018_mlScoreUpdated.fls \
    --runScript          misc/runPython.tpl.sh \
    --cfg                misc/cfg/mlanalsis_signal_2018.tpl.cfg \
    --dest               results/analysis/analysisSkims/mlAnalysis_Signal2018_$version/ \
    --jn                 10000 \
    --fn                 1 \
    --maxEvt             -1 \
    --tag                mlAnalysisSignal2018_$version \
    --maxMeterialize     250
