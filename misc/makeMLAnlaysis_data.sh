version=6p0a1

for era in A B C D;
    do 
        echo Making for Era $i ;
        ./misc/condorJobMakerGeneric.py \
            --exe                /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntupleMaker_ML_AnalysisSkim.py \
            --fsrc               fileList/ntuples/data2018${era}_v2_mlupd.fls \
            --runScript          misc/runPython.tpl.sh \
            --cfg                misc/cfg/mlanalsis_data_2018.tpl.cfg \
            --dest               results/analysis/analysisSkims/mlAnalysis_data2018${era}_$version/ \
            --jn                 10000 \
            --fn                 1 \
            --maxEvt             -1 \
            --tag                mlAnalysisData2018${era}_$version \
            --maxMeterialize     250
    done

ls Condor/*mlAnalysisData2018*_${version}*/*.sub
