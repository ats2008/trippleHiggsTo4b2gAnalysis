
version=0p1

for era in A B C D;
    do 
        echo Making for Era $i ;
        ./misc/condorJobMakerGeneric.py \
            --exe                /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntupleMaker_analysisSkim_legacy.py \
            --fsrc               fileList/ntuples/data2018${era}_v2.fls \
            --runScript          misc/runPython.tpl.sh \
            --cfg                misc/cfg/analsis_data_2018.tpl.cfg \
            --dest               results/analysis/analysisSkims/analysis_data2018${era}_$version/ \
            --jn                 1000 \
            --fn                 1 \
            --maxEvt             -1 \
            --tag                analysisData2018${era}_$version \
            --maxMeterialize     250
    done

ls Condor/*analysisData2018*_${version}*/*.sub
