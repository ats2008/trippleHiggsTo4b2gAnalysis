
version=5p2

for era in A B C D ;
    do
       echo Doing $era
      ./misc/condorJobMakerGeneric.py \
           --exe                /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntupleMaker_mlScoreAdder.py \
           --fsrc               2018${era}.fls \
           --runScript          misc/runPython.tpl.sh \
           --cfg                misc/cfg/ntupleMLScoreUpdater_data_2018${era}.tpl.cfg \
           --dest               results/data/ntuples/UL18/ntupleUpdatedWithMLScore_data2018${era}_$version/ \
           --jn                 10000 \
           --fn                 1 \
           --maxEvt             -1 \
           --tag                mlUpd1_data2018${era}_$version \
           --maxMeterialize     250
    done

ls -ltrh Condor/*mlUpd1_data2018*_$version/*.sub 

#           --fsrc               fileList/ntuples/data2018${era}_v2_upd.fls \
