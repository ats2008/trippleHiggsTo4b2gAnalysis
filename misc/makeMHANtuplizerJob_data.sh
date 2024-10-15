version=6p0

for era in A B C D;
    do
    ./misc/condorJobMakerGeneric.py \
           --exe                /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntuple_forMHA_Background.py \
           --fsrc               fileList/ntuples/data2018${era}_v2_upd.fls \
           --runScript          misc/runPython.tpl.sh \
           --cfg                misc/ntuplizerData.tpl.cfg \
           --dest               results/data/ntuples/UL18/ml_data2018${era}_$version/ \
           --jn                 10000 \
           --fn                 1 \
           --maxEvt             -1 \
           --tag                data2018${era}_$version \
           --maxMeterialize     250  ;
     done;
#
#./misc/condorJobMakerGeneric.py \
#       --exe                /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntuple_forMHA_Background.py \
#       --fsrc               fileList/ntuples/data2018B_v2.fls \
#       --runScript          misc/runPython.tpl.sh \
#       --cfg                misc/ntuplizerData.tpl.cfg \
#       --dest               results/data/ntuples/UL18/ml_data2018B_$version/ \
#       --jn                 10000 \
#       --fn                 1 \
#       --maxEvt             -1 \
#       --tag                data2018B_$version \
#       --maxMeterialize     250
#
#./misc/condorJobMakerGeneric.py \
#       --exe                /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntuple_forMHA_Background.py \
#       --fsrc               fileList/ntuples/data2018C_v2.fls \
#       --runScript          misc/runPython.tpl.sh \
#       --cfg                misc/ntuplizerData.tpl.cfg \
#       --dest               results/data/ntuples/UL18/ml_data2018C_$version/ \
#       --jn                 10000 \
#       --fn                 1 \
#       --maxEvt             -1 \
#       --tag                data2018C_$version \
#       --maxMeterialize     250
#
#./misc/condorJobMakerGeneric.py \
#       --exe                /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntuple_forMHA_Background.py \
#       --fsrc               fileList/ntuples/data2018D_v2.fls \
#       --runScript          misc/runPython.tpl.sh \
#       --cfg                misc/ntuplizerData.tpl.cfg \
#       --dest               results/data/ntuples/UL18/ml_data2018D_$version/ \
#       --jn                 10000 \
#       --fn                 1 \
#       --maxEvt             -1 \
#       --tag                data2018D_$version \
#       --maxMeterialize     250
#
