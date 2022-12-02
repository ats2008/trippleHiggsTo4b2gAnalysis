version=0p1

for era in A B C D;
    do
    ./misc/condorJobMakerGeneric.py \
           --exe                /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntupleMaker_eventIDUpdater.py \
           --fsrc               fileList/ntuples/data2018${era}_v2.fls \
           --runScript          misc/runPython.tpl.sh \
           --cfg                misc/cfg/ntuplizerDataEvtUpdater.tpl.cfg \
           --dest               results/data/ntuples/UL18/data2018${era}_upd_$version/ \
           --jn                 1000 \
           --fn                 1 \
           --offsetStep         20000 \
           --maxEvt             -1 \
           --tag                data2018${era}_$version \
           --maxMeterialize     250  ;
     done;
ll Condor/Jo*data2018${era}_${version}*/*.sub
