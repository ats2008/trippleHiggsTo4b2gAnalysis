version=1p1_bkg

#./misc/condorJobMakerGeneric.py \
#        --exe            /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntuple_forMHA_Signal.py \
#        --fsrc           fileList/ntuples/ggHHH_UL2018_updated.fls \
#        --runScript      misc/runPython.tpl.sh \
#        --cfg            misc/ntuplizerForML_Signal_allEvents.tpl.cfg \
#        --dest           results/mc/ntuples/UL18/ml_ggHHH_allEvts_$version/ \
#        --jn             10000 \
#        --fn             1 \
#        --maxEvt         -1 \
#        --tag            ggHHH_allEvents_$version \
#        --maxMeterialize 250
./misc/condorJobMakerGeneric.py \
        --exe            /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntuple_forMHA_Signal.py \
        --fsrc           fileList/ntuples/ggHHH_UL2018_official.fls \
        --runScript      misc/runPython.tpl.sh \
        --cfg            misc/ntuplizerForML_Signal_allEvents.tpl.cfg \
        --dest           results/mc/ntuples/UL18/ml_ggHHH_official_allEvts_$version/ \
        --jn             100 \
        --fn             1 \
        --maxEvt         -1 \
        --tag            ggHHH_allEvents_$version \
        --maxMeterialize 250


