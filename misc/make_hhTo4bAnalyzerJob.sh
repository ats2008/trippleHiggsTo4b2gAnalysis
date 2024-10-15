version=h3sin61_4p0

./misc/condorJobMakerGeneric.py \
       --exe            /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntupleMaker_h24bseletion.py \
       --fsrc           fileList/ntuples/ggHHH_UL2018_mlScoreUpdated.fls \
       --runScript      misc/runPython.tpl.sh \
       --cfg            misc/cfg/ntupleHHTo4b_Signal.tpl.cfg  \
       --dest           results/analysis/hhTo4bSelection/ggHHH_$version \
       --jn             100 \
       --fn             1 \
       --maxEvt         -1 \
       --tag            ggHHH_hhTo4b_$version \
       --maxMeterialize 250 

