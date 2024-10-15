version=3p01

./misc/condorJobMakerGeneric.py \
       --exe                /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntupleMaker_mlScoreAdder.py \
       --fsrc               fileList/ntuples/DiPhotonJetsBox_MGG-80toInf_13TeV_mlScoreUpdated.fls \
       --runScript          misc/runPython.tpl.sh \
       --cfg                misc/cfg/ntupleMLScoreUpdater_Bkg_qcdInc.tpl.cfg \
       --dest               results/mc/ntuples/UL18/ntupleUpdatedWithMLScore_diPhoton_$version/ \
       --jn                 1000 \
       --fn                 1 \
       --maxEvt            -1 \
       --tag                MLScoreUpdatedNtuples_diPhoton_$version \
       --maxMeterialize     250

./misc/condorJobMakerGeneric.py \
       --exe              /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntupleMaker_mlScoreAdder.py  \
       --fsrc             fileList/ntuples/DiPhotonJetsBox1BJet_MGG-80toInf_13TeV_mlScoreUpdated.fls \
       --runScript        misc/runPython.tpl.sh \
       --cfg              misc/cfg/ntupleMLScoreUpdater_Bkg_qcd1b.tpl.cfg \
       --dest             results/mc/ntuples/UL18/ntupleUpdatedWithMLScore_diPhoton1BJets_$version/ \
       --jn               1000 \
       --fn               1 \
       --maxEvt           -1 \
       --tag              MLScoreUpdatedNtuples_diPhoton1BJets_$version \
       --maxMeterialize   250

./misc/condorJobMakerGeneric.py \
       --exe              /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntupleMaker_mlScoreAdder.py \
       --fsrc             fileList/ntuples/DiPhotonJetsBox2BJets_MGG-80toInf_13TeV_mlScoreUpdated.fls \
       --runScript        misc/runPython.tpl.sh \
       --cfg              misc/cfg/ntupleMLScoreUpdater_Bkg_qcd2b.tpl.cfg \
       --dest             results/mc/ntuples/UL18/ntupleUpdatedWithMLScore_diPhoton2BJets_$version/ \
       --jn               1000 \
       --fn               1 \
       --maxEvt           -1 \
       --tag              MLScoreUpdatedNtuples_diPhoton2BJets_$version \
       --maxMeterialize   250

