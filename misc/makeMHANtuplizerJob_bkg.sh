version=6p0

./misc/condorJobMakerGeneric.py \
     --exe               /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntuple_forMHA_Background.py \
     --fsrc              fileList/ntuples/DiPhotonJetsBox_MGG-80toInf_13TeV.fls \
     --runScript         misc/runPython.tpl.sh \
     --cfg               misc/ntuplizerBkg.tpl.cfg \
     --dest              results/mc/ntuples/UL18/ml_diPhoton_$version/ \
     --jn                10000 \
     --fn                1 \
     --maxEvt            -1 \
     --tag               diPhoton_$version \
     --maxMeterialize    250

./misc/condorJobMakerGeneric.py \
      --exe              /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntuple_forMHA_Background.py \
      --fsrc             fileList/ntuples/DiPhotonJetsBox1BJet_MGG-80toInf_13TeV.fls \
      --runScript        misc/runPython.tpl.sh \
      --cfg              misc/ntuplizerBkg.tpl.cfg \
      --dest             results/mc/ntuples/UL18/ml_diPhoton1BJets_$version/ \
      --jn               10000 \
      --fn               1 \
      --maxEvt           -1 \
      --tag              diPhoton1BJets_$version \
      --maxMeterialize   250

./misc/condorJobMakerGeneric.py \
      --exe              /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntuple_forMHA_Background.py \
      --fsrc             fileList/ntuples/DiPhotonJetsBox2BJets_MGG-80toInf_13TeV.fls \
      --runScript        misc/runPython.tpl.sh \
      --cfg              misc/ntuplizerBkg.tpl.cfg \
      --dest             results/mc/ntuples/UL18/ml_diPhoton2BJets_$version/ \
      --jn               10000 \
      --fn               1 \
      --maxEvt           -1 \
      --tag              diPhoton2BJets_$version \
      --maxMeterialize   250


