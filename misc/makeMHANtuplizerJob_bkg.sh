version=1p0

./misc/condorJobMakerGeneric.py \
       /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntuple_forMHA_Background.py \
       fileList/ntuples/DiPhotonJetsBox_MGG-80toInf_13TeV.fls \
       misc/runPython.tpl.sh \
       misc/ntuplizerBkg.tpl.cfg \
       results/mc/ntuples/UL18/ml_diPhoton_$version/ \
       10000 \
       1 \
       -1 \
       diPhoton_$version \
       250

./misc/condorJobMakerGeneric.py \
       /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntuple_forMHA_Background.py \
       fileList/ntuples/DiPhotonJetsBox1BJet_MGG-80toInf_13TeV.fls \
       misc/runPython.tpl.sh \
       misc/ntuplizerBkg.tpl.cfg \
       results/mc/ntuples/UL18/ml_diPhoton1BJets_$version/ \
       10000 \
       1 \
       -1 \
       diPhoton1BJets_$version \
       250

./misc/condorJobMakerGeneric.py \
       /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/flashgg/CMSSW_10_6_29/analysis/python/ntuple_forMHA_Background.py \
       fileList/ntuples/DiPhotonJetsBox2BJets_MGG-80toInf_13TeV.fls \
       misc/runPython.tpl.sh \
       misc/ntuplizerBkg.tpl.cfg \
       results/mc/ntuples/UL18/ml_diPhoton2BJets_$version/ \
       10000 \
       1 \
       -1 \
       diPhoton2BJets_$version \
       250


