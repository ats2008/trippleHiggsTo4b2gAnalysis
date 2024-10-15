./misc/condorJobMakerGeneric.py \
      --exe /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/fastSimStudy/CMSSW_10_2_22/genAnalysis/python/genAnalyzer.py \
      --fsrc fileList/ntuples/c3_1_c4_1_HHHto4b2gamma_fullSim.fls \
      --runScript misc/runPython.tpl.sh \
      --cfg misc/recoAnalyzer.tpl.cfg \
      --dest results/mc/genAnalysis/v1p4/c3_1_c4_1_HHHto4b2gamma_fullSim/ \
      --jn 10000 \
      --fn 1 \
      --maxevents -1 \
      --tag GenAly_c3_1_c4_1_FS \
      --maxMeterialize 250


