python -W ignore python/ntupleMaker.py analysis_cfgs/signal_hhhTo4b2g.cfg > logs/signal_hhhTo4b2g_ntup.log &
python -W ignore python/ntupleMaker.py analysis_cfgs/ntup_diphotonInclusive.cfg  > logs/diphotonInclusive_ntup.log &
python -W ignore python/ntupleMaker.py analysis_cfgs/ntup_diphotonJetBox1.cfg  > logs/diphotonJetBox1_ntup.log &
python -W ignore python/ntupleMaker.py analysis_cfgs/ntup_diphotonJetBox2.cfg > logs/diphotonJetBox2_ntup.log &
python -W ignore python/ntupleMaker.py analysis_cfgs/ntup_recoAnalyzer.cfg > logs/data_ntup.log &
python -W ignore python/ntupleMaker.py analysis_cfgs/ntup_TGJets.cfg > logs/TGJets_ntup.log &
python -W ignore python/ntupleMaker.py analysis_cfgs/ntup_TTGG.cfg > logs/TTGG_ntup.log &
python -W ignore python/ntupleMaker.py analysis_cfgs/ntup_TTGJets.cfg > logs/TTGJets_ntup.log &
python -W ignore python/ntupleMaker.py analysis_cfgs/ntup_photonInclusive20To40.cfg > logs/photonInclusive20To40_ntup.log &
python -W ignore python/ntupleMaker.py analysis_cfgs/ntup_photonInclusive40ToInf.cfg > logs/photonInclusive40ToInf_ntup.log &