# trippleHiggsTo4b2gAnalysis
Analysis scrpts for trippleH--> 4b2g analysis

- `python/ntupleMaker_analysisSkim_legacy.py` 
    - Does the hh->4b Selection using DHH criteria
    - Uses the Non Resonant MVA
        - MVA cut can be varied 
    - Has the option to use Peaking MVA values too
        - Might have to tune it properly
        - Default MVA file might be not the correct one
    -`misc/makeLegacyAnlaysis_data.sh` and `misc/makeLegacyAnlaysis_Signal.sh`

- `python/ntupleMaker_eventIDUpdater.py`
    - Modifies `event` number of the events .. old `event` moved to a new `eventIdx` branch
    - `misc/makeEventIDUpdater_data.sh`

- `python/ntupleMaker_mlScoreAdder.py`
    - Adds the ML model score from the `pickle`-ed score values
    - Helper class : `python/trippleHiggsMLInterface.py`
    - _depricated_ : misc/makeScoreUpdater_data.sh
    - `misc/addMLScoreToDset.py`

- `python/ntuple_forBDT.py`
    - Makes the Ntuples for BDT Training And testing
    - `misc/makeBDTNtuplizerJob.py`

- `python/variableAnalyzer.py`
    - Make the histograms of all the variables required
    - Compiled scalbale workflow from `trippleHiggsUtils.py` and `trippleHiggsSelector.py`
    - `misc/makeVariableDistributions.py`

- `python/genAnalyzer.py`
    - Makes the kinematic as well as other plots from gen Only information
    - Works on the files made using `genNtuplizer`
    - Template cfg : `misc/cfg/genAnalyzer.cfg`

- `python/deriveScaleFactors.py`
    - Derives the Scale factors , makes basic check-me plots 
    - Stores the results in the destination
    - Based on local wflow
    - Uses mplhep, rootdataframes
