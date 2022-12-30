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
    - misc/makeEventIDUpdater_data.sh

- `python/ntupleMaker_mlScoreAdder.py`
    - Adds the ML model score from the `pickle`-ed score values
    - Helper class : `python/trippleHiggsMLInterface.py`
    - misc/makeScoreUpdater_data.sh
    
