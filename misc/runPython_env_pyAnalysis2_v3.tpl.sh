#!/bin/bash
source /grid_mnt/t3home/athachay/.bashrc
source /cvmfs/cms.cern.ch/cmsset_default.sh 
set -x
export HOME=@@HOME
export X509_USER_PROXY=@@proxy_path
cd @@DIRNAME
set +x
#eval `scramv1 runtime -sh`
conda activate /home/athachay/miniconda3/envs/pyAnalysis2
set -x

if [[ -z "${_CONDOR_SCRATCH_DIR}" ]]; then
    _CONDOR_SCRATCH_DIR=`mktemp -d`
    echo _CONDOR_SCRATCH_DIR was not set, setting it to $_CONDOR_SCRATCH_DIR
else
    echo CONDOR_SCRATCH_DIR is $_CONDOR_SCRATCH_DIR
fi
cd $_CONDOR_SCRATCH_DIR

cp @@DIRNAME/@@CFGFILENAME .
mv @@RUNSCRIPT @@RUNSCRIPT.busy 
python @@EXECUTABLE --cfg @@CFGFILENAME
if [ $? -eq 0 ]; then 
    cp *.root @@DESTINATION
    if [ $? -eq 0 ] ; then
        mv @@RUNSCRIPT.busy @@RUNSCRIPT.sucess 
        echo OK
    else
        mv @@RUNSCRIPT.busy @@RUNSCRIPT
        echo FAIL
    fi
else
    mv @@RUNSCRIPT.busy @@RUNSCRIPT 
    echo FAIL
fi
