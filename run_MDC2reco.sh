#!/bin/bash

export HOME=/sphenix/u/${LOGNAME}
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.322


nEvents=0

inputFiles="{"
for fileList in $@
do
  inputFiles+="\"${fileList}\","
done
inputFiles=${inputFiles::-1}
inputFiles+="}"
echo running: run_HFreco.sh $*
root.exe -q -b Fun4All_MDC2reco.C\(${inputFiles},$nEvents\)
echo Script done
