# Run 3 trigger performance
# ------------------------------

cmsrel CMSSW_13_0_13

cd CMSSW_13_0_13/src

cmsenv

git cms-init

git clone -b run2022and2023_splitL1HLT_13_0_13 --single-branch git@github.com:crovelli/Run3TriggerPerf.git 

cd Run3TriggerPerf/Ntuplizer

