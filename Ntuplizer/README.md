# Trigger studies for Run3

cfgs are for 2022 run

# Setup
cmsrel CMSSW_12_4_3

cd $CMSSW_BASE/src

cmsenv

git clone git@github.com:crovelli/Run3TriggerPerf.git

git checkout -b run2022

cd Run3TriggerPerf/Ntuplizer/

scram b -j 8


# Code
* Double Electron using doubleEG L1: 

plugins/NanoAnalyzer.cc

nanoanalyzercrab_cfg_22_DoubleEle.py

crab_doubleEle.py


* Double Electron using singleEG L1: 

plugins/NanoAnalyzerSingleEle.cc

nanoanalyzercrab_cfg_22_DEusingSE.py

crab_singleEle.py


* Double Muon

plugins/NanoAnalyzerDoubleMu.cc 

nanoanalyzercrab_cfg_22_DM.py

crab_doubleMu.py


* Double Muon, MC

plugins/NanoAnalyzerDoubleMuMC.cc

nanoanalyzercrab_cfg_22_DoubleMuonMC.py

crab_doubleMuMC.py


# To submit crab jobs 
crab submit crab.py
