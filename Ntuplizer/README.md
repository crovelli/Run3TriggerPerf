# Setup
```bash
cmsrel CMSSW_12_4_3
cd $CMSSW_BASE/src
cmsenv
git clone git@github.com:DiElectronX/Run3TriggerPerf.git
cd Run3TriggerPerf/Ntuplizer/
scram b -j 8

```
# Local Run
for Double Electron
```bash
cmsRun nanoanalyzercrab_cfg_22_DoubleEle.py
```
The file to edit for Double Electron is plugins/NanoAnalyzer.cc



for Single Muon
```bash
cmsRun nanoanalyzercrab_cfg_22_Single.py
```
The file to edit for Single Mu is plugins/NanoAnalyzerSingleMu.cc

# CRAB job
To run crab job just submit crab.py file using
```bash
crab submit crab.py
```
As the output is a flat ntuple it must be stored using
```python
config.Site.storageSite = 'T3_CH_CERNBOX'
config.Data.outLFNDirBase = '/store/user/jodedra/Run3SingleMunotrigmatch/20220831/'
```
The outLFNDirBase will be the location where it is stored.
Crab jobs must be submitted for each PD


# MC SAMPLES LOCATIONS ZERO PU
BTOjpsiKEE
/BuTOjpsiKEE20220831fiftyMbettersplitting/jodedra-SUMMER22_MINIAOD-d5db235e2a58bcae594a314d29cbde75/USER
https://cmsweb.cern.ch/das/request?input=%2FBuTOjpsiKEE20220831fiftyMbettersplitting%2Fjodedra-SUMMER22_MINIAOD-d5db235e2a58bcae594a314d29cbde75%2FUSER&instance=prod%2Fphys03 
BTOKEE
/BuTOKEE20220826bettersplitting/jodedra-SUMMER22_MINIAOD-d5db235e2a58bcae594a314d29cbde75/USER
https://cmsweb.cern.ch/das/request?input=%2FBuTOKEE20220826bettersplitting%2Fjodedra-SUMMER22_MINIAOD-d5db235e2a58bcae594a314d29cbde75%2FUSER&instance=prod%2Fphys03
BTOPSI2SKEE
/BuTOpsi2sKEE20220831fiftyMbettersplitting/jodedra-SUMMER22_MINIAOD-d5db235e2a58bcae594a314d29cbde75/USER
https://cmsweb.cern.ch/das/request?input=%2FBuTOpsi2sKEE20220831fiftyMbettersplitting%2Fjodedra-SUMMER22_MINIAOD-d5db235e2a58bcae594a314d29cbde75%2FUSER&instance=prod%2Fphys03