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
