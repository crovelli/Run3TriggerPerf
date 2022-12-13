from CRABClient.UserUtilities import config 
config = config()

config.General.requestName = ''
config.General.transferLogs = True
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nanoanalyzercrab_cfg_22_DEusingSE.py'
config.JobType.maxMemoryMB = 2500
config.JobType.numCores = 1

config.Data.inputDataset = ''
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2

config.General.workArea = 'crab_singleele8_lowmass5/crab_singleele8_lowmass'

name='/ParkingDoubleElectronLowMass5/Run2022F-PromptReco-v1/MINIAOD/SingleEleTest'
config.General.requestName = name.replace('/','_')
config.Data.inputDataset = '/ParkingDoubleElectronLowMass5/Run2022F-PromptReco-v1/MINIAOD'
config.Data.outputDatasetTag = 'Run2022F-PromptReco-v1_SingleEleTest_LowMass5'
config.Data.outLFNDirBase = '/store/user/crovelli/HLT_SingleEle8_Last/%s' % (config.General.workArea)


config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

