from CRABClient.UserUtilities import config 
config = config()

config.General.requestName = ''
config.General.transferLogs = True
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nanoanalyzercrab_cfg_22_DoubleEle.py'
config.JobType.maxMemoryMB = 4000
config.JobType.numCores = 1

config.Data.inputDataset = ''
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2

config.General.workArea = 'crab_elejet_lastversion/crab_elejet_lastversion'

#name='/EGamma/Run2022C-PromptReco-v1/MINIAOD/DoubleEleTest'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/EGamma/Run2022C-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022C-PromptReco-v1_DoubleEleTest'
#config.Data.outLFNDirBase = '/store/user/crovelli/EGamma2022Last/%s' % (config.General.workArea)

#name='/EGamma/Run2022D-PromptReco-v1/MINIAOD/DoubleEleTest'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/EGamma/Run2022D-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022D-PromptReco-v1_DoubleEleTest'
#config.Data.outLFNDirBase = '/store/user/crovelli/EGamma2022Last/%s' % (config.General.workArea)

#name='/EGamma/Run2022D-PromptReco-v2/MINIAOD/DoubleEleTest'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/EGamma/Run2022D-PromptReco-v2/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022D-PromptReco-v2_DoubleEleTest'
#config.Data.outLFNDirBase = '/store/user/crovelli/EGamma2022Last/%s' % (config.General.workArea)
##config.Data.outLFNDirBase = '/store/group/phys_bphys/crovelli/triggerRun3/EGamma2022/%s' % (config.General.workArea)   

#name='/EGamma/Run2022E-PromptReco-v1/MINIAOD/DoubleEleTest'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/EGamma/Run2022E-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022E-PromptReco-v1_DoubleEleTest'
#config.Data.outLFNDirBase = '/store/user/crovelli/EGamma2022Last/%s' % (config.General.workArea)

#name='/EGamma/Run2022F-PromptReco-v1/MINIAOD/DoubleEleTest'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/EGamma/Run2022F-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022F-PromptReco-v1_DoubleEleTest'
#config.Data.outLFNDirBase = '/store/user/crovelli/EGamma2022Last/%s' % (config.General.workArea)

name='/EGamma/Run2022G-PromptReco-v1/MINIAOD/DoubleEleTest'
config.General.requestName = name.replace('/','_')
config.Data.inputDataset = '/EGamma/Run2022G-PromptReco-v1/MINIAOD'
config.Data.outputDatasetTag = 'Run2022G-PromptReco-v1_DoubleEleTest'
config.Data.outLFNDirBase = '/store/user/crovelli/EGamma2022Last/%s' % (config.General.workArea)


config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_IT_Rome'
#config.Site.storageSite = 'T2_CH_CERN'

