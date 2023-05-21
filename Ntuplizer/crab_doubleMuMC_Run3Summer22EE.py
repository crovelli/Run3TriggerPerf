from CRABClient.UserUtilities import config 
config = config()

config.General.requestName = ''
config.General.transferLogs = True
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nanoanalyzercrab_cfg_22_DoubleMuonMC_Run3Summer22EE.py'
config.JobType.numCores = 1

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
config.General.workArea = 'crab_doubleMuMC_Run3SummerEE'

name='/JPsito2Mu_MCTnP'
config.General.requestName = name.replace('/','_')
config.Data.inputDataset = '/JPsito2Mu_JPsiFilter_2MuFilter_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2/MINIAODSIM'
config.Data.outputDatasetTag = 'JPsito2Mu_MCTnP'
config.Data.outLFNDirBase = '/store/user/crovelli/DoubleMuMC/%s' % (config.General.workArea)

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

