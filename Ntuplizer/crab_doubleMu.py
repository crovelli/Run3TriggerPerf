from CRABClient.UserUtilities import config 
config = config()

config.General.requestName = ''
config.General.transferLogs = True
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'nanoanalyzercrab_cfg_22_DM.py'
config.JobType.numCores = 1

config.Data.inputDataset = ''
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Muon.json'
config.General.workArea = 'crab_doubleMu_lowmass0_ter/crab_doubleMu'

#name='/ParkingDoubleMuonLowMass0/Run2022D-PromptReco-v2/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass0/Run2022D-PromptReco-v2/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022D-PromptReco-v2_DoubleMuTnP_LowMass0'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass1/Run2022D-PromptReco-v2/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass1/Run2022D-PromptReco-v2/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022D-PromptReco-v2_DoubleMuTnP_LowMass1'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass2/Run2022D-PromptReco-v2/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass2/Run2022D-PromptReco-v2/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022D-PromptReco-v2_DoubleMuTnP_LowMass2'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass3/Run2022D-PromptReco-v2/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass3/Run2022D-PromptReco-v2/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022D-PromptReco-v2_DoubleMuTnP_LowMass3'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass4/Run2022D-PromptReco-v2/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass4/Run2022D-PromptReco-v2/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022D-PromptReco-v2_DoubleMuTnP_LowMass4'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass5/Run2022D-PromptReco-v2/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass5/Run2022D-PromptReco-v2/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022D-PromptReco-v2_DoubleMuTnP_LowMass5'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass6/Run2022D-PromptReco-v2/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass6/Run2022D-PromptReco-v2/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022D-PromptReco-v2_DoubleMuTnP_LowMass6'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass7/Run2022D-PromptReco-v2/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass7/Run2022D-PromptReco-v2/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022D-PromptReco-v2_DoubleMuTnP_LowMass7'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

# -------------------------------------------------------------------------------

#name='/ParkingDoubleMuonLowMass0/Run2022E-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass0/Run2022E-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022E-PromptReco-v1_DoubleMuTnP_LowMass0'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass1/Run2022E-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass1/Run2022E-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022E-PromptReco-v1_DoubleMuTnP_LowMass1'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass2/Run2022E-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass2/Run2022E-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022E-PromptReco-v1_DoubleMuTnP_LowMass2'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass3/Run2022E-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass3/Run2022E-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022E-PromptReco-v1_DoubleMuTnP_LowMass3'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass4/Run2022E-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass4/Run2022E-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022E-PromptReco-v1_DoubleMuTnP_LowMass4'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass5/Run2022E-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass5/Run2022E-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022E-PromptReco-v1_DoubleMuTnP_LowMass5'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass6/Run2022E-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass6/Run2022E-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022E-PromptReco-v1_DoubleMuTnP_LowMass6'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass7/Run2022E-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass7/Run2022E-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022E-PromptReco-v1_DoubleMuTnP_LowMass7'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

# -------------------------------------------------------------------------

#name='/ParkingDoubleMuonLowMass0/Run2022F-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass0/Run2022F-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022F-PromptReco-v1_DoubleMuTnP_LowMass0'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu20/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass1/Run2022F-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass1/Run2022F-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022F-PromptReco-v1_DoubleMuTnP_LowMass1'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu20/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass2/Run2022F-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass2/Run2022F-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022F-PromptReco-v1_DoubleMuTnP_LowMass2'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu20/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass3/Run2022F-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass3/Run2022F-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022F-PromptReco-v1_DoubleMuTnP_LowMass3'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu20/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass4/Run2022F-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass4/Run2022F-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022F-PromptReco-v1_DoubleMuTnP_LowMass4'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu20/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass5/Run2022F-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass5/Run2022F-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022F-PromptReco-v1_DoubleMuTnP_LowMass5'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu20/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass6/Run2022F-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass6/Run2022F-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022F-PromptReco-v1_DoubleMuTnP_LowMass6'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu20/%s' % (config.General.workArea)

name='/ParkingDoubleMuonLowMass7/Run2022F-PromptReco-v1/MINIAOD/DoubleMuTnP'
config.General.requestName = name.replace('/','_')
config.Data.inputDataset = '/ParkingDoubleMuonLowMass7/Run2022F-PromptReco-v1/MINIAOD'
config.Data.outputDatasetTag = 'Run2022F-PromptReco-v1_DoubleMuTnP_LowMass7'
config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu20/%s' % (config.General.workArea)

# -------------------------------------------------------------------------

#name='/ParkingDoubleMuonLowMass0/Run2022G-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass0/Run2022G-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022G-PromptReco-v1_DoubleMuTnP_LowMass0'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass1/Run2022G-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass1/Run2022G-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022G-PromptReco-v1_DoubleMuTnP_LowMass1'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass2/Run2022G-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass2/Run2022G-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022G-PromptReco-v1_DoubleMuTnP_LowMass2'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass3/Run2022G-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass3/Run2022G-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022G-PromptReco-v1_DoubleMuTnP_LowMass3'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass4/Run2022G-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass4/Run2022G-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022G-PromptReco-v1_DoubleMuTnP_LowMass4'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass5/Run2022G-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass5/Run2022G-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022G-PromptReco-v1_DoubleMuTnP_LowMass5'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass6/Run2022G-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass6/Run2022G-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022G-PromptReco-v1_DoubleMuTnP_LowMass6'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

#name='/ParkingDoubleMuonLowMass7/Run2022G-PromptReco-v1/MINIAOD/DoubleMuTnP'
#config.General.requestName = name.replace('/','_')
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass7/Run2022G-PromptReco-v1/MINIAOD'
#config.Data.outputDatasetTag = 'Run2022G-PromptReco-v1_DoubleMuTnP_LowMass7'
#config.Data.outLFNDirBase = '/store/user/crovelli/HLT_DoubleMu/%s' % (config.General.workArea)

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

