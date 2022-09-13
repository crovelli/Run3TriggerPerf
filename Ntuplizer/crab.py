from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
config = config()

config.General.requestName = ''
config.General.transferLogs = True
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'nanoanalyzercrab_cfg_22_SingleMu.py'
config.JobType.psetName = 'nanoanalyzercrab_cfg_22_DoubleEle.py'
config.JobType.maxMemoryMB = 4000
config.JobType.numCores = 1

config.Data.inputDataset = ''
#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
#config.Data.totalUnits = 400

#config.Data.outLFNDirBase = '/store/user/jodedra/Run3/' #% (getUsernameFromSiteDB())
#config.Data.publication = True
#config.Data.outputDatasetTag = 'winter21'
#config.Site.ignoreGlobalBlacklist = True
#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T2_CH_*']
#config.Site.blacklist = ['T2_US_Purdue']


#config.General.workArea = 'crab_singlemu/crab_singlemumisseddataset'
config.General.workArea = 'crab_doubleele/crab_doubleelemisseddataset'

#name = 'EphemeralZeroBias8'
#name = 'ZeroBias10'

#name='/ParkingSingleMuon2/Run2022D-PromptReco-v2/MINIAODBPeak'
#config.Data.outLFNDirBase = '/store/user/jodedra/Run3SingleMunotrigmatch/20220831/'
#config.Data.inputDataset = '/ParkingSingleMuon2/Run2022D-PromptReco-v2/MINIAOD'

name='/EGamma/Run2022C-PromptReco-v1/MINIAOD/DoubleEleTest'
config.General.requestName = name.replace('/','_')
config.Data.inputDataset = '/EGamma/Run2022C-PromptReco-v1/MINIAOD'
config.Data.outputDatasetTag = 'Run2022C-PromptReco-v1_DoubleEleTest'
config.Data.outLFNDirBase = '/store/user/cquarant/EGamma2022/%s' % (config.General.workArea)

# name='/EGamma/Run2022D-PromptReco-v1/MINIAOD/DoubleEleTest'
# config.General.requestName = name.replace('/','_')
# config.Data.inputDataset = '/EGamma/Run2022D-PromptReco-v1/MINIAOD'
# config.Data.outputDatasetTag = 'Run2022D-PromptReco-v1_DoubleEleTest'
# config.Data.outLFNDirBase = '/store/user/cquarant/EGamma2022/%s' % (config.General.workArea)

# name='/EGamma/Run2022D-PromptReco-v2/MINIAOD/DoubleEleTest'
# config.General.requestName = name.replace('/','_')
# config.Data.inputDataset = '/EGamma/Run2022D-PromptReco-v2/MINIAOD'
# config.Data.outputDatasetTag = 'Run2022D-PromptReco-v2_DoubleEleTest'
# config.Data.outLFNDirBase = '/store/user/cquarant/EGamma2022/%s' % (config.General.workArea)

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_IT_Rome'

#config.Data.inputDataset = '/' + name + '/ytakahas-winter21-30e4efebaefe52740b9ca928c2409cd7/USER'


#for a in [1,2,3,4,5,6,7]:
    
#
#'/EphemeralZeroBias1/Run2018D-PromptReco-v2/MINIAOD'
#'/EphemeralZeroBias2/Run2018D-PromptReco-v2/MINIAOD'
#'/EphemeralZeroBias3/Run2018D-PromptReco-v2/MINIAOD'
#'/EphemeralZeroBias4/Run2018D-PromptReco-v2/MINIAOD'
#'/EphemeralZeroBias5/Run2018D-PromptReco-v2/MINIAOD'
#'/EphemeralZeroBias6/Run2018D-PromptReco-v2/MINIAOD'
#'/EphemeralZeroBias7/Run2018D-PromptReco-v2/MINIAOD'
#'/EphemeralZeroBias8/Run2018D-PromptReco-v2/MINIAOD'
#
#'/EphemeralZeroBias1/Run2018D-v1/RAW'
#'/EphemeralZeroBias2/Run2018D-v1/RAW'
#'/EphemeralZeroBias3/Run2018D-v1/RAW'
#'/EphemeralZeroBias4/Run2018D-v1/RAW'
#'/EphemeralZeroBias5/Run2018D-v1/RAW'
#'/EphemeralZeroBias6/Run2018D-v1/RAW'
#'/EphemeralZeroBias7/Run2018D-v1/RAW'
#'/EphemeralZeroBias8/Run2018D-v1/RAW'




#config.General.requestName = 'EphemeralZeroBias3'
#config.Data.inputDataset = '/EphemeralZeroBias3/Run2018D-PromptReco-v2/MINIAOD'
#config.Data.secondaryInputDataset = '/EphemeralZeroBias3/Run2018D-v1/RAW'



