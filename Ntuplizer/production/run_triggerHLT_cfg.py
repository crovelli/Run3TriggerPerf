from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms
##               ##
#  CUSTOM OPTIONS #
##               ##
options = VarParsing('python')

options.register('isMC', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('Era','2022postEE',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set data taking period : 2022preEE, 2022postEE, 2023preBPix, 2023postBPix"
)
options.register('whichTrigger','HLT',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set the kind of trigger to be studied : L1, HLT"
)
options.register('globalTag', 'NOTSET',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set global tag"
)
options.register('wantSummary', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('wantFullRECO', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('reportEvery', 100,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "report every N events"
)
options.register('skip', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "skip first N events"
)

# set number of events
#options.setDefault('maxEvents', -1) 
options.setDefault('maxEvents', 500)
# set physic process
if not options.whichTrigger :
    which_trigger = 'HLT'
else :
    which_trigger = options.whichTrigger

# set task tag w.r.t. data taking period
if not options.Era : 
    era = '2022postEE'
else :
    era = options.Era
tag = '_'.join([which_trigger, era]) 
options.setDefault('tag', tag)
options.parseArguments()

# set global tag:
global_tags_mc = {
    '2022preEE'     : '130X_mcRun3_2022_realistic_v5',
    '2022postEE'    : '130X_mcRun3_2022_realistic_postEE_v6',
    '2023preBPix'   : '130X_mcRun3_2023_realistic_v14',
    '2023postBPix'  : '130X_mcRun3_2023_realistic_postBPix_v2',
}
global_tags_data = {
    '2022preEE'     : '124X_dataRun3_PromptAnalysis_v1', #era CD use '124X_dataRun3_Prompt_v10' for era E
    '2022postEE'    : '130X_dataRun3_PromptAnalysis_v1', #era FG
    '2023preBPix'   : '130X_dataRun3_PromptAnalysis_v1', #era BC
    '2023postBPix'  : '130X_dataRun3_PromptAnalysis_v1', #era D
}
 
if options._beenSet['globalTag']:
    globaltag = options.globalTag
else:
    globaltag = global_tags_mc[era] if options.isMC else global_tags_data[era] 

extension = {False : 'data', True : 'mc'}
outputFileTrigger = cms.string('TriggerStudyHLT.root')
#outputFileTrigger = cms.string('_'.join(['TriggerStudy', extension[options.isMC], options.tag])+'.root')

if not options.inputFiles :
    if options.isMC :
        options.inputFiles = ['/store/mc/Run3Summer22MiniAODv3/WtoTauNu_Tauto3Mu_TuneCP5_13p6TeV_pythia8/MINIAODSIM/124X_mcRun3_2022_realistic_v12-v2/2820000/0b14e03f-168c-4e39-b441-d1b949ee4890.root'] 
    else :
        options.inputFiles =  ['/store/data/Run2022D/Muon/MINIAOD/PromptReco-v1/000/357/538/00000/5eb5d59a-e3f1-46ba-91ad-9e2340ba0bc6.root']

##                   ##
#  START PRODUCTION   #  
##                   ##
#from Configuration.StandardSequences.Eras import eras
process = cms.Process('NanoAnalyzerDoubleMuHLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents=cms.untracked.uint32(options.skip),
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
)

process.TFileService = cms.Service("TFileService",
                                    fileName = outputFileTrigger
                                   )

# Analyzer
process.nano_ = cms.EDAnalyzer('NanoAnalyzerDoubleMu_HLT',
                              muons = cms.InputTag("slimmedMuons"),
                              vertices = cms.InputTag("offlineSlimmedPrimaryVertices"), 
                              HLT = cms.InputTag("TriggerResults","","HLT"),
                              triggerobjects = cms.InputTag("slimmedPatTrigger"),
                              l1MU = cms.InputTag("gmtStage2Digis", "Muon"),  
)

process.p = cms.Path(process.nano_)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')


from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
