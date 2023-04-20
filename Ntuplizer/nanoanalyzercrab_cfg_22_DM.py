##################################################################################
# Nanoanalyzer configuration file for all years                                  #
# Use for HT Condor and VM                                                       #
# Uncomment & comment relevant lines before you run it                           #
##################################################################################

import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Types as CfgTypes
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.Utilities.FileUtils as FileUtils
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("Nano")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

###################################### GLOBAL TAG ################################
# 2022 data
process.GlobalTag.globaltag = '123X_dataRun3_Prompt_v12'

##################################################################################

# Intialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
                                        
# Set the maximum number of events to be processed here
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#################################### INPUT FILE ##################################
# To test locally or submit batch job using condor, use this:
#fileinPut = FileUtils.loadListFromFile (InputDir)

#/ParkingDoubleElectronLowMass0/Run2022F-PromptReco-v1/MINIAOD  --> LowMass5

inputFiles= [
    '/store/data/Run2022D/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v2/000/357/734/00000/0f35818c-326c-4d85-a190-3bdec4664b2c.root',
    '/store/data/Run2022D/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v2/000/357/734/00000/1fc635c0-d59f-4a5a-90be-de3cfe51ed83.root',
    '/store/data/Run2022D/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v2/000/357/734/00000/21b36205-62c0-4104-a4b8-3f663827caaa.root',
    '/store/data/Run2022D/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v2/000/357/734/00000/2d045926-8dd5-4312-b330-75ea1af8680a.root'
]

process.source = cms.Source("PoolSource",
# for crab
	  fileNames = cms.untracked.vstring (inputFiles),
)

##################################################################################

process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string ('JPsi_DoubleMuon.root')
                                   )

# Process the analyzer
process.nano_ = cms.EDAnalyzer('NanoAnalyzerDoubleMu',
                              muons = cms.InputTag("slimmedMuons"),
                              vertices = cms.InputTag("offlineSlimmedPrimaryVertices"), 
                              HLT = cms.InputTag("TriggerResults","","HLT"),
                              triggerobjects = cms.InputTag("slimmedPatTrigger"),
                              l1MU = cms.InputTag("gmtStage2Digis", "Muon"),  
)

process.p = cms.Path(process.nano_)
