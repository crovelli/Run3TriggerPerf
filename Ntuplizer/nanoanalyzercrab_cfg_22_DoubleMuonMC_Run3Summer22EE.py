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
# 2022 MC Run3Summer22
process.GlobalTag.globaltag = '124X_mcRun3_2022_realistic_postEE_v1'

##################################################################################

# Intialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
                                        
# Set the maximum number of events to be processed here
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#################################### INPUT FILE ##################################
inputFiles= [
    '/store/mc/Run3Summer22EEMiniAODv3/JPsito2Mu_JPsiFilter_2MuFilter_TuneCP5_13p6TeV_pythia8/MINIAODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/2550000/00165039-2147-4f59-ac73-212a30d8c167.root',
]

process.source = cms.Source("PoolSource",
# for crab
	  fileNames = cms.untracked.vstring (inputFiles),
)

##################################################################################

process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string ('JPsi_DoubleMuon_MC.root')
                                   )

# Process the analyzer
process.nano_ = cms.EDAnalyzer('NanoAnalyzerDoubleMuMC',
                              muons = cms.InputTag("slimmedMuons"),
                              vertices = cms.InputTag("offlineSlimmedPrimaryVertices"), 
                              HLT = cms.InputTag("TriggerResults","","HLT"),
                              triggerobjects = cms.InputTag("slimmedPatTrigger"),
                              l1MU = cms.InputTag("gmtStage2Digis", "Muon", "RECO"),     
                              packedgenparticles = cms.InputTag("packedGenParticles"), 
)

process.p = cms.Path(process.nano_)
