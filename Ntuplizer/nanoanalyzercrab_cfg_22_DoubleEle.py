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
#process.load("Configuration.Geometry.GeometryIdeal_cff") # for 2011
#process.load("Configuration.StandardSequences.Geometry_cff") # for 2010
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

################################### VARPASSING ###################################
# Use VarParsing to specify your input directory and output file
# Comment Varparsing if you want to submit job using CRAB
# because CRAB does not support input directory as we put input dataset directly
#options = VarParsing('analysis')


#options.register(
#    "inputDir",
#    "",
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.string,
#    "Input directory with inputs"
#    )
#options.register(
#    "outputName",
#    "",
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.string,
#    "Output name"
#    )
#options.parseArguments()

#if (options.inputDir == ""):
#    sys.exit("Directory to find input file where???")
#else:
    # 2010 VM
    #InputDir = "/home/cms-opendata/CMSSW_4_2_8/NanoAOD/NanoAnalyzer/" + options.inputDir
    # 2011 NAF/HTC
    # InputDir = "/nfs/dust/cms/user/zulaiha/PhD/CMSSW_5_3_32/src/NanoAOD/NanoAnalyzer/" + options.inputDir
    # 2015 NAF/HTC
    #InputDir = "/nfs/dust/cms/user/zulaiha/PhD/CMSSW_7_6_1/src/NanoAOD/NanoAnalyzer/" + options.inputDir

##################################################################################

###################################### GLOBAL TAG ################################
# Change the global tag accordingly
# ParkingBPH 2021 UL
#process.GlobalTag.globaltag = '106X_dataRun2_v35'
# 2021 data
#process.GlobalTag.globaltag = '120X_dataRun3_Prompt_v2'
#process.GlobalTag.globaltag = '121X_dataRun3_v13'
# 2022 data
#process.GlobalTag.globaltag = '123X_dataRun3_Prompt_v8'
#process.GlobalTag.globaltag = '123X_dataRun3_Express_v10'
process.GlobalTag.globaltag = '123X_dataRun3_Prompt_v12'

##################################################################################

# Intialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Nano')
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
                                        
# Set the maximum number of events to be processed here
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#####################################  JSON FILE #################################
# Change the directory and JSON file accordingly
# Only uncomment if you run in Data
# ParkingBPH 2021 UL
#goodJSON = '/nfs/dust/cms/user/yangq2/goodJson/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
# pilot 2021
#goodJSON = './BT21GOOD_withTKon.json'
# 2022A, data taken up to June 3
#goodJSON = './Collisions22AGOOD_withALLON.json'


##################################################################################

# Get the luminosity list from the JSON
# Only uncomment if you run in Data
#myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')

##################################################################################

# Load jet correction services for all jet algoritms
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
#

#################################### INPUT FILE ##################################
# To test locally or submit batch job using condor, use this:
#fileinPut = FileUtils.loadListFromFile (InputDir)

#process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring(*fileinPut)
#)

# To submit batch job using CRAB or test locally (2nd option), use this:

#files = FileUtils.loadListFromFile("/eos/user/j/jodedra/BPARKINGNANOSTUFF/CMSSW_12_4_0_pre3/src/PhysicsTools/BParkingNano/test/singlefileboffifical.txt")
#files.extend(FileUtils.loadListFromFile("/eos/user/j/jodedra/BPARKINGNANOSTUFF/CMSSW_12_4_0_pre3/src/PhysicsTools/BParkingNano/test/singlefileboffifical.txt"))


inputFiles= [
    '/store/data/Run2022D/EGamma/MINIAOD/PromptReco-v2/000/357/899/00000/990a1665-d685-4623-8ddb-803c01962243.root',
]

process.source = cms.Source("PoolSource",
# for crab
	  fileNames = cms.untracked.vstring (inputFiles),
                            #fileNames = cms.untracked.vstring(*files)
)
#print(options.inputFiles)

##################################################################################

# Process the lumi
# Only uncomment if you run in Data
#process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#process.source.lumisToProcess.extend(myLumis)

process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string ('JPsi_ElePlusJet_controlTrigger.root')
                                   )


# Process the analyzer
process.nano_ = cms.EDAnalyzer('NanoAnalyzer',
                              electrons = cms.InputTag("slimmedElectrons"),
                              vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                              packedpfcandidates = cms.InputTag('packedPFCandidates'),
                              HLT = cms.InputTag("TriggerResults","","HLT"),
                              triggerobjects = cms.InputTag("slimmedPatTrigger"),
                              l1EG = cms.InputTag("caloStage2Digis", "EGamma"), 
                              # Change this:
                              # If HTC/VM:
                              #outFile = cms.string(options.outputName),
                              # If interactive:
                              #outFile = cms.string('test.root'), 
                              # If CRAB:
                              # make sure the name is same as the crab config
#                              outFile = cms.string('22ZeroBias136.root'),

                              # Change this:
                              # If MC:
                              #isData = cms.bool(False)
                              # If Data:
#                              isData = cms.bool(True)
)
process.p = cms.Path(process.nano_)
