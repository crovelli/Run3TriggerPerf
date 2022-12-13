// **********************************************************************//
// Prototype implementation of NanoAnalyzer                              //
// for the creation of Run 1 nanoAOD-like ntuple from AOD or RECO input  //
// and corresp. Run 2 reference/validation ntuples from AOD or miniAOD   //
// **********************************************************************//

// ****************************************************
// for implementation history, see testnanoreadme.txt *
// ****************************************************
// Direct contributors:  A. Anuar, A. Bermudez, A. Geiser (coordinator), 
// N.Z. Jomhari, S. Wunsch, Q. Wang, H. Yang, Y. Yang, 2018-2021. 

// ******************************************
// automatically set appropriate CMSSW flag *
// ******************************************
// recognize automatically (no user action needed): 
// CMSSW is tied to particular compiler versions
// 42X taken from 4_2_8, 53X from 5_3_32, 7XX from 7_6_4
#define GCC_VERSION ( 10000 * __GNUC__ + 100 * __GNUC_MINOR__ + __GNUC_PATCHLEVEL__ )
#if GCC_VERSION < 40305
// Early Run 1 legacy CMSSW42X (e.g. 2010, 4_2_8)
#define CMSSW42X
#elif GCC_VERSION < 40703
// Main application is Run 1 legacy CMSSW53X (e.g. 2011/2, 5_3_32)
#define CMSSW53X
#elif GCC_VERSION > 40902
// instead in case CMSSW version is for Run 2, CMSSW7 or higher
// (e.g. 2015-18), for validation purposes
#define CMSSW7plus
#if GCC_VERSION < 50000
// for CMSSW 7_6_X  (GCC_VERSION might need to be changed/sharpened)
// (not clear whether this logic will still work for CMSSW 8 and 9)
#define CMSSW7XX
#endif
#endif
// for ultra-legacy AOD (10_6_4 or higher)
#if GCC_VERSION > 70400
// this one comes on top of the previous
#define CMSSW106plus
#endif
#if GCC_VERSION > 80300
// for Run 3 MC studies
// this one comes on top of the previous
// GCC_VERSION preliminary, might need to be changed/sharpened
#define CMSSW11plus
#endif
#if GCC_VERSION > 90299
// for 2021 pilot data 
// this one comes on top of the previous
#define CMSSW12plus
#endif

// ********************************************************************
// the following are flags which can be (de)activated manually
// defaults are set such that for normal use no manual action is needed
// ********************************************************************

// activate this to achieve maximum compatibility with official Run 2 nanoAOD
// (default is off -> better performance, e.g. use of beam spot constraint)
// (mainly for Run 2 validation - not recommended for Run 1)
// only relevant for Run 2 AOD, turn on for better consistency with miniAOD 
//                        (slightly worse performance)
// nanoext flag can also be steered/overruled via configuration file
//#define Compatibility

// activation flag to check and store compatibility with Golden JSON 
// default is on for 2010 data (not MC, automatic), off otherwise
#ifdef CMSSW42X
#define JSONcheck
#endif

#ifdef CMSSW42X
// activate this to check 2010 Golden JSON setting, i.e. abort upon nonJSON event
// activate this only when Golden JSON is activated in configuration
// (protection against accidental use of wrong or no JSON in configuration,
//  for check and validation, not strictly needed)
#define JSONcheckabort
#endif

// activate this only for data sets for which "plus" part of trigger treatment 
// is already implemented; checks and aborts in case of inconsistency;
// should be activated by default if trigger is implemented for dataset
// (protection against inconsistencies in NanoTrigger implementation)
//#define trigcheckabort

#ifdef CMSSW7plus
// activate this when you read from miniAOD for validation (Run 2/3 only!)
//#define miniAOD
#endif


// turn this on to deactivate code related to jet corrections
//  (e.g. in case of problems with the configuration, and only if you do 
//   *not* use jets in your analysis)
//  default should be flag off
#define noJetCor 

// turn this on to activate code related to charm final states
#define charm

#ifdef charm
// turn this on to activate D meson cuts optimized for large rapidities
// default (on) is for new looser cuts
#define beauty
#endif

// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <string>

#include "../interface/helperfunc.h"
#include "../interface/MyStruct.h"

#ifdef CMSSW42X
#include <boost/unordered_map.hpp>
using boost::unordered_map;
#else
#include <unordered_map>
using std::unordered_map;
#endif

//*****************************
// general user include files *
//*****************************
#include "FWCore/Framework/interface/Frameworkfwd.h"
#ifndef CMSSW12plus
// Run 1 and 2
#include "FWCore/Framework/interface/EDAnalyzer.h"
#else
// Run 3
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#endif
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// for fillDescriptions
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

// ------ EXTRA HEADER FILES--------------------//
#include "FWCore/Framework/interface/EventSetup.h"
#ifndef CMSSW12plus
#include "FWCore/Framework/interface/ESHandle.h"
#endif
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Ref.h"
// the following does not seem to be needed (included through dEdx), but also not to hurt
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
// to get release version
#include "FWCore/Version/interface/GetReleaseVersion.h"
// header for conversion tools
#ifndef CMSSW11plus
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#else
#include "CommonTools/Egamma/interface/ConversionTools.h"
#endif
// effective area for rho
//https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/RecoEgamma/EgammaTools/interface/EffectiveAreas.h
//#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

// for math and file handling
#include "TMath.h"
#include "TFile.h"
#include "Math/VectorUtil.h"

// for ntuple output
#include "TLorentzVector.h"
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#ifdef charm
// for control histograms
#include "TH1.h"
#include "TH2.h"
#endif

//**************************
// for trigger information *
//**************************
#include "FWCore/Common/interface/TriggerNames.h"
// not needed?
//#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/Common/interface/TriggerResults.h"
// for automatic trigger recognition
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
// for trigger objects
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include <cassert> 
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#ifndef CMSSW11plus
// no longer exists in CMSSW_11_2_X, not needed??
#include "HLTrigger/HLTcore/interface/TriggerSummaryAnalyzerAOD.h"
#endif

//***************************
// for tracking information *
//***************************
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
// for dEdx
#include "DataFormats/TrackReco/interface/DeDxData.h" 	

//#ifdef miniAOD
// tracks from PATCandidates
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
// there are three track collections in miniAOD (in addition to the muon 
// and electron collections), embedded into particleflow objects:
// packedPFCandidates allows to rebuild tracks (for pt>0.5 GeV) using
//                    pseudotrack() (object) or besttrack() (pointer)
// PackedPFCandidatesDiscarded presumably contains only discarded duplicate 
//                    muon candidates -> do not use
// lostTracks (high purity only) contains some of the non-vertex tracks 
//#endif

/// for track parametrization 
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/PerigeeTrajectoryParameters.h"

//*************************
// for vertex information *
//*************************
// reconstructed primary is typically within 0.02 cm of true primary in z
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/VertexReco/interface/Vertex.h>

// for vertices refit:
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFinder/interface/AdaptiveVertexReconstructor.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/TrimmedKalmanVertexFinder/interface/KalmanTrimmedVertexFinder.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerAlgorithm.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// for beamspot information (beam pipe radius: 5.8 -> 4.3 cm, 
//                           beam spot at x~0.2, y~0.4, z~0.3 cm,
//                           width x~0.002?, y~0.002?, z~5.5 cm,
// beam-spot constrained vertices: x~0.001 , y~0.001 , z~0.004 cm) 
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//****************************
// for error matrix handling *
//****************************
#include "TMatrixT.h"
#include "Math/SMatrix.h"
#include "Math/StaticCheck.h"
//#include <MatrixRepresentationsStatic.h>
#include "Math/Expression.h"

// set namespaces
   using namespace edm;
   using namespace reco;
   using namespace std;

//***********************
// for muon information *
//***********************


#include "DataFormats/PatCandidates/interface/Muon.h"

//***************************
// for electron information *
//***************************


#include "DataFormats/PatCandidates/interface/Electron.h"



// for photon information
#include "DataFormats/PatCandidates/interface/Photon.h"

//**********************
// for MET information *
#include "DataFormats/PatCandidates/interface/MET.h"

//**********************
// for jet information *
//**********************
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/Jet.h"


//*******************************
// for gen particle information *
//*******************************
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//********************************
// for particle flow information *
//********************************
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

//********************************
// for L1 information *
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"

//***************************
// class member declaration *
//***************************

//***********************************
// main analyzer class (EDAnalyzer) *
//***********************************

#ifndef CMSSW12plus
// Run 1 and 2
class NanoAnalyzerSingleEle : public edm::EDAnalyzer
#else
// Run 3
class NanoAnalyzerSingleEle : public edm::one::EDAnalyzer<>
#endif
{
public:
  explicit NanoAnalyzerSingleEle(const edm::ParameterSet&);
  ~NanoAnalyzerSingleEle();

  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

  // this is the place to define global variables and parameters 
  
  // declare global trigger variables
  //#include "NanoTrigger.h"

private:

  virtual void beginJob(const edm::ParameterSet& iConfig);
  virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iStp);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void endJob();

  void createBranch();

  void reset();
  
  // HLT config for reading the table and its associated process name
  //HLTConfigProvider hlt_cfg;   // superseded above


  EDGetTokenT<pat::ElectronCollection> electronToken_;
  EDGetTokenT<reco::VertexCollection> verticeToken_;
  EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken_;
  EDGetTokenT<edm::TriggerResults> triggerToken_;  
  EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobjectToken_;
  edm::EDGetTokenT<l1t::EGammaBxCollection> l1EG_;

  ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken;

  Handle<pat::ElectronCollection> electrons_;
  Handle< reco::VertexCollection > vertices_;
  Handle< std::vector<pat::PackedCandidate> > packedpfcandidates_   ;
  Handle< edm::TriggerResults> HLTtriggers_;
  Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;


/////////////////////////////////////////////////////////////////////////
////////////////////////// declare tree, file, //////////////////////////
/////////////////////////////////////////////////////////////////////////
  
  edm::Service<TFileService> fs;
  TTree* tree_;

  /// original nanoAOD ///
  int run;
  ULong64_t event;
  int luminosityBlock;
  int nvtx;

  float                JpsiKE_e1_pt      ;
  float                JpsiKE_e1_eta     ;
  float                JpsiKE_e1_phi     ;
  float                JpsiKE_e1_mass     ;
  int                  JpsiKE_e1_q      ;   
  float                JpsiKE_e1_vx       ;
  float                JpsiKE_e1_vy       ;
  float                JpsiKE_e1_vz       ;
  float                JpsiKE_e1_passMVA  ;
  float                JpsiKE_e1_bestL1pt  ;
  float                JpsiKE_e1_bestL1eta ;
  float                JpsiKE_e1_bestL1phi  ;
  float                JpsiKE_e1_bestL1Deta ;
  float                JpsiKE_e1_bestL1Dphi ;
  
  float                JpsiKE_e2_pt      ;
  float                JpsiKE_e2_eta     ;
  float                JpsiKE_e2_phi     ;
  float                JpsiKE_e2_mass     ;
  int                  JpsiKE_e2_alsotag ;
  int                  JpsiKE_e2_q   ;   
  float                JpsiKE_e2_vx       ;
  float                JpsiKE_e2_vy       ;
  float                JpsiKE_e2_vz       ;
  float                JpsiKE_e2_passMVA  ;
  float                JpsiKE_e2_bestL1pt  ;
  float                JpsiKE_e2_bestL1eta ;
  float                JpsiKE_e2_bestL1phi  ;
  float                JpsiKE_e2_bestL1Deta ;
  float                JpsiKE_e2_bestL1Dphi ;

  float                JpsiKE_e1_trgobj_pt      ;
  float                JpsiKE_e1_trgobj_eta     ;
  float                JpsiKE_e1_trgobj_phi     ;
  float                JpsiKE_e1_trgobj_mass     ;
  int                  JpsiKE_e1_trgobj_q      ;   
  float                JpsiKE_e1_trgobj_vx       ;
  float                JpsiKE_e1_trgobj_vy       ;
  float                JpsiKE_e1_trgobj_vz       ;

  float                JpsiKE_e2_trgobj_pt      ;
  float                JpsiKE_e2_trgobj_eta     ;
  float                JpsiKE_e2_trgobj_phi     ;
  float                JpsiKE_e2_trgobj_mass     ;
  int                  JpsiKE_e2_trgobj_q      ;   
  float                JpsiKE_e2_trgobj_vx       ;
  float                JpsiKE_e2_trgobj_vy       ;
  float                JpsiKE_e2_trgobj_vz       ;

  float                JpsiKE_Jpsi_pt      ;
  float                JpsiKE_Jpsi_nonfit_pt;
  float                JpsiKE_Jpsi_eta     ;
  float                JpsiKE_Jpsi_nonfit_eta;
  float                JpsiKE_Jpsi_phi     ;
  float                JpsiKE_Jpsi_nonfit_phi;
  float                JpsiKE_Jpsi_mass       ;
  float                JpsiKE_Jpsi_mass_nofit ;
  float                JpsiKE_Jpsi_vprob    ;
  float                JpsiKE_elesDr    ;

  std::vector<int> DoubleEle_fired{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  std::vector<float> dieleobj1_pt{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> dieleobj1_eta{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> dieleobj1_phi{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> dieleobj2_pt{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> dieleobj2_eta{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> dieleobj2_phi{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  
  std::vector<int> ele1_matchedDiEle{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<int> ele2_matchedDiEle{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<float> ele1_matchedDiEle_pt{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> ele2_matchedDiEle_pt{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> ele1_matchedDiEle_eta{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> ele2_matchedDiEle_eta{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> ele1_matchedDiEle_phi{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> ele2_matchedDiEle_phi{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};

  helperfunc aux;
  float chi = 0.;
  float ndf = 0.;

  TH1F * hist; 

}; // end of class member



//////////////////////////////////////////////////////////////////////////////
//                        set analysis loop parameters                      //
//////////////////////////////////////////////////////////////////////////////

NanoAnalyzerSingleEle::NanoAnalyzerSingleEle(const edm::ParameterSet& iConfig)
{
  electronToken_           = consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"));

  packedpfcandidatesToken_ = consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("packedpfcandidates")); 

  verticeToken_            = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  triggerToken_	      	   = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLT"));
  triggerobjectToken_	   = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerobjects"));

  l1EG_                    = consumes<l1t::EGammaBxCollection>(iConfig.getParameter<edm::InputTag>("l1EG"));

  ttkToken = esConsumes(edm::ESInputTag{"","TransientTrackBuilder"});

  hist = fs->make<TH1F>("cutflow", "cutflow", 10,0,10);
  tree_ = fs->make<TTree>( "tree", "tree" );

  createBranch();
} // end of constructor

NanoAnalyzerSingleEle::~NanoAnalyzerSingleEle() { }

///////////////////////////////////////////////////////////////////////////////
////////////// main analysis loop: method called for each event ///////////////
///////////////////////////////////////////////////////////////////////////////

void NanoAnalyzerSingleEle::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  reset();
  
  iEvent.getByToken(triggerToken_, HLTtriggers_);
  iEvent.getByToken(triggerobjectToken_ , triggerObjects);
  iEvent.getByToken(verticeToken_, vertices_     ); 

  run = (iEvent.id()).run();
  event = (iEvent.id()).event();
  luminosityBlock = (iEvent.id()).luminosityBlock();
  
  nvtx = vertices_->size();

  hist->Fill(1);

  // pT thresholds for di-ele trigger
  std::vector<std::string> pt_thr_v_string{"10", "9p5", "9", "8p5", "8", "7p5", "7", "6p5", "6", "5p5", "5", "4p5", "4"};

  // Loop over HLT paths
  bool isTriggered = false;
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*HLTtriggers_);
  std::string singleeleTriggerName="";
  // 
  for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {
    
    // Check if the single ele trigger path fired
    if(trigNames.triggerName(i).find("HLT_SingleEle8_v")!= std::string::npos){
      //if(trigNames.triggerName(i).find("HLT_SingleEle8_SingleEGL1_v")!= std::string::npos){
      if(HLTtriggers_->accept(i)){
	isTriggered = true;
	singleeleTriggerName=trigNames.triggerName(i);  
      }
    }

    // Check if any of the di-ele trigger paths fired
    for(int j=0; j<int(pt_thr_v_string.size()); j++){
      std::string pt_thr_string = pt_thr_v_string[j];
      std::string DoubleEleTrigName = "HLT_DoubleEle"+pt_thr_string+"_eta1p22_mMax6";
      if(trigNames.triggerName(i).find(DoubleEleTrigName)!= std::string::npos){
	if(HLTtriggers_->accept(i)){
	  DoubleEle_fired[j]=1;
	  // std::cout << "Event = " << event << " : DoubleEle_fired[j] fired for j = " << j << std::endl;
	}
      }
    }
  }

  // Our wanted (singleele) path fired
  if(!isTriggered) return; 
  hist->Fill(2);


  // Take needed collections
  iEvent.getByToken(electronToken_ , electrons_    );

  // Output collections
  std::vector<pat::Electron> electroncollection;                  // collection with offline electrons passing minimal selection
  std::vector<pat::TriggerObjectStandAlone> trg_obj_collection;   // collection with HLT candidate with best match with offline electrons 
                                                                  // (passing matching criteria) - at most one per electron

  std::vector<int> electronmatched;                               // integer with the position in the trg_obj_collection of the HLT object matched to ele
  electroncollection.clear();
  trg_obj_collection.clear();
  electronmatched.clear();

  // std::cout << "event = " << event << ", electrons_->size() = " << electrons_->size() << std::endl;

  // Offline electrons
  for(size_t ielectron = 0; ielectron < electrons_->size(); ++ ielectron){
    
    const pat::Electron & electron = (*electrons_)[ielectron];

    // std::cout << "event = " << event << ", ielectron = " << ielectron << ", electron.pt() =  " << electron.pt() << std::endl;

    // Offline cuts
    if ( electron.pt() < 2.5 ) continue;                 
    if ( fabs(electron.eta()) > 2.4 ) continue;    
    if (!electron.passConversionVeto() ) continue; 
    const reco::GsfTrackRef gsfTrk = electron.gsfTrack();
    if(!gsfTrk.isNonnull()) continue;

    // std::cout << "event = " << event << ", ielectron = " << ielectron << ", electron.pt() =  " << electron.pt() << std::endl;

    /// Trigger matching 
    bool trigObjMatchEle = false;        // this offline electron matches a HLT ele-candidate
    pat::TriggerObjectStandAlone best_match_obj;      // this offline electron matches a HLT ele-candidate and this is the best matched candidate
    Float_t best_match_dR = 9999999;                  // this offline electron matches a HLT ele-candidate and this is the best match DR

    // Loop over trigger objects matching the single ele path
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
      
      // consider only objects which match the single ele path
      obj.unpackPathNames(trigNames);
      obj.unpackFilterLabels(iEvent, *HLTtriggers_);
      std::vector<std::string> pathNamesAll = obj.pathNames(false);
      bool isPathExist = false;
      for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
      	if(pathNamesAll[h]==singleeleTriggerName) isPathExist = true;
      }
      if(!isPathExist) continue;

      int eleObjNumber = -1;
      for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){	
	if(obj.filterLabels()[hh].find("hltSingleEle8ValidHitsFilter") != std::string::npos) {
	  //if(obj.filterLabels()[hh].find("hltSingleEle8SingleEGL1ValidHitsFilter") != std::string::npos) {
	  eleObjNumber = hh;
	}
      }
        
      // here HLT obj vs reco electron candidates
      TVector3 eleTV3, objTV3;
      eleTV3.SetPtEtaPhi( electron.pt(), electron.eta(), electron.phi() );
      objTV3.SetPtEtaPhi( obj.pt(), obj.eta(), obj.phi() );
      Float_t deltaEta = fabs(obj.eta() - electron.eta());
      Float_t deltaPhi = fabs(eleTV3.DeltaPhi(objTV3));
      Float_t deltaR   = fabs(eleTV3.DeltaR(objTV3));
    
      // here HLT-ele candidates
      if (eleObjNumber>=0) {
	//std::cout<< "DeltaEta = " << deltaEta << ", deltaPhi = " << deltaPhi << ", deltaR = " << deltaR << std::endl;
	//std::cout << "This is an electron HLT candidate" << endl;
	if(deltaPhi < 0.2 && deltaEta < 0.07){ 
	  trigObjMatchEle = true;
	  if (deltaR < best_match_dR){
	    best_match_dR = deltaR;
	    best_match_obj = obj;
	  }
	  //std::cout << "This object is matched eith ele: deltaPhi = " << deltaPhi << ", deltaEta = " << deltaEta << ", best_match_dR = " << best_match_dR << std::endl;
	  //std::cout << "Offline: " << electron.pt() << " " << electron.eta() << " " << electron.phi() << std::endl;
	  //std::cout << "HLT: " << obj.pt() << " " << obj.eta() << " " << obj.phi() << std::endl;
	}
      }
    
      // std::cout << "event = " << event << " : eleObjNumber = " << eleObjNumber << std::endl;
      // std::cout << "trigObjMatchEle = " << trigObjMatchEle << std::endl;
      // std::cout << "In objects loop, best_match_obj => " << best_match_obj.pt() << " " << best_match_obj.eta() << " " << best_match_obj.phi() << std::endl;
      
    }  // Loop over trigger objects
  
    // std::cout << "After objects loop, best_match_obj => " << best_match_obj.pt() << " " << best_match_obj.eta() << " " << best_match_obj.phi() << std::endl;


    // This is to further process 
    electroncollection.push_back(electron);
    // ele/HLT
    if(trigObjMatchEle){ 
      trg_obj_collection.push_back(best_match_obj);
      electronmatched.push_back(int(trg_obj_collection.size())-1);
    }
    else{
      electronmatched.push_back(-999);	
    }

  } // Loop over offline electrons

  // std::cout << "event = " << event << ", electroncollection.size() = " << electroncollection.size() << ", electronmatched.size() = " << electronmatched.size() << std::endl;

  if (electroncollection.size() < 2) return; 
  hist->Fill(3);

  // NB: electroncollection and electronmatched must have the same size
  if (electroncollection.size() != electronmatched.size() ) return;
  hist->Fill(4);

  // Tools for tracks
  const TransientTrackBuilder* builder = &iSetup.getData(ttkToken);

  // Prepare ele pairs
  float jpsi_max_pt = -1;
  int mcidx_e1 = -1;
  int mcidx_e2 = -1;
  int mcidx_trgobj1 = -1;
  int mcidx_trgobj2 = -1;
  TLorentzVector jpsi_tlv_highest;

  for(int ie = 0; ie < (int)electroncollection.size(); ie++){
    for(int je = ie+1; je < (int)electroncollection.size(); je++){
      
      // at least 1 electron of the jpsi pair should match a trigger ele  
      if ( electronmatched[ie]<0 && electronmatched[je]<0 ) continue;
      
      // std::cout << "ie = " << ie << ", electronmatched[ie] = " << electronmatched[ie] << ", je = " << je << ", electronmatched[je] = " << electronmatched[je] << std::endl; 

      // match with trigger ele
      int this_mcidx_trgobj1 = -1;
      int this_mcidx_trgobj2 = -1;
      if (electronmatched[ie]>-999){
	this_mcidx_trgobj1 = electronmatched[ie];
	// std::cout << "mcidx_trgobj1 = " << this_mcidx_trgobj1 << std::endl; 
      }
      if (electronmatched[je]>-999){
	this_mcidx_trgobj2 = electronmatched[je];
	// std::cout << "mcidx_trgobj2 = " << this_mcidx_trgobj2 << std::endl; 
      }
      
      const pat::Electron e1 = electroncollection[ie];
      const pat::Electron e2 = electroncollection[je];

      /*
	std::cout << "Ele1 " << ie << ", pt = " << e1.pt() << ", eta = " << e1.eta() << ", phi = " << e1.phi() << ", id = " << e1.pdgId() << std::endl;
	
	if (this_mcidx_trgobj1>=0) std::cout << "TrgObj1 pt = " << trg_obj_collection[this_mcidx_trgobj1].pt() << ", eta = " << trg_obj_collection[this_mcidx_trgobj1].eta() << ", phi = " << trg_obj_collection[this_mcidx_trgobj1].phi() << std::endl; 
	else 
	  std::cout << "TrgObj1 not found" << std::endl;

	std::cout << "Ele2 " << je << ", pt = " << e2.pt() << ", eta = " << e2.eta() << ", phi = " << e2.phi() << ", id = " << e2.pdgId() << std::endl;

	if (this_mcidx_trgobj2>=0) std::cout << "TrgObj2 pt = " << trg_obj_collection[this_mcidx_trgobj2].pt() << ", eta = " << trg_obj_collection[this_mcidx_trgobj2].eta() << ", phi = " << trg_obj_collection[this_mcidx_trgobj2].phi() << std::endl; 
	else 
	  std::cout << "TrgObj2 not found" << std::endl;
      */

      TLorentzVector tlv_e1;
      TLorentzVector tlv_e2;
      tlv_e1.SetPtEtaPhiM(e1.pt(), e1.eta(), e1.phi(), aux.mass_electron);
      tlv_e2.SetPtEtaPhiM(e2.pt(), e2.eta(), e2.phi(), aux.mass_electron);

      TLorentzVector tlv_jpsi = (tlv_e1 + tlv_e2);
      float jpsi_mass = tlv_jpsi.M();
      float jpsi_pt = tlv_jpsi.Pt();

      // std::cout << "event = " << event << ", jpsi_mass = " << jpsi_mass << std::endl;

      if (e1.charge() + e2.charge() !=0) continue;
      if (jpsi_mass < 2.0) continue; 
      if (jpsi_mass > 4.0) continue;
      // std::cout << "event = " << event << ": jpsi_mass = " << jpsi_mass << ", jpsi_max_pt = " << jpsi_max_pt << ", jpsi_pt = " << jpsi_pt << std::endl;

      if(jpsi_max_pt < jpsi_pt){
	jpsi_max_pt = jpsi_pt;
	mcidx_e1 = ie;
	mcidx_e2 = je;
	mcidx_trgobj1 = this_mcidx_trgobj1;
	mcidx_trgobj2 = this_mcidx_trgobj2;
	jpsi_tlv_highest = tlv_jpsi;
      }
      // std::cout << "jpsi_mass = " << jpsi_mass << ", jpsi_max_pt = " << jpsi_max_pt << ", jpsi_pt = " << jpsi_pt << std::endl;      
      // std::cout << "event = " << event << "; In the loop: mcidx_e1 = " << mcidx_e1 << ", mcidx_e2 = " << mcidx_e2 << std::endl;
    }
  }

  // std::cout << "event = " << event << ", jpsi_max_pt = " << jpsi_max_pt << std::endl;

  // At least 1 reco J/psi
  if(jpsi_max_pt == -1) return;  
  hist->Fill(5);

  // std::cout << "event = " << event << ", jpsi_max_pt = " << jpsi_max_pt << " => passed " << std::endl;


  // Inputs to kin fit
  const reco::GsfTrackRef track1_electron = electroncollection[mcidx_e1].gsfTrack();
  const reco::GsfTrackRef track2_electron = electroncollection[mcidx_e2].gsfTrack();
  reco::TransientTrack tt1_electron = (*builder).build(track1_electron);
  reco::TransientTrack tt2_electron = (*builder).build(track2_electron);
  KinematicParticleFactoryFromTransientTrack pFactory;
  std::vector<RefCountedKinematicParticle> electronParticles;
  electronParticles.push_back(pFactory.particle(tt1_electron, aux.electron_mass, chi, ndf, aux.electron_sigma));
  electronParticles.push_back(pFactory.particle(tt2_electron, aux.electron_mass, chi, ndf, aux.electron_sigma));

  // Kinematic fit of the two electrons to a common vtx
  RefCountedKinematicParticle jpsi_part;
  RefCountedKinematicVertex jpsi_vertex;
  RefCountedKinematicTree jpTree;
  Bool_t jpsifit_flag;
  std::tie(jpsifit_flag, jpsi_part, jpsi_vertex, jpTree) = aux.KinematicFit(electronParticles, -1, -1);

  // Successfull kin fit
  if(!jpsifit_flag) return;
  hist->Fill(6);
  
  
  // Loop over di-ele paths, when fired
  for(int j=0; j<int(pt_thr_v_string.size()); j++){

    if (DoubleEle_fired[j]==1) {   

      std::string pt_thr_string = pt_thr_v_string[j];
      std::string DoubleEleTrigName = "HLT_DoubleEle"+pt_thr_string+"_eta1p22_mMax6";
      std::string DoubleEleObjName  = "hltDoubleEle"+pt_thr_string+"eta1p22mMax6ValidHitsFilter"; 
      
      // Loop over trigger objects matching the diele path 
      for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
	
	// check if this obj comes from the wanted di-ele path
	obj.unpackPathNames(trigNames);
	obj.unpackFilterLabels(iEvent, *HLTtriggers_);
	std::vector<std::string> pathNamesAll = obj.pathNames(false);
	bool isPathExist = false;
	for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
	  // std::cout << "event " << event << " => path " << h << " => " << pathNamesAll[h] << std::endl;
	  if(pathNamesAll[h].find(DoubleEleTrigName)!= std::string::npos) {
	    isPathExist = true;   
	    // std::cout << "event " << event << " => found" << std::endl;
	  }
	}
	if(!isPathExist) continue;
	
	// std::cout << "ok, this object matches our diele path " << DoubleEleTrigName << std::endl;
	
	// check if the object is from the correct filter
	for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){	
	  // std::cout << "event " << event << " => hh = " << hh << ", obj.filterLabels()[hh] = " << obj.filterLabels()[hh] << std::endl;
	  if(obj.filterLabels()[hh].find(DoubleEleObjName) != std::string::npos) {
	    if (dieleobj1_pt[j]<0) {
	      dieleobj1_pt[j]  = obj.pt();
	      dieleobj1_eta[j] = obj.eta();
	      dieleobj1_phi[j] = obj.phi();
	    } else if (dieleobj2_pt[j]<0) {
	      dieleobj2_pt[j]  = obj.pt();
	      dieleobj2_eta[j] = obj.eta();
	      dieleobj2_phi[j] = obj.phi();
	    }
	  }
	}
      }
    } // path fired
  } // Loop over paths


  // Swap 1<->2 so that: 1 is always tag, 2 is always probe, and additional variable tells if 2 is also tag
  int mcidx_tag = -99;
  int mcidx_pro = -99;
  int mcidx_trgobj_tag = -99;
  int mcidx_trgobj_pro = -99;
  int is_probe_also_tag = -99;

  if (mcidx_trgobj1>=0) {
    mcidx_pro = mcidx_e2;
    mcidx_trgobj_pro = mcidx_trgobj2;
    mcidx_tag = mcidx_e1;
    mcidx_trgobj_tag = mcidx_trgobj1;
    if (mcidx_trgobj2<0)  is_probe_also_tag = 0;
    if (mcidx_trgobj2>=0) is_probe_also_tag = 1;
  } else {  
    mcidx_pro = mcidx_e1;
    mcidx_trgobj_pro = mcidx_trgobj1;
    mcidx_tag = mcidx_e2;
    mcidx_trgobj_tag = mcidx_trgobj2;
    is_probe_also_tag = 0;
  }

  mcidx_e1 = mcidx_tag;
  mcidx_e2 = mcidx_pro;
  mcidx_trgobj1 = mcidx_trgobj_tag;
  mcidx_trgobj2 = mcidx_trgobj_pro;

  // Offline electrons
  float reco_e1_pt  = electroncollection[mcidx_e1].pt(); 
  float reco_e1_eta = electroncollection[mcidx_e1].eta(); 
  float reco_e1_phi = electroncollection[mcidx_e1].phi(); 
  float reco_e2_pt  = electroncollection[mcidx_e2].pt(); 
  float reco_e2_eta = electroncollection[mcidx_e2].eta(); 
  float reco_e2_phi = electroncollection[mcidx_e2].phi(); 
  TVector3 ele1TV3, ele2TV3;
  ele1TV3.SetPtEtaPhi( reco_e1_pt, reco_e1_eta, reco_e1_phi );
  ele2TV3.SetPtEtaPhi( reco_e2_pt, reco_e2_eta, reco_e2_phi );
  //std::cout << "Offline1: " << reco_e1_pt << " " << reco_e1_eta << " " << reco_e1_phi << std::endl;
  //std::cout << "Offline2: " << reco_e2_pt << " " << reco_e2_eta << " " << reco_e2_phi << std::endl;


  // Match HLT / offline
  for(int j=0; j<int(pt_thr_v_string.size()); j++){
    
    if (dieleobj1_pt[j]>=0) {
      float obj1_pt  = dieleobj1_pt[j];
      float obj1_eta = dieleobj1_eta[j];
      float obj1_phi = dieleobj1_phi[j];
      TVector3 obj1TV3;
      obj1TV3.SetPtEtaPhi( obj1_pt, obj1_eta, obj1_phi );
      Float_t deltaEta11 = fabs(reco_e1_eta - obj1_eta);
      Float_t deltaPhi11 = fabs(ele1TV3.DeltaPhi(obj1TV3));
      Float_t deltaEta21 = fabs(reco_e2_eta - obj1_eta);
      Float_t deltaPhi21 = fabs(ele2TV3.DeltaPhi(obj1TV3));
      if (deltaEta11<0.07 && deltaPhi11<0.2) {
	ele1_matchedDiEle[j] = 1;      
	ele1_matchedDiEle_pt[j]  = dieleobj1_pt[j];
	ele1_matchedDiEle_eta[j] = dieleobj1_eta[j];
	ele1_matchedDiEle_phi[j] = dieleobj1_phi[j];
	//std::cout << "Ok match offline1 / online1 " << std::endl;
      }
      if (deltaEta21<0.07 && deltaPhi21<0.2) {
	ele2_matchedDiEle[j] = 1; 
	ele2_matchedDiEle_pt[j]  = dieleobj1_pt[j];
	ele2_matchedDiEle_eta[j] = dieleobj1_eta[j];
	ele2_matchedDiEle_phi[j] = dieleobj1_phi[j];
	//std::cout << "Ok match offline2 / online1 " << std::endl;
      }
    }

    if (dieleobj2_pt[j]>=0) {
      float obj2_pt  = dieleobj2_pt[j];
      float obj2_eta = dieleobj2_eta[j];
      float obj2_phi = dieleobj2_phi[j];
      TVector3 obj2TV3;
      obj2TV3.SetPtEtaPhi( obj2_pt, obj2_eta, obj2_phi );
      Float_t deltaEta12 = fabs(reco_e1_eta - obj2_eta);
      Float_t deltaPhi12 = fabs(ele1TV3.DeltaPhi(obj2TV3));
      Float_t deltaEta22 = fabs(reco_e2_eta - obj2_eta);
      Float_t deltaPhi22 = fabs(ele2TV3.DeltaPhi(obj2TV3));
      if (deltaEta12<0.07 && deltaPhi12<0.2) {
	ele1_matchedDiEle[j] = 1; 
	ele1_matchedDiEle_pt[j]  = dieleobj2_pt[j];
	ele1_matchedDiEle_eta[j] = dieleobj2_eta[j];
	ele1_matchedDiEle_phi[j] = dieleobj2_phi[j];
	//std::cout << "Ok match offline1 / online2 " << std::endl;
      }
      if (deltaEta22<0.07 && deltaPhi22<0.2) {
	ele2_matchedDiEle[j] = 1; 
	ele2_matchedDiEle_pt[j]  = dieleobj2_pt[j];
	ele2_matchedDiEle_eta[j] = dieleobj2_eta[j];
	ele2_matchedDiEle_phi[j] = dieleobj2_phi[j];
	//std::cout << "Ok match offline2 / online2 " << std::endl;
      }
    }
  }

  // Match L1 / offline
  float bestMatchE1_eta  = -99.;
  float bestMatchE1_phi  = -99.;
  float bestMatchE1_pt   = -99.;
  float bestMatchE1_Deta = -99.;
  float bestMatchE1_Dphi = -99.;
  float bestMatchE2_eta  = -99.;
  float bestMatchE2_phi  = -99.;
  float bestMatchE2_pt   = -99.;
  float bestMatchE2_Deta = -99.;
  float bestMatchE2_Dphi = -99.;
  float bestMatchDrE1 = 999;
  float bestMatchDrE2 = 999;
  const auto &l1EG = iEvent.get(l1EG_);
  for (l1t::EGammaBxCollection::const_iterator it = l1EG.begin(0); it != l1EG.end(0); it++) {
    pat::TriggerObjectStandAlone l1obj(it->p4());
    
    TVector3 l1objTV3;
    l1objTV3.SetPtEtaPhi( it->pt(), it->eta(), it->phi() );
    Float_t deltaRE1 = fabs(ele1TV3.DeltaR(l1objTV3));
    Float_t deltaRE2 = fabs(ele2TV3.DeltaR(l1objTV3));

    if (deltaRE1<bestMatchDrE1 && deltaRE1<0.25){
      bestMatchDrE1 = deltaRE1;
      bestMatchE1_eta  = it->eta(); 
      bestMatchE1_phi  = it->phi(); 
      bestMatchE1_pt   = it->pt();    
      bestMatchE1_Deta = fabs(it->eta() - electroncollection[mcidx_e1].eta());
      bestMatchE1_Dphi = fabs(ele1TV3.DeltaPhi(l1objTV3));
    }
    if (deltaRE2<bestMatchDrE2 && deltaRE2<0.25){
      bestMatchDrE2 = deltaRE2;
      bestMatchE2_eta  = it->eta(); 
      bestMatchE2_phi  = it->phi(); 
      bestMatchE2_pt   = it->pt(); 
      bestMatchE2_Deta = fabs(it->eta() - electroncollection[mcidx_e2].eta());
      bestMatchE2_Dphi = fabs(ele2TV3.DeltaPhi(l1objTV3));
    }
  }

  // Infos about JPsi candidate
  JpsiKE_e1_pt   = electroncollection[mcidx_e1].pt();
  JpsiKE_e1_eta  = electroncollection[mcidx_e1].eta();
  JpsiKE_e1_phi  = electroncollection[mcidx_e1].phi();
  JpsiKE_e1_mass = electroncollection[mcidx_e1].mass();
  JpsiKE_e1_q    = electroncollection[mcidx_e1].charge();
  JpsiKE_e1_vx   = electroncollection[mcidx_e1].vx();
  JpsiKE_e1_vy   = electroncollection[mcidx_e1].vy();
  JpsiKE_e1_vz   = electroncollection[mcidx_e1].vz();
  JpsiKE_e1_passMVA = electroncollection[mcidx_e1].electronID("mvaEleID-Fall17-noIso-V2-wpLoose");
  JpsiKE_e1_bestL1pt   = bestMatchE1_pt;
  JpsiKE_e1_bestL1eta  = bestMatchE1_eta;
  JpsiKE_e1_bestL1phi  = bestMatchE1_phi;
  JpsiKE_e1_bestL1Deta = bestMatchE1_Deta;  
  JpsiKE_e1_bestL1Dphi = bestMatchE1_Dphi;  

  JpsiKE_e2_pt   = electroncollection[mcidx_e2].pt();
  JpsiKE_e2_eta  = electroncollection[mcidx_e2].eta();
  JpsiKE_e2_phi  = electroncollection[mcidx_e2].phi();
  JpsiKE_e2_mass = electroncollection[mcidx_e2].mass();
  JpsiKE_e2_alsotag = is_probe_also_tag;
  JpsiKE_e2_q    = electroncollection[mcidx_e2].charge();
  JpsiKE_e2_vx   = electroncollection[mcidx_e2].vx();
  JpsiKE_e2_vy   = electroncollection[mcidx_e2].vy();
  JpsiKE_e2_vz   = electroncollection[mcidx_e2].vz();
  JpsiKE_e2_passMVA = electroncollection[mcidx_e2].electronID("mvaEleID-Fall17-noIso-V2-wpLoose");
  JpsiKE_e2_bestL1pt   = bestMatchE2_pt;
  JpsiKE_e2_bestL1eta  = bestMatchE2_eta;
  JpsiKE_e2_bestL1phi  = bestMatchE2_phi;
  JpsiKE_e2_bestL1Deta = bestMatchE2_Deta;  
  JpsiKE_e2_bestL1Dphi = bestMatchE2_Dphi;  

  if (mcidx_trgobj1>=0) {
    JpsiKE_e1_trgobj_pt   = trg_obj_collection[mcidx_trgobj1].pt();
    JpsiKE_e1_trgobj_eta  = trg_obj_collection[mcidx_trgobj1].eta();
    JpsiKE_e1_trgobj_phi  = trg_obj_collection[mcidx_trgobj1].phi();
    JpsiKE_e1_trgobj_mass = trg_obj_collection[mcidx_trgobj1].mass();
    JpsiKE_e1_trgobj_q    = trg_obj_collection[mcidx_trgobj1].charge();
    JpsiKE_e1_trgobj_vx   = trg_obj_collection[mcidx_trgobj1].vx();
    JpsiKE_e1_trgobj_vy   = trg_obj_collection[mcidx_trgobj1].vy();
    JpsiKE_e1_trgobj_vz   = trg_obj_collection[mcidx_trgobj1].vz();
  }

  if (mcidx_trgobj2>=0) {
    JpsiKE_e2_trgobj_pt   = trg_obj_collection[mcidx_trgobj2].pt();
    JpsiKE_e2_trgobj_eta  = trg_obj_collection[mcidx_trgobj2].eta();
    JpsiKE_e2_trgobj_phi  = trg_obj_collection[mcidx_trgobj2].phi();
    JpsiKE_e2_trgobj_mass = trg_obj_collection[mcidx_trgobj2].mass();
    JpsiKE_e2_trgobj_q    = trg_obj_collection[mcidx_trgobj2].charge();
    JpsiKE_e2_trgobj_vx   = trg_obj_collection[mcidx_trgobj2].vx();
    JpsiKE_e2_trgobj_vy   = trg_obj_collection[mcidx_trgobj2].vy();
    JpsiKE_e2_trgobj_vz   = trg_obj_collection[mcidx_trgobj2].vz();
  }

  JpsiKE_Jpsi_pt  = jpsi_part->currentState().globalMomentum().perp();
  JpsiKE_Jpsi_eta = jpsi_part->currentState().globalMomentum().eta();
  JpsiKE_Jpsi_phi = jpsi_part->currentState().globalMomentum().phi();
  JpsiKE_Jpsi_nonfit_pt  = jpsi_tlv_highest.Pt();
  JpsiKE_Jpsi_nonfit_eta = jpsi_tlv_highest.Eta();
  JpsiKE_Jpsi_nonfit_phi = jpsi_tlv_highest.Phi();
  JpsiKE_Jpsi_mass = jpsi_part->currentState().mass();
  JpsiKE_Jpsi_mass_nofit = jpsi_tlv_highest.M();
  JpsiKE_Jpsi_vprob = TMath::Prob(jpsi_part->chiSquared(), jpsi_part->degreesOfFreedom());
  JpsiKE_elesDr = ele1TV3.DeltaR(ele2TV3); 

  tree_->Fill();

  return;

  } //NanoAnalyzerSingleEle::analyze ends


//**************************************************
//************* additional methods *****************
//**************************************************

// ------------ method called once each job just before starting event loop  ------------
void NanoAnalyzerSingleEle::beginJob(const edm::ParameterSet& iConfig)
{ }

// ------------ method called once each job just before starting run  ------------
void NanoAnalyzerSingleEle::beginRun(const edm::Run &iRun, const edm::EventSetup &iStp)
{ }

void NanoAnalyzerSingleEle::fillDescriptions(edm::ConfigurationDescriptions & descriptions) 
{ }

// ------------ method called when ending the processing of a run  ------------  Qun below
void NanoAnalyzerSingleEle::endRun(edm::Run const&, edm::EventSetup const&)
{ }

void 
NanoAnalyzerSingleEle::endJob()
{ }

//include methods for special trigger variables
//#include "NanoTrigger.cc.forinclude"

//define this as a plug-in

// branch title creation
void NanoAnalyzerSingleEle::createBranch() { 

  tree_->Branch("run", &run, "run/i");
  tree_->Branch("event", &event, "event/l");
  tree_->Branch("luminosityBlock", &luminosityBlock, "luminosityBlock/i");
  tree_->Branch("nvtx", &nvtx, "nvtx/i");

  tree_->Branch("JpsiKE_e1_pt", &JpsiKE_e1_pt );
  tree_->Branch("JpsiKE_e1_eta", &JpsiKE_e1_eta );
  tree_->Branch("JpsiKE_e1_phi", &JpsiKE_e1_phi );
  tree_->Branch("JpsiKE_e1_mass", &JpsiKE_e1_mass );
  tree_->Branch("JpsiKE_e1_q", &JpsiKE_e1_q );
  tree_->Branch("JpsiKE_e1_vx"   , &JpsiKE_e1_vx    );
  tree_->Branch("JpsiKE_e1_vy"   , &JpsiKE_e1_vy    );
  tree_->Branch("JpsiKE_e1_vz"   , &JpsiKE_e1_vz    );
  tree_->Branch("JpsiKE_e1_passMVA"   , &JpsiKE_e1_passMVA    );
  tree_->Branch("JpsiKE_e1_bestL1pt",   &JpsiKE_e1_bestL1pt   );
  tree_->Branch("JpsiKE_e1_bestL1eta" , &JpsiKE_e1_bestL1eta  );
  tree_->Branch("JpsiKE_e1_bestL1phi" , &JpsiKE_e1_bestL1phi  );
  tree_->Branch("JpsiKE_e1_bestL1Deta", &JpsiKE_e1_bestL1Deta );
  tree_->Branch("JpsiKE_e1_bestL1Dphi", &JpsiKE_e1_bestL1Dphi );

  tree_->Branch("JpsiKE_e1_Ele10_match",  &ele1_matchedDiEle[0] );
  tree_->Branch("JpsiKE_e1_Ele9p5_match", &ele1_matchedDiEle[1] );
  tree_->Branch("JpsiKE_e1_Ele9_match",   &ele1_matchedDiEle[2] );
  tree_->Branch("JpsiKE_e1_Ele8p5_match", &ele1_matchedDiEle[3] );
  tree_->Branch("JpsiKE_e1_Ele8_match",   &ele1_matchedDiEle[4] );
  tree_->Branch("JpsiKE_e1_Ele7p5_match", &ele1_matchedDiEle[5] );
  tree_->Branch("JpsiKE_e1_Ele7_match",   &ele1_matchedDiEle[6] );
  tree_->Branch("JpsiKE_e1_Ele6p5_match", &ele1_matchedDiEle[7] );
  tree_->Branch("JpsiKE_e1_Ele6_match",   &ele1_matchedDiEle[8] );
  tree_->Branch("JpsiKE_e1_Ele5p5_match", &ele1_matchedDiEle[9] );
  tree_->Branch("JpsiKE_e1_Ele5_match",   &ele1_matchedDiEle[10] );
  tree_->Branch("JpsiKE_e1_Ele4p5_match", &ele1_matchedDiEle[11] );
  tree_->Branch("JpsiKE_e1_Ele4_match",   &ele1_matchedDiEle[12] );

  tree_->Branch("JpsiKE_e1_Ele10_pt",  &ele1_matchedDiEle_pt[0] );
  tree_->Branch("JpsiKE_e1_Ele9p5_pt", &ele1_matchedDiEle_pt[1] );
  tree_->Branch("JpsiKE_e1_Ele9_pt",   &ele1_matchedDiEle_pt[2] );
  tree_->Branch("JpsiKE_e1_Ele8p5_pt", &ele1_matchedDiEle_pt[3] );
  tree_->Branch("JpsiKE_e1_Ele8_pt",   &ele1_matchedDiEle_pt[4] );
  tree_->Branch("JpsiKE_e1_Ele7p5_pt", &ele1_matchedDiEle_pt[5] );
  tree_->Branch("JpsiKE_e1_Ele7_pt",   &ele1_matchedDiEle_pt[6] );
  tree_->Branch("JpsiKE_e1_Ele6p5_pt", &ele1_matchedDiEle_pt[7] );
  tree_->Branch("JpsiKE_e1_Ele6_pt",   &ele1_matchedDiEle_pt[8] );
  tree_->Branch("JpsiKE_e1_Ele5p5_pt", &ele1_matchedDiEle_pt[9] );
  tree_->Branch("JpsiKE_e1_Ele5_pt",   &ele1_matchedDiEle_pt[10] );
  tree_->Branch("JpsiKE_e1_Ele4p5_pt", &ele1_matchedDiEle_pt[11] );
  tree_->Branch("JpsiKE_e1_Ele4_pt",   &ele1_matchedDiEle_pt[12] );

  tree_->Branch("JpsiKE_e1_Ele10_eta",  &ele1_matchedDiEle_eta[0] );
  tree_->Branch("JpsiKE_e1_Ele9p5_eta", &ele1_matchedDiEle_eta[1] );
  tree_->Branch("JpsiKE_e1_Ele9_eta",   &ele1_matchedDiEle_eta[2] );
  tree_->Branch("JpsiKE_e1_Ele8p5_eta", &ele1_matchedDiEle_eta[3] );
  tree_->Branch("JpsiKE_e1_Ele8_eta",   &ele1_matchedDiEle_eta[4] );
  tree_->Branch("JpsiKE_e1_Ele7p5_eta", &ele1_matchedDiEle_eta[5] );
  tree_->Branch("JpsiKE_e1_Ele7_eta",   &ele1_matchedDiEle_eta[6] );
  tree_->Branch("JpsiKE_e1_Ele6p5_eta", &ele1_matchedDiEle_eta[7] );
  tree_->Branch("JpsiKE_e1_Ele6_eta",   &ele1_matchedDiEle_eta[8] );
  tree_->Branch("JpsiKE_e1_Ele5p5_eta", &ele1_matchedDiEle_eta[9] );
  tree_->Branch("JpsiKE_e1_Ele5_eta",   &ele1_matchedDiEle_eta[10] );
  tree_->Branch("JpsiKE_e1_Ele4p5_eta", &ele1_matchedDiEle_eta[11] );
  tree_->Branch("JpsiKE_e1_Ele4_eta",   &ele1_matchedDiEle_eta[12] );

  tree_->Branch("JpsiKE_e1_Ele10_phi",  &ele1_matchedDiEle_phi[0] );
  tree_->Branch("JpsiKE_e1_Ele9p5_phi", &ele1_matchedDiEle_phi[1] );
  tree_->Branch("JpsiKE_e1_Ele9_phi",   &ele1_matchedDiEle_phi[2] );
  tree_->Branch("JpsiKE_e1_Ele8p5_phi", &ele1_matchedDiEle_phi[3] );
  tree_->Branch("JpsiKE_e1_Ele8_phi",   &ele1_matchedDiEle_phi[4] );
  tree_->Branch("JpsiKE_e1_Ele7p5_phi", &ele1_matchedDiEle_phi[5] );
  tree_->Branch("JpsiKE_e1_Ele7_phi",   &ele1_matchedDiEle_phi[6] );
  tree_->Branch("JpsiKE_e1_Ele6p5_phi", &ele1_matchedDiEle_phi[7] );
  tree_->Branch("JpsiKE_e1_Ele6_phi",   &ele1_matchedDiEle_phi[8] );
  tree_->Branch("JpsiKE_e1_Ele5p5_phi", &ele1_matchedDiEle_phi[9] );
  tree_->Branch("JpsiKE_e1_Ele5_phi",   &ele1_matchedDiEle_phi[10] );
  tree_->Branch("JpsiKE_e1_Ele4p5_phi", &ele1_matchedDiEle_phi[11] );
  tree_->Branch("JpsiKE_e1_Ele4_phi",   &ele1_matchedDiEle_phi[12] );

  tree_->Branch("JpsiKE_e2_pt", &JpsiKE_e2_pt );
  tree_->Branch("JpsiKE_e2_eta", &JpsiKE_e2_eta );
  tree_->Branch("JpsiKE_e2_phi", &JpsiKE_e2_phi );
  tree_->Branch("JpsiKE_e2_mass", &JpsiKE_e2_mass );
  tree_->Branch("JpsiKE_e2_alsotag", &JpsiKE_e2_alsotag );
  tree_->Branch("JpsiKE_e2_q", &JpsiKE_e2_q );
  tree_->Branch("JpsiKE_e2_vx"   , &JpsiKE_e2_vx    );
  tree_->Branch("JpsiKE_e2_vy"   , &JpsiKE_e2_vy    );
  tree_->Branch("JpsiKE_e2_vz"   , &JpsiKE_e2_vz    );
  tree_->Branch("JpsiKE_e2_passMVA"   , &JpsiKE_e2_passMVA    );
  tree_->Branch("JpsiKE_e2_bestL1pt",   &JpsiKE_e2_bestL1pt );
  tree_->Branch("JpsiKE_e2_bestL1eta" , &JpsiKE_e2_bestL1eta );
  tree_->Branch("JpsiKE_e2_bestL1phi" , &JpsiKE_e2_bestL1phi );
  tree_->Branch("JpsiKE_e2_bestL1Deta" , &JpsiKE_e2_bestL1Deta );
  tree_->Branch("JpsiKE_e2_bestL1Dphi" , &JpsiKE_e2_bestL1Dphi );

  tree_->Branch("JpsiKE_e2_Ele10_match",  &ele2_matchedDiEle[0] );
  tree_->Branch("JpsiKE_e2_Ele9p5_match", &ele2_matchedDiEle[1] );
  tree_->Branch("JpsiKE_e2_Ele9_match",   &ele2_matchedDiEle[2] );
  tree_->Branch("JpsiKE_e2_Ele8p5_match", &ele2_matchedDiEle[3] );
  tree_->Branch("JpsiKE_e2_Ele8_match",   &ele2_matchedDiEle[4] );
  tree_->Branch("JpsiKE_e2_Ele7p5_match", &ele2_matchedDiEle[5] );
  tree_->Branch("JpsiKE_e2_Ele7_match",   &ele2_matchedDiEle[6] );
  tree_->Branch("JpsiKE_e2_Ele6p5_match", &ele2_matchedDiEle[7] );
  tree_->Branch("JpsiKE_e2_Ele6_match",   &ele2_matchedDiEle[8] );
  tree_->Branch("JpsiKE_e2_Ele5p5_match", &ele2_matchedDiEle[9] );
  tree_->Branch("JpsiKE_e2_Ele5_match",   &ele2_matchedDiEle[10] );
  tree_->Branch("JpsiKE_e2_Ele4p5_match", &ele2_matchedDiEle[11] );
  tree_->Branch("JpsiKE_e2_Ele4_match",   &ele2_matchedDiEle[12] );

  tree_->Branch("JpsiKE_e2_Ele10_pt",  &ele2_matchedDiEle_pt[0] );
  tree_->Branch("JpsiKE_e2_Ele9p5_pt", &ele2_matchedDiEle_pt[1] );
  tree_->Branch("JpsiKE_e2_Ele9_pt",   &ele2_matchedDiEle_pt[2] );
  tree_->Branch("JpsiKE_e2_Ele8p5_pt", &ele2_matchedDiEle_pt[3] );
  tree_->Branch("JpsiKE_e2_Ele8_pt",   &ele2_matchedDiEle_pt[4] );
  tree_->Branch("JpsiKE_e2_Ele7p5_pt", &ele2_matchedDiEle_pt[5] );
  tree_->Branch("JpsiKE_e2_Ele7_pt",   &ele2_matchedDiEle_pt[6] );
  tree_->Branch("JpsiKE_e2_Ele6p5_pt", &ele2_matchedDiEle_pt[7] );
  tree_->Branch("JpsiKE_e2_Ele6_pt",   &ele2_matchedDiEle_pt[8] );
  tree_->Branch("JpsiKE_e2_Ele5p5_pt", &ele2_matchedDiEle_pt[9] );
  tree_->Branch("JpsiKE_e2_Ele5_pt",   &ele2_matchedDiEle_pt[10] );
  tree_->Branch("JpsiKE_e2_Ele4p5_pt", &ele2_matchedDiEle_pt[11] );
  tree_->Branch("JpsiKE_e2_Ele4_pt",   &ele2_matchedDiEle_pt[12] );

  tree_->Branch("JpsiKE_e2_Ele10_eta",  &ele2_matchedDiEle_eta[0] );
  tree_->Branch("JpsiKE_e2_Ele9p5_eta", &ele2_matchedDiEle_eta[1] );
  tree_->Branch("JpsiKE_e2_Ele9_eta",   &ele2_matchedDiEle_eta[2] );
  tree_->Branch("JpsiKE_e2_Ele8p5_eta", &ele2_matchedDiEle_eta[3] );
  tree_->Branch("JpsiKE_e2_Ele8_eta",   &ele2_matchedDiEle_eta[4] );
  tree_->Branch("JpsiKE_e2_Ele7p5_eta", &ele2_matchedDiEle_eta[5] );
  tree_->Branch("JpsiKE_e2_Ele7_eta",   &ele2_matchedDiEle_eta[6] );
  tree_->Branch("JpsiKE_e2_Ele6p5_eta", &ele2_matchedDiEle_eta[7] );
  tree_->Branch("JpsiKE_e2_Ele6_eta",   &ele2_matchedDiEle_eta[8] );
  tree_->Branch("JpsiKE_e2_Ele5p5_eta", &ele2_matchedDiEle_eta[9] );
  tree_->Branch("JpsiKE_e2_Ele5_eta",   &ele2_matchedDiEle_eta[10] );
  tree_->Branch("JpsiKE_e2_Ele4p5_eta", &ele2_matchedDiEle_eta[11] );
  tree_->Branch("JpsiKE_e2_Ele4_eta",   &ele2_matchedDiEle_eta[12] );

  tree_->Branch("JpsiKE_e2_Ele10_phi",  &ele2_matchedDiEle_phi[0] );
  tree_->Branch("JpsiKE_e2_Ele9p5_phi", &ele2_matchedDiEle_phi[1] );
  tree_->Branch("JpsiKE_e2_Ele9_phi",   &ele2_matchedDiEle_phi[2] );
  tree_->Branch("JpsiKE_e2_Ele8p5_phi", &ele2_matchedDiEle_phi[3] );
  tree_->Branch("JpsiKE_e2_Ele8_phi",   &ele2_matchedDiEle_phi[4] );
  tree_->Branch("JpsiKE_e2_Ele7p5_phi", &ele2_matchedDiEle_phi[5] );
  tree_->Branch("JpsiKE_e2_Ele7_phi",   &ele2_matchedDiEle_phi[6] );
  tree_->Branch("JpsiKE_e2_Ele6p5_phi", &ele2_matchedDiEle_phi[7] );
  tree_->Branch("JpsiKE_e2_Ele6_phi",   &ele2_matchedDiEle_phi[8] );
  tree_->Branch("JpsiKE_e2_Ele5p5_phi", &ele2_matchedDiEle_phi[9] );
  tree_->Branch("JpsiKE_e2_Ele5_phi",   &ele2_matchedDiEle_phi[10] );
  tree_->Branch("JpsiKE_e2_Ele4p5_phi", &ele2_matchedDiEle_phi[11] );
  tree_->Branch("JpsiKE_e2_Ele4_phi",   &ele2_matchedDiEle_phi[12] );

  tree_->Branch("DoubleEle10_fired", &DoubleEle_fired[0] );
  tree_->Branch("DoubleEle9p5_fired", &DoubleEle_fired[1] );
  tree_->Branch("DoubleEle9_fired", &DoubleEle_fired[2] );
  tree_->Branch("DoubleEle8p5_fired", &DoubleEle_fired[3] );
  tree_->Branch("DoubleEle8_fired", &DoubleEle_fired[4] );
  tree_->Branch("DoubleEle7p5_fired", &DoubleEle_fired[5] );
  tree_->Branch("DoubleEle7_fired", &DoubleEle_fired[6] );
  tree_->Branch("DoubleEle6p5_fired", &DoubleEle_fired[7] );
  tree_->Branch("DoubleEle6_fired", &DoubleEle_fired[8] );
  tree_->Branch("DoubleEle5p5_fired", &DoubleEle_fired[9] );
  tree_->Branch("DoubleEle5_fired", &DoubleEle_fired[10] );
  tree_->Branch("DoubleEle4p5_fired", &DoubleEle_fired[11] );
  tree_->Branch("DoubleEle4_fired", &DoubleEle_fired[12] );

  tree_->Branch("JpsiKE_e1_trgobj_pt", &JpsiKE_e1_trgobj_pt );
  tree_->Branch("JpsiKE_e1_trgobj_eta", &JpsiKE_e1_trgobj_eta );
  tree_->Branch("JpsiKE_e1_trgobj_phi", &JpsiKE_e1_trgobj_phi );
  tree_->Branch("JpsiKE_e1_trgobj_mass", &JpsiKE_e1_trgobj_mass );
  tree_->Branch("JpsiKE_e1_trgobj_q", &JpsiKE_e1_trgobj_q );
  tree_->Branch("JpsiKE_e1_trgobj_vx"   , &JpsiKE_e1_trgobj_vx    );
  tree_->Branch("JpsiKE_e1_trgobj_vy"   , &JpsiKE_e1_trgobj_vy    );
  tree_->Branch("JpsiKE_e1_trgobj_vz"   , &JpsiKE_e1_trgobj_vz    );

  tree_->Branch("JpsiKE_e2_trgobj_pt", &JpsiKE_e2_trgobj_pt );
  tree_->Branch("JpsiKE_e2_trgobj_eta", &JpsiKE_e2_trgobj_eta );
  tree_->Branch("JpsiKE_e2_trgobj_phi", &JpsiKE_e2_trgobj_phi );
  tree_->Branch("JpsiKE_e2_trgobj_mass", &JpsiKE_e2_trgobj_mass );
  tree_->Branch("JpsiKE_e2_trgobj_q", &JpsiKE_e2_trgobj_q );
  tree_->Branch("JpsiKE_e2_trgobj_vx"   , &JpsiKE_e2_trgobj_vx    );
  tree_->Branch("JpsiKE_e2_trgobj_vy"   , &JpsiKE_e2_trgobj_vy    );
  tree_->Branch("JpsiKE_e2_trgobj_vz"   , &JpsiKE_e2_trgobj_vz    );

  tree_->Branch("JpsiKE_Jpsi_pt", &JpsiKE_Jpsi_pt );
  tree_->Branch("JpsiKE_Jpsi_nonfit_pt", &JpsiKE_Jpsi_nonfit_pt );
  tree_->Branch("JpsiKE_Jpsi_eta", &JpsiKE_Jpsi_eta );
  tree_->Branch("JpsiKE_Jpsi_nonfit_eta", &JpsiKE_Jpsi_nonfit_eta );
  tree_->Branch("JpsiKE_Jpsi_phi", &JpsiKE_Jpsi_phi );
  tree_->Branch("JpsiKE_Jpsi_nonfit_phi", &JpsiKE_Jpsi_nonfit_phi);
  tree_->Branch("JpsiKE_Jpsi_mass", &JpsiKE_Jpsi_mass );
  tree_->Branch("JpsiKE_Jpsi_mass_nofit", &JpsiKE_Jpsi_mass_nofit );
  tree_->Branch("JpsiKE_Jpsi_vprob", &JpsiKE_Jpsi_vprob );
  tree_->Branch("JpsiKE_elesDr", &JpsiKE_elesDr );
}

void NanoAnalyzerSingleEle::reset(void){

  run = -1;
  event = -1;
  luminosityBlock = -1;
  nvtx = -1;

  JpsiKE_e1_pt = -99;
  JpsiKE_e1_eta = -99;
  JpsiKE_e1_phi = -99;
  JpsiKE_e1_mass = -99;
  JpsiKE_e1_q = -99;
  JpsiKE_e1_vx = -99;
  JpsiKE_e1_vy = -99;
  JpsiKE_e1_vz = -99;
  JpsiKE_e1_passMVA = -99;
  JpsiKE_e1_bestL1pt = -99;
  JpsiKE_e1_bestL1eta  = -99;
  JpsiKE_e1_bestL1phi  = -99;
  JpsiKE_e1_bestL1Deta = -99;
  JpsiKE_e1_bestL1Dphi = -99;

  JpsiKE_e2_pt = -99;
  JpsiKE_e2_eta = -99;
  JpsiKE_e2_phi = -99;
  JpsiKE_e2_mass = -99;
  JpsiKE_e2_alsotag = -99;
  JpsiKE_e2_q = -99;
  JpsiKE_e2_vx = -99;
  JpsiKE_e2_vy = -99;
  JpsiKE_e2_vz = -99;
  JpsiKE_e2_passMVA = -99;
  JpsiKE_e2_bestL1pt = -99;
  JpsiKE_e2_bestL1eta = -99;
  JpsiKE_e2_bestL1phi = -99;
  JpsiKE_e2_bestL1Deta = -99;
  JpsiKE_e2_bestL1Dphi = -99;

  JpsiKE_e1_trgobj_pt = -99;
  JpsiKE_e1_trgobj_eta = -99;
  JpsiKE_e1_trgobj_phi = -99;
  JpsiKE_e1_trgobj_mass = -99;
  JpsiKE_e1_trgobj_q = -99;
  JpsiKE_e1_trgobj_vx = -99;
  JpsiKE_e1_trgobj_vy = -99;
  JpsiKE_e1_trgobj_vz = -99;

  JpsiKE_e2_trgobj_pt = -99;
  JpsiKE_e2_trgobj_eta = -99;
  JpsiKE_e2_trgobj_phi = -99;
  JpsiKE_e2_trgobj_mass = -99;
  JpsiKE_e2_trgobj_q = -99;
  JpsiKE_e2_trgobj_vx = -99;
  JpsiKE_e2_trgobj_vy = -99;
  JpsiKE_e2_trgobj_vz = -99;

  JpsiKE_Jpsi_pt = -99;
  JpsiKE_Jpsi_nonfit_pt = -99;
  JpsiKE_Jpsi_eta = -99;
  JpsiKE_Jpsi_nonfit_eta = -99;
  JpsiKE_Jpsi_phi = -99;
  JpsiKE_Jpsi_nonfit_phi = -99;
  JpsiKE_Jpsi_mass = -99;
  JpsiKE_Jpsi_mass_nofit = -99;
  JpsiKE_Jpsi_vprob = -99;
  JpsiKE_elesDr = -99;

  for (int ii=0; ii<13; ii++) DoubleEle_fired[ii] = 0.;

  for (int ii=0; ii<13; ii++) dieleobj1_pt[ii]  = -999.;
  for (int ii=0; ii<13; ii++) dieleobj1_eta[ii] = -999.;
  for (int ii=0; ii<13; ii++) dieleobj1_phi[ii] = -999.;
  for (int ii=0; ii<13; ii++) dieleobj2_pt[ii]  = -999.;
  for (int ii=0; ii<13; ii++) dieleobj2_eta[ii] = -999.;
  for (int ii=0; ii<13; ii++) dieleobj2_phi[ii] = -999.;

  for (int ii=0; ii<13; ii++) ele1_matchedDiEle[ii] = 0.;
  for (int ii=0; ii<13; ii++) ele2_matchedDiEle[ii] = 0.;
  for (int ii=0; ii<13; ii++) ele1_matchedDiEle_pt[ii]  = 0.;
  for (int ii=0; ii<13; ii++) ele2_matchedDiEle_pt[ii]  = 0.;
  for (int ii=0; ii<13; ii++) ele1_matchedDiEle_eta[ii] = 0.;
  for (int ii=0; ii<13; ii++) ele2_matchedDiEle_eta[ii] = 0.;
  for (int ii=0; ii<13; ii++) ele1_matchedDiEle_phi[ii] = 0.;
  for (int ii=0; ii<13; ii++) ele2_matchedDiEle_phi[ii] = 0.;
}

DEFINE_FWK_MODULE(NanoAnalyzerSingleEle);
