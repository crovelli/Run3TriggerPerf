#define GCC_VERSION ( 10000 * __GNUC__ + 100 * __GNUC_MINOR__ + __GNUC_PATCHLEVEL__ )
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

// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include "Math/VectorUtil.h"

// Root
#include "TMath.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"

// Utilities
#include "../interface/helperfunc.h"
#include "../interface/MyStruct.h"


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

// FWCore
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

// ------ EXTRA HEADER FILES--------------------//
#include "FWCore/Framework/interface/EventSetup.h"
#ifndef CMSSW12plus
#include "FWCore/Framework/interface/ESHandle.h"
#endif

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"

//**************************
// for trigger information *
//**************************
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

//***************************
// for tracking information *
//***************************
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

//*************************
// for vertex information *
//*************************
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//***********************
// for muon information *
//***********************
#include "DataFormats/PatCandidates/interface/Muon.h"

//*******************************
// for gen particle information *
//*******************************
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// set namespaces
using namespace edm;
using namespace reco;
using namespace std;


//***********************************
// main analyzer class (EDAnalyzer) *
//***********************************

#ifndef CMSSW12plus
// Run 1 and 2
class NanoAnalyzerDoubleMu : public edm::EDAnalyzer
#else
// Run 3
class NanoAnalyzerDoubleMu : public edm::one::EDAnalyzer<>
#endif
{
public:
  explicit NanoAnalyzerDoubleMu(const edm::ParameterSet&);
  ~NanoAnalyzerDoubleMu();

  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

private:

  virtual void beginJob(const edm::ParameterSet& iConfig);
  virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iStp);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void endJob();

  void createBranch();

  void reset();

  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;
  EDGetTokenT<reco::VertexCollection> verticeToken_;
  EDGetTokenT<pat::MuonCollection> muonToken_;
  EDGetTokenT<edm::TriggerResults> triggerToken_;  
  EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobjectToken_;
  edm::EDGetTokenT<BXVector<l1t::Muon>> l1MU_;

  Handle<pat::MuonCollection> muons_;
  Handle< reco::VertexCollection > vertices_;
  Handle< edm::TriggerResults> HLTtriggers_;
  Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;


/////////////////////////////////////////////////////////////////////////
////////////////////////// declare tree, file, //////////////////////////
/////////////////////////////////////////////////////////////////////////
  
  edm::Service<TFileService> fs;
  TTree* tree_;

  /// original nanos
  UInt_t run;
  ULong64_t event;
  UInt_t luminosityBlock;
  int nvtx;

  float  Jpsi_m1_pt;
  float  Jpsi_m1_eta;
  float  Jpsi_m1_phi;
  float  Jpsi_m1_mass;
  int    Jpsi_m1_q;   
  float  Jpsi_m1_looseId;   
  float  Jpsi_m1_bestL1pt;
  float  Jpsi_m1_bestL1eta;
  float  Jpsi_m1_bestL1phi;
  float  Jpsi_m1_bestL1dR;
  //
  float  Jpsi_m2_pt      ;
  float  Jpsi_m2_eta     ;
  float  Jpsi_m2_phi     ;
  float  Jpsi_m2_mass     ;
  int    Jpsi_m2_alsotag ;
  int    Jpsi_m2_q   ;   
  float  Jpsi_m2_looseId  ;
  float  Jpsi_m2_bestL1pt ;
  float  Jpsi_m2_bestL1eta ;
  float  Jpsi_m2_bestL1phi ;
  float  Jpsi_m2_bestL1dR ;
  //
  float  Jpsi_m1_trgobj_pt      ;
  float  Jpsi_m1_trgobj_eta     ;
  float  Jpsi_m1_trgobj_phi     ;
  float  Jpsi_m1_trgobj_mass     ;
  int    Jpsi_m1_trgobj_q      ;   
  float  Jpsi_m1_trgobj_dR      ;
  // 
  float  Jpsi_m2_trgobj_pt ;
  float  Jpsi_m2_trgobj_eta ;
  float  Jpsi_m2_trgobj_phi ;
  float  Jpsi_m2_trgobj_mass ;
  int    Jpsi_m2_trgobj_q ;   
  float  Jpsi_m2_trgobj_dR ;
  //
  int DoubleMu_fired;
  float mu1_matchedDiMu_dR  = 0;
  float mu2_matchedDiMu_dR  = 0;
  float mu1_matchedDiMu_pt  = 0;
  float mu2_matchedDiMu_pt  = 0;
  float mu1_matchedDiMu_eta = 0;
  float mu2_matchedDiMu_eta = 0;
  float mu1_matchedDiMu_phi = 0;
  float mu2_matchedDiMu_phi = 0;
  //
  int L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4;
  int L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6;
  int L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6;
  int L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5;
  int L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4;
  int L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4;
  int L1_DoubleMu4p5_SQ_OS_dR_Max1p2;
  int L1_DoubleMu4_SQ_OS_dR_Max1p2;
  int L1_OR;
  // 
  float  Jpsi_fit_pt;
  float  Jpsi_nonfit_pt;
  float  Jpsi_fit_eta;
  float  Jpsi_nonfit_eta;
  float  Jpsi_fit_phi;
  float  Jpsi_nonfit_phi;
  float  Jpsi_fit_mass;
  float  Jpsi_nonfit_mass;
  float  Jpsi_fit_vprob;
  float  Jpsi_muonsDr;

  // Not for tree, other variables declaration
  helperfunc aux;
  float chi = 0.;
  float ndf = 0.;

  float dimuobj1_pt  = -999.;
  float dimuobj1_eta = -999.;
  float dimuobj1_phi = -999.;
  float dimuobj2_pt  = -999.;
  float dimuobj2_eta = -999.;
  float dimuobj2_phi = -999.;
  TH1F * hist; 

}; // end of class member


NanoAnalyzerDoubleMu::NanoAnalyzerDoubleMu(const edm::ParameterSet& iConfig): 
  bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>())
{
  muonToken_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));

  verticeToken_       = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));

  triggerToken_	      = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLT"));
  triggerobjectToken_ = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerobjects"));

  l1MU_ = consumes<BXVector<l1t::Muon>>(iConfig.getParameter<edm::InputTag>("l1MU"));

  hist = fs->make<TH1F>("cutflow", "cutflow", 10,0,10);
  tree_ = fs->make<TTree>( "tree", "tree" );

  createBranch();

} // end of constructor

NanoAnalyzerDoubleMu::~NanoAnalyzerDoubleMu() { }

void NanoAnalyzerDoubleMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  reset();

  const auto& bField = iSetup.getData(bFieldToken_);

  iEvent.getByToken(verticeToken_, vertices_ );
  iEvent.getByToken(triggerToken_, HLTtriggers_);
  iEvent.getByToken(triggerobjectToken_ , triggerObjects); 

  nvtx = vertices_->size();

  run = (iEvent.id()).run();
  event = (iEvent.id()).event();
  luminosityBlock = (iEvent.id()).luminosityBlock();

  // check if the reference trigger and/or the analysis trigger fired
  bool isTriggered = false;
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*HLTtriggers_);
  std::string refTriggerName="";    
  std::string finalTriggerName="";  

  for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {

    ///////std::cout << "Loop over all triggers: i = " << i << ", name = " << trigNames.triggerName(i) << ", fired = " << HLTtriggers_->accept(i) << std::endl;
    
    // Check if the reference trigger path fired
    if(trigNames.triggerName(i).find("HLT_Mu8_v")!= std::string::npos){             // chiara
      // if(trigNames.triggerName(i).find("HLT_Mu20_v")!= std::string::npos){ 
      if(HLTtriggers_->accept(i)){
	isTriggered = true;
	refTriggerName=trigNames.triggerName(i);  
      }
    }
    
    // Check if the analysis trigger path firef
    if(trigNames.triggerName(i).find("HLT_DoubleMu4_3_LowMass_v")!= std::string::npos){ 
      if(HLTtriggers_->accept(i)){
	DoubleMu_fired = true;
	finalTriggerName=trigNames.triggerName(i); 
      }
    }
  }

  /*
  std::cout << std::endl;
  std::cout << "Ref trigger? " << isTriggered << std::endl;
  std::cout << "AN trigger? "  << DoubleMu_fired << std::endl;
  std::cout << "finalTriggerName = " << finalTriggerName << std::endl;
  std::cout << std::endl;
  */

  // Our reference path fired
  if(!isTriggered) return;    
  hist->Fill(2);


  // Take the muons collection
  iEvent.getByToken(muonToken_ , muons_ );

  // Output collections
  std::vector<pat::Muon> muoncollection;                          // collection with offline muons passing minimal selection
  std::vector<pat::TriggerObjectStandAlone> trg_obj_collection;   // collection with HLT candidate with best match with offline muons
                                                                  // (passing matching criteria) - at most one per muon
  std::vector<int> muonmatched;                                   // integer with the position in the trg_obj_collection of the HLT object matched to muon
  // 
  muoncollection.clear();
  trg_obj_collection.clear();
  muonmatched.clear();

  // std::cout << "event = " << event << ", muons_->size() = " << muons_->size() << std::endl;


  // Offline muons
  for(size_t imuon = 0; imuon < muons_->size(); imuon++){
    const pat::Muon & muon = (*muons_)[imuon];

    // Offline cuts
    if(muon.pt() < 2.5) continue;
    if(fabs(muon.eta()) > 2.6) continue;
    if(!(muon.track().isNonnull())) continue;
    
    // std::cout << "event = " << event << ", imuon = " << imuon << ", muon.pt() =  " << muon.pt() << std::endl;

    // Trigger matching
    bool trigObjMatchMu = false;                      // this offline muon matches a HLT mu-candidate
    pat::TriggerObjectStandAlone best_match_obj;      // this offline muon matches a HLT mu-candidate and this is the best matched candidate
    Float_t best_match_dR = 9999999;                  // this offline muon matches a HLT mu-candidate and this is the best match DR

    // Loop over trigger objects matching the reference path
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
      
      // consider only objects which match the ref path    
      obj.unpackPathNames(trigNames);
      obj.unpackFilterLabels(iEvent, *HLTtriggers_);
      std::vector<std::string> pathNamesAll  = obj.pathNames(false);
      bool isPathExist = false;
      for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
	if(pathNamesAll[h]==refTriggerName) isPathExist = true;     
      }
      if(!isPathExist) continue;
      
      int muObjNumber = -1;
      for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){	

	//std::cout << "Event: Filter " << hh << " => " << obj.filterLabels()[hh] << " ";
	//std::cout << "" << std::endl;
	  

	if(obj.filterLabels()[hh].find("hltL3fL1sMu5L1f0L2f5L3Filtered8") != std::string::npos) {  // chiara  HLT_Mu8
	  // if(obj.filterLabels()[hh].find("hltL3fL1sMu18L1f0L2f10QL3Filtered20Q") != std::string::npos) {  // chiara HLT_Mu20
	  muObjNumber = hh;
	}
      }
        
      // here HLT obj vs reco muon candidates
      TVector3 muTV3, objTV3;
      muTV3.SetPtEtaPhi( muon.pt(), muon.eta(), muon.phi() );
      objTV3.SetPtEtaPhi( obj.pt(), obj.eta(), obj.phi() );
      Float_t deltaR = fabs(muTV3.DeltaR(objTV3));
    
      // here HLT-muon candidates
      if (muObjNumber>=0) {
	//std::cout<< "DeltaR = " << deltaR << std::endl;
	//std::cout << "This is a muon HLT candidate" << endl;
	if(deltaR < 0.1){    // chiara
	  trigObjMatchMu = true;
	  if (deltaR < best_match_dR){
	    best_match_dR = deltaR;
	    best_match_obj = obj;
	  }
	  //std::cout << "This object is matched with muon: deltaPhi = " << best_match_dR = " << best_match_dR << std::endl;
	  //std::cout << "Offline: " << muon.pt() << " " << muon.eta() << " " << muon.phi() << std::endl;
	  //std::cout << "HLT: " << obj.pt() << " " << obj.eta() << " " << obj.phi() << std::endl;
	}
      }

      //std::cout << "event = " << event << " : muObjNumber = " << muObjNumber << std::endl;
      //std::cout << "trigObjMatchMu = " << trigObjMatchMu << std::endl;
      //std::cout << "In objects loop, best_match_obj => " << best_match_obj.pt() << " " << best_match_obj.eta() << " " << best_match_obj.phi() << std::endl;
      
    } // Loop over trigger object
    
    // std::cout << "After objects loop, best_match_obj => " << best_match_obj.pt() << " " << best_match_obj.eta() << " " << best_match_obj.phi() << std::endl;

    
    // This is to further process 
    muoncollection.push_back(muon);
    if(trigObjMatchMu){ 
      trg_obj_collection.push_back(best_match_obj);
      muonmatched.push_back(int(trg_obj_collection.size())-1);
    }
    else{
      muonmatched.push_back(-999);	
    }

  } // Loop over offline muons

  // std::cout << "event = " << event << ", muoncollection.size() = " << muoncollection.size() << ", muonmatched.size() = " << muonmatched.size() << std::endl;
  
  if (muoncollection.size() < 2) return; 
  hist->Fill(3);

  // NB: muoncollection and muonmatched must have the same size
  if (muoncollection.size() != muonmatched.size() ) return;
  hist->Fill(4);


  // Prepare offline muon pairs
  float jpsi_max_pt = -1;
  int mcidx_m1 = -1;
  int mcidx_m2 = -1;
  int mcidx_trgobj1 = -1;
  int mcidx_trgobj2 = -1;
  TLorentzVector jpsi_tlv_highest;

  for(int im = 0; im < (int)muoncollection.size(); im++){
    for(int jm = im+1; jm < (int)muoncollection.size(); jm++){
      
      // at least 1 muon of the jpsi pair should match a trigger muon
      if ( muonmatched[im]<0 && muonmatched[jm]<0 ) continue;
      
      // std::cout << "im = " << im << ", muonmatched[im] = " << muonmatched[im] << ", jm = " << jm << ", muonmatched[jm] = " << muonmatched[jm] << std::endl; 

      // match with trigger muon
      int this_mcidx_trgobj1 = -1;
      int this_mcidx_trgobj2 = -1;
      if (muonmatched[im]>-999){
	this_mcidx_trgobj1 = muonmatched[im];
	// std::cout << "mcidx_trgobj1 = " << this_mcidx_trgobj1 << std::endl; 
      }
      if (muonmatched[jm]>-999){
	this_mcidx_trgobj2 = muonmatched[jm];
	// std::cout << "mcidx_trgobj2 = " << this_mcidx_trgobj2 << std::endl; 
      }
      
      const pat::Muon mu1 = muoncollection[im];
      const pat::Muon mu2 = muoncollection[jm];


      /*
      std::cout << "Mu1 " << im << ", pt = " << mu1.pt() << ", eta = " << mu1.eta() << ", phi = " << mu1.phi() << ", id = " << mu1.pdgId() << std::endl;
      
      if (this_mcidx_trgobj1>=0) std::cout << "TrgObj1 pt = " << trg_obj_collection[this_mcidx_trgobj1].pt() << ", eta = " << trg_obj_collection[this_mcidx_trgobj1].eta() << ", phi = " << trg_obj_collection[this_mcidx_trgobj1].phi() << std::endl; 
      else 
      std::cout << "TrgObj1 not found" << std::endl;
      
      std::cout << "Mu2 " << jm << ", pt = " << mu2.pt() << ", eta = " << mu2.eta() << ", phi = " << mu2.phi() << ", id = " << mu2.pdgId() << std::endl;
      
      if (this_mcidx_trgobj2>=0) std::cout << "TrgObj2 pt = " << trg_obj_collection[this_mcidx_trgobj2].pt() << ", eta = " << trg_obj_collection[this_mcidx_trgobj2].eta() << ", phi = " << trg_obj_collection[this_mcidx_trgobj2].phi() << std::endl; 
      else 
      std::cout << "TrgObj2 not found" << std::endl;
      */


      TLorentzVector tlv_m1;
      TLorentzVector tlv_m2;
      tlv_m1.SetPtEtaPhiM(mu1.pt(), mu1.eta(), mu1.phi(), aux.mass_muon);
      tlv_m2.SetPtEtaPhiM(mu2.pt(), mu2.eta(), mu2.phi(), aux.mass_muon);

      TLorentzVector tlv_jpsi = (tlv_m1 + tlv_m2);
      float jpsi_mass = tlv_jpsi.M();
      float jpsi_pt = tlv_jpsi.Pt();

      // std::cout << "event = " << event << ", jpsi_mass = " << jpsi_mass << std::endl;

      if (mu1.charge() + mu2.charge() !=0) continue;
      if (jpsi_mass < 2.0) continue; 
      if (jpsi_mass > 4.0) continue;
      // std::cout << "event = " << event << ": jpsi_mass = " << jpsi_mass << ", jpsi_max_pt = " << jpsi_max_pt << ", jpsi_pt = " << jpsi_pt << std::endl;
      
      if(jpsi_max_pt < jpsi_pt){
	jpsi_max_pt = jpsi_pt;
	mcidx_m1 = im;
	mcidx_m2 = jm;
	mcidx_trgobj1 = this_mcidx_trgobj1;
	mcidx_trgobj2 = this_mcidx_trgobj2;
	jpsi_tlv_highest = tlv_jpsi;
      }
      // std::cout << "jpsi_mass = " << jpsi_mass << ", jpsi_max_pt = " << jpsi_max_pt << ", jpsi_pt = " << jpsi_pt << std::endl;      
      // std::cout << "event = " << event << "; In the loop: mcidx_m1 = " << mcidx_m1 << ", mcidx_m2 = " << mcidx_m2 << std::endl;
    }
  }

  // std::cout << "event = " << event << ", jpsi_max_pt = " << jpsi_max_pt << std::endl;

  // At least 1 reco J/psi
  if(jpsi_max_pt == -1) return;  
  hist->Fill(5);

  // std::cout << "event = " << event << ", jpsi_max_pt = " << jpsi_max_pt << " => passed " << std::endl;


  // Offline muons
  float reco_m1_pt  = muoncollection[mcidx_m1].pt(); 
  float reco_m1_eta = muoncollection[mcidx_m1].eta(); 
  float reco_m1_phi = muoncollection[mcidx_m1].phi(); 
  float reco_m2_pt  = muoncollection[mcidx_m2].pt(); 
  float reco_m2_eta = muoncollection[mcidx_m2].eta(); 
  float reco_m2_phi = muoncollection[mcidx_m2].phi(); 
  TVector3 mu1TV3, mu2TV3;
  mu1TV3.SetPtEtaPhi( reco_m1_pt, reco_m1_eta, reco_m1_phi );
  mu2TV3.SetPtEtaPhi( reco_m2_pt, reco_m2_eta, reco_m2_phi );
  // std::cout << "Offline1: " << reco_m1_pt << " " << reco_m1_eta << " " << reco_m1_phi << std::endl;
  // std::cout << "Offline2: " << reco_m2_pt << " " << reco_m2_eta << " " << reco_m2_phi << std::endl;


  // -------------------------------------------------------
  // Loop over di-mu path, when fired
  if (DoubleMu_fired==1) {   

    std::string DoubleMuTrigName = "HLT_DoubleMu4_3_LowMass_v";       // chiara
    std::string DoubleMuObjName  = "hltDoubleMu43LowMassL3Filtered";  // chiara
    
    // Loop over trigger objects matching the dimuon path 
    int theHLTobj=-1;
    vector<int> idxHLTobj;
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
      theHLTobj++;

      // check if this obj comes from the wanted di-muon path
      obj.unpackPathNames(trigNames);
      obj.unpackFilterLabels(iEvent, *HLTtriggers_);
      std::vector<std::string> pathNamesAll = obj.pathNames(false);
      bool isPathExist = false;
      for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
	// std::cout << "event " << event << " => path " << h << " => " << pathNamesAll[h] << std::endl;
	if(pathNamesAll[h].find(DoubleMuTrigName)!= std::string::npos) {
	  isPathExist = true;   
	  // std::cout << "event " << event << " => found" << std::endl;
	}
      }
      if(!isPathExist) continue;
	
      // vector with objects coming from the correct HLT filter
      bool okMatched=false;
      for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){	
	if(obj.filterLabels()[hh].find(DoubleMuObjName) != std::string::npos) { okMatched=true; }
      }
      if (!okMatched) continue;
      
      idxHLTobj.push_back(theHLTobj);
    }	

    // Build online object pairs which satisfy the trigger requirements
    int theHLTobj1=-1;
    for (pat::TriggerObjectStandAlone obj1 : *triggerObjects) {
      theHLTobj1++;

      int theHLTobj2=-1;
      for (pat::TriggerObjectStandAlone obj2 : *triggerObjects) {
	theHLTobj2++;

	// can not be the same
	if (theHLTobj1==theHLTobj2) continue;
	
	// both selected above
	bool found1 = std::find(idxHLTobj.begin(), idxHLTobj.end(), theHLTobj1) != idxHLTobj.end();
	bool found2 = std::find(idxHLTobj.begin(), idxHLTobj.end(), theHLTobj2) != idxHLTobj.end();
	if (!found1 || !found2) continue;
	  	  
	TLorentzVector tlv_obj1;
	TLorentzVector tlv_obj2;
	tlv_obj1.SetPtEtaPhiM(obj1.pt(), obj1.eta(), obj1.phi(), aux.mass_muon);
	tlv_obj2.SetPtEtaPhiM(obj2.pt(), obj2.eta(), obj2.phi(), aux.mass_muon);

	TVector3 tv3_obj1;
	TVector3 tv3_obj2;
	tv3_obj1.SetPtEtaPhi(obj1.pt(), obj1.eta(), obj1.phi());
	tv3_obj2.SetPtEtaPhi(obj2.pt(), obj2.eta(), obj2.phi());

	float objmass = (tlv_obj1+tlv_obj2).M();
	float objpt   = (tlv_obj1+tlv_obj2).Pt();

	if (objpt<4.9) continue;
	if (objmass<0.2 || objmass>8.5) continue;

	// match with offline muons
	float deltaR11 = fabs(mu1TV3.DeltaR(tv3_obj1));
	float deltaR21 = fabs(mu2TV3.DeltaR(tv3_obj1));
	float deltaR12 = fabs(mu1TV3.DeltaR(tv3_obj2));
	float deltaR22 = fabs(mu2TV3.DeltaR(tv3_obj2));
	
	if (deltaR11<0.1 && deltaR11<deltaR12 && deltaR11<mu1_matchedDiMu_dR && deltaR22<0.1 && deltaR22<deltaR21 && deltaR22<mu2_matchedDiMu_dR) { 
	  mu1_matchedDiMu_dR  = deltaR11;
	  mu2_matchedDiMu_dR  = deltaR22;
	  mu1_matchedDiMu_pt  = obj1.pt();
	  mu2_matchedDiMu_pt  = obj2.pt();
	  mu1_matchedDiMu_eta = obj1.eta();
	  mu2_matchedDiMu_eta = obj2.eta();
	  mu1_matchedDiMu_phi = obj1.phi();
	  mu2_matchedDiMu_phi = obj2.phi();
	}
	if (deltaR12<0.1 && deltaR12<deltaR11 && deltaR12<mu1_matchedDiMu_dR && deltaR21<0.1 && deltaR21<deltaR22 && deltaR21<mu2_matchedDiMu_dR) {
	  mu1_matchedDiMu_dR  = deltaR12;
	  mu2_matchedDiMu_dR  = deltaR21;
	  mu1_matchedDiMu_pt  = obj2.pt();
	  mu2_matchedDiMu_pt  = obj1.pt();
	  mu1_matchedDiMu_eta = obj2.eta();
	  mu2_matchedDiMu_eta = obj1.eta();
	  mu1_matchedDiMu_phi = obj2.phi();
	  mu2_matchedDiMu_phi = obj1.phi();
	}
	
      } // Loop over trigger objects
    }   // Loop over trigger objects   
    
  } // path fired


  // -----------------------------------------------
  // Match L1 / offline
  float bestMatchM1_eta = -99.;
  float bestMatchM1_phi = -99.;
  float bestMatchM1_pt  = -99.;
  float bestMatchM1_dR  = 999.;
  float bestMatchM2_eta = -99.;
  float bestMatchM2_phi = -99.;
  float bestMatchM2_pt  = -99.;
  float bestMatchM2_dR  = 999.;


  edm::Handle<BXVector<l1t::Muon> > gmuons;   
  iEvent.getByToken(l1MU_, gmuons);   

  for (auto it1 = gmuons->begin(0); it1 != gmuons->end(0); ++it1) {
    for (auto it2 = gmuons->begin(0); it2 != gmuons->end(0); ++it2) {

      // can not be the same
      if (it1==it2) continue;

      // L1 candidate quality
      int obj1qual = it1->hwQual();
      int obj2qual = it2->hwQual();
      if (obj1qual<12) continue;
      if (obj2qual<12) continue;

      // infos propagated at vtx (default is at mu chambers)
      TLorentzVector tlv_obj1;
      TLorentzVector tlv_obj2;
      tlv_obj1.SetPtEtaPhiM(it1->pt(), it1->etaAtVtx(), it1->phiAtVtx(), aux.mass_muon);
      tlv_obj2.SetPtEtaPhiM(it2->pt(), it2->etaAtVtx(), it2->phiAtVtx(), aux.mass_muon);

      TVector3 tv3_obj1;
      TVector3 tv3_obj2;
      tv3_obj1.SetPtEtaPhi(it1->pt(), it1->etaAtVtx(), it1->phiAtVtx());
      tv3_obj2.SetPtEtaPhi(it2->pt(), it2->etaAtVtx(), it2->phiAtVtx());
      
      float q_obj1 = it1->charge();
      float q_obj2 = it2->charge();

      float deltaR   = fabs(tv3_obj1.DeltaR(tv3_obj2));
      float deltaEta = fabs(tv3_obj1.Eta()-(tv3_obj2.Eta()));

      // Bool for each L1 seed considered in the OR 
      // L1SeedsLogicalExpression = cms.string( 
      // "L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4 OR 
      // L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6 OR 
      // L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6 OR 
      // L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5 OR
      // L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 OR 
      // L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 OR 
      // L1_DoubleMu4p5_SQ_OS_dR_Max1p2 OR L1_DoubleMu4_SQ_OS_dR_Max1p2" ),

      bool theL1_DoubleMu3er2p0_SQ_OS_dR_Max1p4   = (it1->pt()>3) && (it2->pt()>3) && (fabs(it1->etaAtVtx())<2.0) && (fabs(it2->etaAtVtx())<2.0) && (q_obj1*q_obj2<0) && deltaR<1.4;
      bool theL1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6 = (fabs(it1->etaAtVtx())<2.0) && (fabs(it2->etaAtVtx())<2.0) && (q_obj1*q_obj2<0) && deltaEta<1.6;
      bool theL1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6 = (fabs(it1->etaAtVtx())<1.4) && (fabs(it2->etaAtVtx())<1.4) && (q_obj1*q_obj2<0) && deltaEta<1.6;
      bool theL1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5 = (fabs(it1->etaAtVtx())<2.0) && (fabs(it2->etaAtVtx())<2.0) && (q_obj1*q_obj2<0) && deltaR<1.5;
      bool theL1_DoubleMu0er1p4_SQ_OS_dR_Max1p4   = (fabs(it1->etaAtVtx())<1.4) && (fabs(it2->etaAtVtx())<1.4) && (q_obj1*q_obj2<0) && deltaR<1.4;
      bool theL1_DoubleMu0er1p5_SQ_OS_dR_Max1p4   = (fabs(it1->etaAtVtx())<1.5) && (fabs(it2->etaAtVtx())<1.5) && (q_obj1*q_obj2<0) && deltaR<1.4;
      bool theL1_DoubleMu4p5_SQ_OS_dR_Max1p2      = (it1->pt()>4.5) && (it2->pt()>4.5) && (q_obj1*q_obj2<0) && deltaR<1.2;
      bool theL1_DoubleMu4_SQ_OS_dR_Max1p2        = (it1->pt()>4.0) && (it2->pt()>4.0) && (q_obj1*q_obj2<0) && deltaR<1.2;
      bool theL1_OR = theL1_DoubleMu3er2p0_SQ_OS_dR_Max1p4 || theL1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6 || theL1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6 || theL1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5 || theL1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 || theL1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 || theL1_DoubleMu4p5_SQ_OS_dR_Max1p2 || theL1_DoubleMu4_SQ_OS_dR_Max1p2;
      
      // Pass one of the seeds
      if (theL1_OR!=1) continue;

      // match with offline muons
      float deltaR11 = fabs(mu1TV3.DeltaR(tv3_obj1));
      float deltaR21 = fabs(mu2TV3.DeltaR(tv3_obj1));
      float deltaR12 = fabs(mu1TV3.DeltaR(tv3_obj2));
      float deltaR22 = fabs(mu2TV3.DeltaR(tv3_obj2));

      if (deltaR11<0.3 && deltaR11<deltaR12 && deltaR11<bestMatchM1_dR && deltaR22<0.3 && deltaR22<deltaR21 && deltaR22<bestMatchM2_dR) { 
	bestMatchM1_dR  = deltaR11;
	bestMatchM2_dR  = deltaR22;
	bestMatchM1_pt  = it1->pt();
	bestMatchM2_pt  = it2->pt();
	bestMatchM1_eta = it1->etaAtVtx();
        bestMatchM2_eta = it2->etaAtVtx();
	bestMatchM1_phi = it1->phiAtVtx();
	bestMatchM2_phi = it2->phiAtVtx();
	L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4   = theL1_DoubleMu3er2p0_SQ_OS_dR_Max1p4;
	L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6 = theL1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6;
	L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6 = theL1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6;
	L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5 = theL1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5;
	L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4   = theL1_DoubleMu0er1p4_SQ_OS_dR_Max1p4;
	L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4   = theL1_DoubleMu0er1p5_SQ_OS_dR_Max1p4;
	L1_DoubleMu4p5_SQ_OS_dR_Max1p2      = theL1_DoubleMu4p5_SQ_OS_dR_Max1p2;
	L1_DoubleMu4_SQ_OS_dR_Max1p2        = theL1_DoubleMu4_SQ_OS_dR_Max1p2;
	L1_OR = theL1_OR;
      }
      if (deltaR12<0.3 && deltaR12<deltaR11 && deltaR12<bestMatchM1_dR && deltaR21<0.3 && deltaR21<deltaR22 && deltaR21<bestMatchM2_dR) {
	bestMatchM1_dR  = deltaR12;
        bestMatchM2_dR  = deltaR21;
	bestMatchM1_pt  = it2->pt();
        bestMatchM2_pt  = it1->pt();
	bestMatchM1_eta = it2->etaAtVtx();
	bestMatchM2_eta = it1->etaAtVtx();
	bestMatchM1_phi = it2->phiAtVtx();
	bestMatchM2_phi = it1->phiAtVtx();
	L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4   = theL1_DoubleMu3er2p0_SQ_OS_dR_Max1p4;
	L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6 = theL1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6;
	L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6 = theL1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6;
	L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5 = theL1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5;
	L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4   = theL1_DoubleMu0er1p4_SQ_OS_dR_Max1p4;
	L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4   = theL1_DoubleMu0er1p5_SQ_OS_dR_Max1p4;
	L1_DoubleMu4p5_SQ_OS_dR_Max1p2      = theL1_DoubleMu4p5_SQ_OS_dR_Max1p2;
	L1_DoubleMu4_SQ_OS_dR_Max1p2        = theL1_DoubleMu4_SQ_OS_dR_Max1p2;
	L1_OR = theL1_OR;
      }
    } // Loop over L1 candidates
  }   // Loop over L1 candidates


  // ---------------------------------------------------------------
  // Swap 1<->2 so that: 1 is always tag, 2 is always probe, and additional variable tells if 2 is also tag
  int mcidx_tag = -99;
  int mcidx_pro = -99;
  int mcidx_trgobj_tag = -99;
  int mcidx_trgobj_pro = -99;
  int is_probe_also_tag = -99;

  if (mcidx_trgobj1>=0) {
    mcidx_pro = mcidx_m2;
    mcidx_trgobj_pro = mcidx_trgobj2;
    mcidx_tag = mcidx_m1;
    mcidx_trgobj_tag = mcidx_trgobj1;
    if (mcidx_trgobj2<0)  is_probe_also_tag = 0;
    if (mcidx_trgobj2>=0) is_probe_also_tag = 1;
  } else {  
    mcidx_pro = mcidx_m1;
    mcidx_trgobj_pro = mcidx_trgobj1;
    mcidx_tag = mcidx_m2;
    mcidx_trgobj_tag = mcidx_trgobj2;
    is_probe_also_tag = 0;
  }

  mcidx_m1 = mcidx_tag;
  mcidx_m2 = mcidx_pro;
  mcidx_trgobj1 = mcidx_trgobj_tag;
  mcidx_trgobj2 = mcidx_trgobj_pro;



  // ------------------------------------------------------
  // Kin fit
  const reco::TransientTrack muon1TT((*(muoncollection[mcidx_m1].bestTrack())),&bField);  
  const reco::TransientTrack muon2TT((*(muoncollection[mcidx_m2].bestTrack())),&bField);
  KinematicParticleFactoryFromTransientTrack pFactory;
  std::vector<RefCountedKinematicParticle> muonParticles;
  muonParticles.push_back(pFactory.particle(muon1TT, aux.muon_mass, chi, ndf, aux.muon_sigma));
  muonParticles.push_back(pFactory.particle(muon2TT, aux.muon_mass, chi, ndf, aux.muon_sigma));

  // Kinematic fit of the two muons to a common vtx
  RefCountedKinematicParticle jpsi_part;
  RefCountedKinematicVertex jpsi_vertex;
  RefCountedKinematicTree jpTree;
  Bool_t jpsifit_flag;
  std::tie(jpsifit_flag, jpsi_part, jpsi_vertex, jpTree) = aux.KinematicFit(muonParticles, -1, -1);
  
  // Successfull kin fit
  if(!jpsifit_flag) return;
  hist->Fill(6);

  Jpsi_fit_pt    = jpsi_part->currentState().globalMomentum().perp();
  Jpsi_fit_eta   = jpsi_part->currentState().globalMomentum().eta();
  Jpsi_fit_phi   = jpsi_part->currentState().globalMomentum().phi();
  Jpsi_fit_mass  = jpsi_part->currentState().mass();
  Jpsi_fit_vprob = TMath::Prob(jpsi_part->chiSquared(), jpsi_part->degreesOfFreedom());



  // Infos about JPsi candidate
  Jpsi_m1_pt   = muoncollection[mcidx_m1].pt();
  Jpsi_m1_eta  = muoncollection[mcidx_m1].eta();
  Jpsi_m1_phi  = muoncollection[mcidx_m1].phi();
  Jpsi_m1_mass = muoncollection[mcidx_m1].mass();
  Jpsi_m1_q    = muoncollection[mcidx_m1].charge();
  Jpsi_m1_looseId = muoncollection[mcidx_m1].isLooseMuon();
  Jpsi_m1_bestL1pt   = bestMatchM1_pt;
  Jpsi_m1_bestL1eta  = bestMatchM1_eta;
  Jpsi_m1_bestL1phi  = bestMatchM1_phi;
  Jpsi_m1_bestL1dR   = bestMatchM1_dR;  

  Jpsi_m2_pt   = muoncollection[mcidx_m2].pt();
  Jpsi_m2_eta  = muoncollection[mcidx_m2].eta();
  Jpsi_m2_phi  = muoncollection[mcidx_m2].phi();
  Jpsi_m2_mass = muoncollection[mcidx_m2].mass();
  Jpsi_m2_alsotag = is_probe_also_tag;
  Jpsi_m2_q    = muoncollection[mcidx_m2].charge();
  Jpsi_m2_looseId = muoncollection[mcidx_m2].isLooseMuon();
  Jpsi_m2_bestL1pt  = bestMatchM2_pt;
  Jpsi_m2_bestL1eta = bestMatchM2_eta;
  Jpsi_m2_bestL1phi = bestMatchM2_phi;
  Jpsi_m2_bestL1dR  = bestMatchM2_dR;  
  
  if (mcidx_trgobj1>=0) {
    
    Jpsi_m1_trgobj_pt   = trg_obj_collection[mcidx_trgobj1].pt();
    Jpsi_m1_trgobj_eta  = trg_obj_collection[mcidx_trgobj1].eta();
    Jpsi_m1_trgobj_phi  = trg_obj_collection[mcidx_trgobj1].phi();
    Jpsi_m1_trgobj_mass = trg_obj_collection[mcidx_trgobj1].mass();
    Jpsi_m1_trgobj_q    = trg_obj_collection[mcidx_trgobj1].charge();
    TVector3 trgobj1TV3;
    trgobj1TV3.SetPtEtaPhi( Jpsi_m1_trgobj_pt, Jpsi_m1_trgobj_eta, Jpsi_m1_trgobj_phi);
    Jpsi_m1_trgobj_dR = fabs(mu1TV3.DeltaR(trgobj1TV3));
  }

  if (mcidx_trgobj2>=0) {

    Jpsi_m2_trgobj_pt   = trg_obj_collection[mcidx_trgobj2].pt();
    Jpsi_m2_trgobj_eta  = trg_obj_collection[mcidx_trgobj2].eta();
    Jpsi_m2_trgobj_phi  = trg_obj_collection[mcidx_trgobj2].phi();
    Jpsi_m2_trgobj_mass = trg_obj_collection[mcidx_trgobj2].mass();
    Jpsi_m2_trgobj_q    = trg_obj_collection[mcidx_trgobj2].charge();
    TVector3 trgobj2TV3;
    trgobj2TV3.SetPtEtaPhi( Jpsi_m2_trgobj_pt, Jpsi_m2_trgobj_eta, Jpsi_m2_trgobj_phi);
    Jpsi_m2_trgobj_dR = fabs(mu2TV3.DeltaR(trgobj2TV3));
  }

  Jpsi_nonfit_pt   = jpsi_tlv_highest.Pt();
  Jpsi_nonfit_eta  = jpsi_tlv_highest.Eta();
  Jpsi_nonfit_phi  = jpsi_tlv_highest.Phi();
  Jpsi_nonfit_mass = jpsi_tlv_highest.M(); 

  Jpsi_muonsDr = mu1TV3.DeltaR(mu2TV3); 

  tree_->Fill();

  return;

} //NanoAnalyzer::analyze ends



//**************************************************
//************* additional methods *****************
//**************************************************

void NanoAnalyzerDoubleMu::beginJob(const edm::ParameterSet& iConfig) { }

void NanoAnalyzerDoubleMu::beginRun(const edm::Run &iRun, const edm::EventSetup &iStp) { }

void NanoAnalyzerDoubleMu::fillDescriptions(edm::ConfigurationDescriptions & descriptions) { }

void NanoAnalyzerDoubleMu::endRun(edm::Run const&, edm::EventSetup const&) { }

void NanoAnalyzerDoubleMu::endJob() { }

//define this as a plug-in

// branch title creation
void NanoAnalyzerDoubleMu::createBranch() { 

  tree_->Branch("run", &run, "run/i");
  tree_->Branch("event", &event, "event/l");
  tree_->Branch("luminosityBlock", &luminosityBlock, "luminosityBlock/i");
  tree_->Branch("nvtx", &nvtx, "nvtx/i");

  // Offline
  tree_->Branch("Jpsi_m1_pt", &Jpsi_m1_pt );
  tree_->Branch("Jpsi_m1_eta", &Jpsi_m1_eta );
  tree_->Branch("Jpsi_m1_phi", &Jpsi_m1_phi );
  tree_->Branch("Jpsi_m1_mass", &Jpsi_m1_mass );
  tree_->Branch("Jpsi_m1_q", &Jpsi_m1_q );
  tree_->Branch("Jpsi_m1_looseId"   , &Jpsi_m1_looseId    );
  tree_->Branch("Jpsi_m1_bestL1pt",   &Jpsi_m1_bestL1pt   );
  tree_->Branch("Jpsi_m1_bestL1eta" , &Jpsi_m1_bestL1eta  );
  tree_->Branch("Jpsi_m1_bestL1phi" , &Jpsi_m1_bestL1phi  );
  tree_->Branch("Jpsi_m1_bestL1dR",   &Jpsi_m1_bestL1dR );

  tree_->Branch("Jpsi_m2_pt", &Jpsi_m2_pt );
  tree_->Branch("Jpsi_m2_eta", &Jpsi_m2_eta );
  tree_->Branch("Jpsi_m2_phi", &Jpsi_m2_phi );
  tree_->Branch("Jpsi_m2_mass", &Jpsi_m2_mass );
  tree_->Branch("Jpsi_m2_alsotag", &Jpsi_m2_alsotag );
  tree_->Branch("Jpsi_m2_q", &Jpsi_m2_q );
  tree_->Branch("Jpsi_m2_looseId"   , &Jpsi_m2_looseId    );
  tree_->Branch("Jpsi_m2_bestL1pt",   &Jpsi_m2_bestL1pt );
  tree_->Branch("Jpsi_m2_bestL1eta" , &Jpsi_m2_bestL1eta );
  tree_->Branch("Jpsi_m2_bestL1phi" , &Jpsi_m2_bestL1phi );
  tree_->Branch("Jpsi_m2_bestL1dR" ,  &Jpsi_m2_bestL1dR );

  // Analysis double-mu trigger
  tree_->Branch("Jpsi_m1_Dimu_dR",     &mu1_matchedDiMu_dR );
  tree_->Branch("Jpsi_m1_Dimu_pt",     &mu1_matchedDiMu_pt );
  tree_->Branch("Jpsi_m1_Dimu_eta",    &mu1_matchedDiMu_eta );
  tree_->Branch("Jpsi_m1_Dimu_phi",    &mu1_matchedDiMu_phi );
  tree_->Branch("Jpsi_m2_Dimu_dR",     &mu2_matchedDiMu_dR );
  tree_->Branch("Jpsi_m2_Dimu_pt",     &mu2_matchedDiMu_pt );
  tree_->Branch("Jpsi_m2_Dimu_eta",    &mu2_matchedDiMu_eta );
  tree_->Branch("Jpsi_m2_Dimu_phi",    &mu2_matchedDiMu_phi );
  tree_->Branch("DoubleMu_fired", &DoubleMu_fired );

  // Reference single-mu trigger
  tree_->Branch("Jpsi_m1_trgobj_pt",   &Jpsi_m1_trgobj_pt );
  tree_->Branch("Jpsi_m1_trgobj_eta",  &Jpsi_m1_trgobj_eta );
  tree_->Branch("Jpsi_m1_trgobj_phi",  &Jpsi_m1_trgobj_phi );
  tree_->Branch("Jpsi_m1_trgobj_mass", &Jpsi_m1_trgobj_mass );
  tree_->Branch("Jpsi_m1_trgobj_q",    &Jpsi_m1_trgobj_q );
  tree_->Branch("Jpsi_m1_trgobj_dR",   &Jpsi_m1_trgobj_dR );
  tree_->Branch("Jpsi_m2_trgobj_pt",   &Jpsi_m2_trgobj_pt );
  tree_->Branch("Jpsi_m2_trgobj_eta",  &Jpsi_m2_trgobj_eta );
  tree_->Branch("Jpsi_m2_trgobj_phi",  &Jpsi_m2_trgobj_phi );
  tree_->Branch("Jpsi_m2_trgobj_mass", &Jpsi_m2_trgobj_mass );
  tree_->Branch("Jpsi_m2_trgobj_q",    &Jpsi_m2_trgobj_q );
  tree_->Branch("Jpsi_m2_trgobj_dR",   &Jpsi_m2_trgobj_dR );

  // L1
  tree_->Branch("L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4", &L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4 );
  tree_->Branch("L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6", &L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6 );
  tree_->Branch("L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6", &L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6 );
  tree_->Branch("L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5", &L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5 );
  tree_->Branch("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 );
  tree_->Branch("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 );
  tree_->Branch("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2 );
  tree_->Branch("L1_DoubleMu4_SQ_OS_dR_Max1p2", &L1_DoubleMu4_SQ_OS_dR_Max1p2 );
  tree_->Branch("L1_OR", &L1_OR );

  // JPsi and daughters
  tree_->Branch("Jpsi_fit_pt",      &Jpsi_fit_pt );
  tree_->Branch("Jpsi_nonfit_pt",   &Jpsi_nonfit_pt );
  tree_->Branch("Jpsi_fit_eta",     &Jpsi_fit_eta );
  tree_->Branch("Jpsi_nonfit_eta",  &Jpsi_nonfit_eta );
  tree_->Branch("Jpsi_fit_phi",     &Jpsi_fit_phi );
  tree_->Branch("Jpsi_nonfit_phi",  &Jpsi_nonfit_phi);
  tree_->Branch("Jpsi_fit_mass",    &Jpsi_fit_mass );
  tree_->Branch("Jpsi_nonfit_mass", &Jpsi_nonfit_mass );
  tree_->Branch("Jpsi_fit_vprob",   &Jpsi_fit_vprob );
  tree_->Branch("Jpsi_muonsDr",     &Jpsi_muonsDr );
}

void NanoAnalyzerDoubleMu::reset(void){

  run = -1;
  event = -1;
  luminosityBlock = -1;
  nvtx = -1;

  Jpsi_m1_pt = -99;
  Jpsi_m1_eta = -99;
  Jpsi_m1_phi = -99;
  Jpsi_m1_mass = -99;
  Jpsi_m1_q = -99;
  Jpsi_m1_looseId = -99;
  Jpsi_m1_bestL1pt = -99;
  Jpsi_m1_bestL1eta = -99;
  Jpsi_m1_bestL1phi = -99;
  Jpsi_m1_bestL1dR = -99;

  Jpsi_m2_pt = -99;
  Jpsi_m2_eta = -99;
  Jpsi_m2_phi = -99;
  Jpsi_m2_mass = -99;
  Jpsi_m2_alsotag = -99;
  Jpsi_m2_q = -99;
  Jpsi_m2_looseId = -99;
  Jpsi_m2_bestL1pt = -99;
  Jpsi_m2_bestL1eta = -99;
  Jpsi_m2_bestL1phi = -99;
  Jpsi_m2_bestL1dR = -99;

  Jpsi_m1_trgobj_pt   = -99;
  Jpsi_m1_trgobj_eta  = -99;
  Jpsi_m1_trgobj_phi  = -99;
  Jpsi_m1_trgobj_mass = -99;
  Jpsi_m1_trgobj_q  = -99;
  Jpsi_m1_trgobj_dR = -99;

  Jpsi_m2_trgobj_pt   = -99;
  Jpsi_m2_trgobj_eta  = -99;
  Jpsi_m2_trgobj_phi  = -99;
  Jpsi_m2_trgobj_mass = -99;
  Jpsi_m2_trgobj_q  = -99;
  Jpsi_m2_trgobj_dR = -99;

  DoubleMu_fired = 0.;
  mu1_matchedDiMu_dR  = 99.;    
  mu2_matchedDiMu_dR  = 99.;
  mu1_matchedDiMu_pt  = 99.;
  mu2_matchedDiMu_pt  = 99.;
  mu1_matchedDiMu_eta = 99.;
  mu2_matchedDiMu_eta = 99.;
  mu1_matchedDiMu_phi = 99.;
  mu2_matchedDiMu_phi = 99.;

  L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4   = 99;
  L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6 = 99;
  L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6 = 99;
  L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5 = 99;
  L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4   = 99;
  L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4   = 99;
  L1_DoubleMu4p5_SQ_OS_dR_Max1p2      = 99;
  L1_DoubleMu4_SQ_OS_dR_Max1p2        = 99;
  L1_OR = 99;

  Jpsi_fit_pt      = -99;
  Jpsi_nonfit_pt   = -99;
  Jpsi_fit_eta     = -99;
  Jpsi_nonfit_eta  = -99;
  Jpsi_fit_phi     = -99;
  Jpsi_nonfit_phi  = -99;
  Jpsi_fit_mass    = -99;
  Jpsi_nonfit_mass = -99;
  Jpsi_fit_vprob   = -99;
  Jpsi_muonsDr     = -99;

  // Not for tree
  dimuobj1_pt  = -999.;
  dimuobj1_eta = -999.;
  dimuobj1_phi = -999.;
  dimuobj2_pt  = -999.;
  dimuobj2_eta = -999.;
  dimuobj2_phi = -999.;
}

DEFINE_FWK_MODULE(NanoAnalyzerDoubleMu);
