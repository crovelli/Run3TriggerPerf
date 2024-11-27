// system include files
#include <memory>
#include <iostream>
#include <vector>

// Root
#include "TMath.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"

// Utilities
#include "../interface/helperfunc.h"


//*****************************
// general user include files *
//*****************************
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"

//**************************
// for trigger information *
//**************************
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

//***************************
// for tracking information *
//***************************
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//*************************
// for vertex information *
//*************************
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/VertexReco/interface/Vertex.h>

//***********************
// for muon information *
//***********************
#include "DataFormats/PatCandidates/interface/Muon.h"

//***********************
// extra                * 
//***********************   
#include "KinVtxFitter.h"   


// set namespaces
using namespace edm;
using namespace reco;
using namespace std;


//***********************************
class NanoAnalyzerDoubleMu_tagPerL1 : public edm::one::EDAnalyzer<>
{
public:
  explicit NanoAnalyzerDoubleMu_tagPerL1(const edm::ParameterSet&);
  ~NanoAnalyzerDoubleMu_tagPerL1();
  
  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
  
private:

  virtual void beginJob(const edm::ParameterSet& iConfig);
  virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iStp);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void endJob();

  void createBranch();

  void reset();

  EDGetTokenT<reco::VertexCollection> verticeToken_;
  EDGetTokenT<pat::MuonCollection> muonToken_;
  EDGetTokenT<edm::TriggerResults> triggerToken_;  
  EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobjectToken_;
  
  Handle<pat::MuonCollection> muons_;
  Handle< reco::VertexCollection > vertices_;
  Handle< edm::TriggerResults> HLTtriggers_;
  Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  
/////////////////////////////////////////////////////////////////////////
  edm::Service<TFileService> fs;
  TTree* tree_;

  /// original nanos
  UInt_t run;
  ULong64_t event;
  int nvtx;

  float  Jpsi_m1_pt;
  float  Jpsi_m1_eta;
  float  Jpsi_m1_mediumId;
  float  Jpsi_m1_tightId;

  float  Jpsi_m2_pt      ;
  float  Jpsi_m2_eta     ;
  float  Jpsi_m2_mediumId;
  float  Jpsi_m2_tightId;

  float  Jpsi_m1_trgobj_pt      ;
  float  Jpsi_m1_trgobj_eta     ;
  float  Jpsi_m1_trgobj_dR      ;
  float  Jpsi_m2_trgobj_pt ;
  float  Jpsi_m2_trgobj_eta ;
  float  Jpsi_m2_trgobj_dR ;

  int tagPathFired;
  
  float  Jpsi_nonfit_mass;
  float  Jpsi_muonsDr;
  float  Jpsi_muonsDz;
  
  // Not for tree, other variables declaration
  helperfunc aux;
  float chi = 0.;
  float ndf = 0.;

  TH1F * hist; 

}; // end of class member


NanoAnalyzerDoubleMu_tagPerL1::NanoAnalyzerDoubleMu_tagPerL1(const edm::ParameterSet& iConfig)
{
  muonToken_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  
  verticeToken_       = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));

  triggerToken_	      = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLT"));
  triggerobjectToken_ = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerobjects"));

  hist = fs->make<TH1F>("cutflow", "cutflow", 10,0,10);
  tree_ = fs->make<TTree>( "tree", "tree" );
  
  createBranch();

} // end of constructor

NanoAnalyzerDoubleMu_tagPerL1::~NanoAnalyzerDoubleMu_tagPerL1() { }

void NanoAnalyzerDoubleMu_tagPerL1::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  reset();

  iEvent.getByToken(verticeToken_, vertices_ );
  iEvent.getByToken(triggerToken_, HLTtriggers_);
  iEvent.getByToken(triggerobjectToken_ , triggerObjects); 
  
  nvtx = vertices_->size();
  run = (iEvent.id()).run();
  event = (iEvent.id()).event();

  // 1st PV 
  const reco::Vertex &PV = vertices_->front(); 
  
  // check if the reference trigger fired
  tagPathFired = false;
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*HLTtriggers_);
  std::string tagTriggerName="";    
  for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {
    if(trigNames.triggerName(i).find("HLT_Mu8_v")!= std::string::npos){    
      if(HLTtriggers_->accept(i)){
	tagPathFired = true;
	tagTriggerName=trigNames.triggerName(i);  
      }
    }
  }

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

  // Offline muons
  for(size_t imuon = 0; imuon < muons_->size(); imuon++){
    const pat::Muon & muon = (*muons_)[imuon];

    // Offline cuts
    if (fabs(muon.eta()) > 2.4) continue;
    if (!(muon.track().isNonnull())) continue;
    
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
	if(pathNamesAll[h]==tagTriggerName) isPathExist = true;     
      }
      if(!isPathExist) continue;
      
      int muObjNumber = -1;
      for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){	
	if(obj.filterLabels()[hh].find("hltL3fL1sMu5L1f0L2f5L3Filtered8") != std::string::npos) {  // chiara  HLT_Mu8
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
	if(deltaR < 0.3){    // chiara
	  trigObjMatchMu = true;
	  if (deltaR < best_match_dR){
	    best_match_dR = deltaR;
	    best_match_obj = obj;
	  }
	}
      }
    } // Loop over trigger object
    
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

  // At least 2 offline muons
  if (muoncollection.size() < 2) return; 
  hist->Fill(1);
  
  // Prepare offline muon pairs
  float jpsi_max_pt = -1;
  int mcidx_m1 = -1;
  int mcidx_m2 = -1;
  int mcidx_trgobj1 = -1;
  int mcidx_trgobj2 = -1;
  TLorentzVector jpsi_tlv_highest;

  for(int im = 0; im < (int)muoncollection.size(); im++){
    for(int jm = im+1; jm < (int)muoncollection.size(); jm++){
      
      // match with trigger muon
      int this_mcidx_trgobj1 = -1;
      int this_mcidx_trgobj2 = -1;
      if (muonmatched[im]>-999) 
	this_mcidx_trgobj1 = muonmatched[im];
      if (muonmatched[jm]>-999)
	this_mcidx_trgobj2 = muonmatched[jm];
      
      const pat::Muon mu1 = muoncollection[im];
      const pat::Muon mu2 = muoncollection[jm];

      TLorentzVector tlv_m1;
      TLorentzVector tlv_m2;
      tlv_m1.SetPtEtaPhiM(mu1.pt(), mu1.eta(), mu1.phi(), aux.mass_muon);
      tlv_m2.SetPtEtaPhiM(mu2.pt(), mu2.eta(), mu2.phi(), aux.mass_muon);
      
      TLorentzVector tlv_jpsi = (tlv_m1 + tlv_m2);
      float jpsi_mass = tlv_jpsi.M();
      float jpsi_pt = tlv_jpsi.Pt();

      // if (mu1.charge() + mu2.charge() !=0) continue;
      // if (jpsi_mass < 2.0) continue; 
      // if (jpsi_mass > 4.0) continue;
      
      if(jpsi_max_pt < jpsi_pt){
	jpsi_max_pt = jpsi_pt;
	mcidx_m1 = im;
	mcidx_m2 = jm;
	mcidx_trgobj1 = this_mcidx_trgobj1;
	mcidx_trgobj2 = this_mcidx_trgobj2;
	jpsi_tlv_highest = tlv_jpsi;
      }
    }
  }

  // At least 1 reco J/psi
  if(jpsi_max_pt == -1) return;  
  hist->Fill(2);

  // Swap 1<->2 so that: 1 is always tag, 2 is always probe
  int mcidx_tag = -99;
  int mcidx_pro = -99;
  int mcidx_trgobj_tag = -99;
  int mcidx_trgobj_pro = -99;
  if (mcidx_trgobj1>=0) {
    mcidx_pro = mcidx_m2;
    mcidx_trgobj_pro = mcidx_trgobj2;
    mcidx_tag = mcidx_m1;
    mcidx_trgobj_tag = mcidx_trgobj1;
  } else if (mcidx_trgobj2>=0) { 
    mcidx_pro = mcidx_m1;
    mcidx_trgobj_pro = mcidx_trgobj1;
    mcidx_tag = mcidx_m2;
    mcidx_trgobj_tag = mcidx_trgobj2;
  } else {
    mcidx_pro = mcidx_m1;
    mcidx_tag = mcidx_m2;
  }
  
  mcidx_m1 = mcidx_tag;
  mcidx_m2 = mcidx_pro;
  mcidx_trgobj1 = mcidx_trgobj_tag;
  mcidx_trgobj2 = mcidx_trgobj_pro;

  
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
  Jpsi_muonsDr   = mu1TV3.DeltaR(mu2TV3);
  Jpsi_muonsDz   = fabs(muoncollection[mcidx_m1].vz() - muoncollection[mcidx_m2].vz() );
  
  // Infos about JPsi candidate
  Jpsi_m1_pt   = muoncollection[mcidx_m1].pt();
  Jpsi_m1_eta  = muoncollection[mcidx_m1].eta();
  Jpsi_m1_mediumId = muoncollection[mcidx_m1].isMediumMuon();
  Jpsi_m1_tightId  = muoncollection[mcidx_m1].isTightMuon(PV);

  Jpsi_m2_pt   = muoncollection[mcidx_m2].pt();
  Jpsi_m2_eta  = muoncollection[mcidx_m2].eta();
  Jpsi_m2_mediumId = muoncollection[mcidx_m2].isMediumMuon();
  Jpsi_m2_tightId  = muoncollection[mcidx_m2].isTightMuon(PV);

  if (mcidx_trgobj1>=0) {
    
    Jpsi_m1_trgobj_pt   = trg_obj_collection[mcidx_trgobj1].pt();
    Jpsi_m1_trgobj_eta  = trg_obj_collection[mcidx_trgobj1].eta();
    float Jpsi_m1_trgobj_phi = trg_obj_collection[mcidx_trgobj1].phi();
    TVector3 trgobj1TV3;
    trgobj1TV3.SetPtEtaPhi( Jpsi_m1_trgobj_pt, Jpsi_m1_trgobj_eta, Jpsi_m1_trgobj_phi);
    Jpsi_m1_trgobj_dR = fabs(mu1TV3.DeltaR(trgobj1TV3));
  }

  if (mcidx_trgobj2>=0) {

    Jpsi_m2_trgobj_pt   = trg_obj_collection[mcidx_trgobj2].pt();
    Jpsi_m2_trgobj_eta  = trg_obj_collection[mcidx_trgobj2].eta();
    float Jpsi_m2_trgobj_phi = trg_obj_collection[mcidx_trgobj2].phi();
    TVector3 trgobj2TV3;
    trgobj2TV3.SetPtEtaPhi( Jpsi_m2_trgobj_pt, Jpsi_m2_trgobj_eta, Jpsi_m2_trgobj_phi);
    Jpsi_m2_trgobj_dR = fabs(mu2TV3.DeltaR(trgobj2TV3));
  }
  
  Jpsi_nonfit_mass = jpsi_tlv_highest.M(); 

  tree_->Fill();

  return;

} //NanoAnalyzer::analyze ends



//**************************************************
//************* additional methods *****************
//**************************************************

void NanoAnalyzerDoubleMu_tagPerL1::beginJob(const edm::ParameterSet& iConfig) { }

void NanoAnalyzerDoubleMu_tagPerL1::beginRun(const edm::Run &iRun, const edm::EventSetup &iStp) { }

void NanoAnalyzerDoubleMu_tagPerL1::fillDescriptions(edm::ConfigurationDescriptions & descriptions) { }

void NanoAnalyzerDoubleMu_tagPerL1::endRun(edm::Run const&, edm::EventSetup const&) { }

void NanoAnalyzerDoubleMu_tagPerL1::endJob() { }

//define this as a plug-in

// branch title creation
void NanoAnalyzerDoubleMu_tagPerL1::createBranch() { 

  tree_->Branch("run",   &run,   "run/i");
  tree_->Branch("event", &event, "event/l");
  tree_->Branch("nvtx",  &nvtx,  "nvtx/i");

  // Offline
  tree_->Branch("Jpsi_m1_pt",  &Jpsi_m1_pt );
  tree_->Branch("Jpsi_m1_eta", &Jpsi_m1_eta );
  tree_->Branch("Jpsi_m1_mediumId", &Jpsi_m1_mediumId   );
  tree_->Branch("Jpsi_m1_tightId" , &Jpsi_m1_tightId   );
  tree_->Branch("Jpsi_m2_pt",  &Jpsi_m2_pt );
  tree_->Branch("Jpsi_m2_eta", &Jpsi_m2_eta );
  tree_->Branch("Jpsi_m2_mediumId", &Jpsi_m2_mediumId   );
  tree_->Branch("Jpsi_m2_tightId" , &Jpsi_m2_tightId   );

  // Single-mu trigger
  tree_->Branch("tagPathFired",        &tagPathFired );   
  tree_->Branch("Jpsi_m1_trgobj_pt",   &Jpsi_m1_trgobj_pt );
  tree_->Branch("Jpsi_m1_trgobj_eta",  &Jpsi_m1_trgobj_eta );
  tree_->Branch("Jpsi_m1_trgobj_dR",   &Jpsi_m1_trgobj_dR );
  tree_->Branch("Jpsi_m2_trgobj_pt",   &Jpsi_m2_trgobj_pt );
  tree_->Branch("Jpsi_m2_trgobj_eta",  &Jpsi_m2_trgobj_eta );
  tree_->Branch("Jpsi_m2_trgobj_dR",   &Jpsi_m2_trgobj_dR );
  
  // JPsi and daughters
  tree_->Branch("Jpsi_nonfit_mass", &Jpsi_nonfit_mass );
  tree_->Branch("Jpsi_muonsDr",     &Jpsi_muonsDr );
  tree_->Branch("Jpsi_muonsDz",     &Jpsi_muonsDz );
}

void NanoAnalyzerDoubleMu_tagPerL1::reset(void){

  run = -1;
  event = -1;
  nvtx = -1;

  Jpsi_m1_pt = -99;
  Jpsi_m1_eta = -99;
  Jpsi_m1_mediumId = -99;
  Jpsi_m1_tightId = -99;

  Jpsi_m2_pt = -99;
  Jpsi_m2_eta = -99;
  Jpsi_m2_mediumId = -99;
  Jpsi_m2_tightId = -99;

  tagPathFired = 0.;   
  
  Jpsi_m1_trgobj_pt   = -99;
  Jpsi_m1_trgobj_eta  = -99;
  Jpsi_m1_trgobj_dR = -99;

  Jpsi_m2_trgobj_pt   = -99;
  Jpsi_m2_trgobj_eta  = -99;
  Jpsi_m2_trgobj_dR = -99;

  Jpsi_nonfit_mass = -99;
  Jpsi_muonsDr     = -99;
  Jpsi_muonsDz     = -99;
}

DEFINE_FWK_MODULE(NanoAnalyzerDoubleMu_tagPerL1);
