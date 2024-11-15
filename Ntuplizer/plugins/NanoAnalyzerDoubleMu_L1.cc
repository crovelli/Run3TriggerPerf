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
class NanoAnalyzerDoubleMu_L1 : public edm::one::EDAnalyzer<>
{
public:
  explicit NanoAnalyzerDoubleMu_L1(const edm::ParameterSet&);
  ~NanoAnalyzerDoubleMu_L1();
  
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
  edm::EDGetTokenT<l1t::MuonBxCollection> l1MU_;      
  
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
  float  Jpsi_m1_abseta;
  float  Jpsi_m1_phi;
  int    Jpsi_m1_q;   
  float  Jpsi_m1_looseId;   
  float  Jpsi_m1_mediumId;
  float  Jpsi_m1_tightId;
  //
  float  Jpsi_m1_bestL1pt;
  float  Jpsi_m1_bestL1eta;
  float  Jpsi_m1_bestL1phi;
  float  Jpsi_m1_bestL1dR;

  float  Jpsi_m2_pt      ;
  float  Jpsi_m2_eta     ;
  float  Jpsi_m2_abseta;
  float  Jpsi_m2_phi     ;
  int    Jpsi_m2_alsotag ;
  int    Jpsi_m2_q   ;   
  float  Jpsi_m2_looseId  ;
  float  Jpsi_m2_mediumId;
  float  Jpsi_m2_tightId;
  // 
  float  Jpsi_m2_bestL1pt ;
  float  Jpsi_m2_bestL1eta ;
  float  Jpsi_m2_bestL1phi ;
  float  Jpsi_m2_bestL1dR ;

  float  Jpsi_m1_trgobj_pt      ;
  float  Jpsi_m1_trgobj_eta     ;
  float  Jpsi_m1_trgobj_phi     ;
  int    Jpsi_m1_trgobj_q      ;   
  float  Jpsi_m1_trgobj_dR      ;
  float  Jpsi_m2_trgobj_pt ;
  float  Jpsi_m2_trgobj_eta ;
  float  Jpsi_m2_trgobj_phi ;
  int    Jpsi_m2_trgobj_q ;   
  float  Jpsi_m2_trgobj_dR ;
  //
  int probePathFired;
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
  float  Jpsi_muonsDz;
  
  // Not for tree, other variables declaration
  helperfunc aux;
  float chi = 0.;
  float ndf = 0.;

  TH1F * hist; 

}; // end of class member


NanoAnalyzerDoubleMu_L1::NanoAnalyzerDoubleMu_L1(const edm::ParameterSet& iConfig): 
  bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>())
{
  muonToken_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));

  verticeToken_       = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));

  triggerToken_	      = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLT"));
  triggerobjectToken_ = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerobjects"));

  l1MU_ = consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("l1MU"));

  hist = fs->make<TH1F>("cutflow", "cutflow", 10,0,10);
  tree_ = fs->make<TTree>( "tree", "tree" );

  createBranch();

} // end of constructor

NanoAnalyzerDoubleMu_L1::~NanoAnalyzerDoubleMu_L1() { }

void NanoAnalyzerDoubleMu_L1::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  reset();

  const auto& bField = iSetup.getData(bFieldToken_);

  iEvent.getByToken(verticeToken_, vertices_ );
  iEvent.getByToken(triggerToken_, HLTtriggers_);
  iEvent.getByToken(triggerobjectToken_ , triggerObjects); 
  
  nvtx = vertices_->size();
  run = (iEvent.id()).run();
  event = (iEvent.id()).event();

  // 1st PV 
  const reco::Vertex &PV = vertices_->front(); 
  
  // check if the reference trigger and/or the analysis trigger fired
  bool tagPathFired = false;
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*HLTtriggers_);
  std::string tagTriggerName="";    
  std::string probeTriggerName="";  

  for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {

    ///////std::cout << "Loop over all triggers: i = " << i << ", name = " << trigNames.triggerName(i) << ", fired = " << HLTtriggers_->accept(i) << std::endl;
    
    // Check if the tag trigger path fired
    if(trigNames.triggerName(i).find("HLT_Mu8_v")!= std::string::npos){             // chiara
      if(HLTtriggers_->accept(i)){
	tagPathFired = true;
	tagTriggerName=trigNames.triggerName(i);  
      }
    }
    
    // Check if the wanted trigger path firef
    if(trigNames.triggerName(i).find("HLT_Mu0_L1DoubleMu_v")!= std::string::npos){ 
      if(HLTtriggers_->accept(i)){
	probePathFired = true;
	probeTriggerName=trigNames.triggerName(i); 
      }
    }
  }

  /*
  std::cout << std::endl;
  std::cout << "Ref trigger? " << tagPathFired << std::endl;
  std::cout << "AN trigger? "  << probePathFired << std::endl;
  std::cout << "probeTriggerName = " << probeTriggerName << std::endl;
  std::cout << std::endl;
  */

  // Reference path fired - skip otherwise
  if(!tagPathFired) return;    
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
    if (fabs(muon.eta()) > 2.4) continue;
    if (!(muon.track().isNonnull())) continue;
    
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
	if(pathNamesAll[h]==tagTriggerName) isPathExist = true;     
      }
      if(!isPathExist) continue;
      
      int muObjNumber = -1;
      for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){	

	//std::cout << "Event: Filter " << hh << " => " << obj.filterLabels()[hh] << " ";
	//std::cout << "" << std::endl;
	  
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
	//std::cout<< "DeltaR = " << deltaR << std::endl;
	//std::cout << "This is a muon HLT candidate" << endl;
	if(deltaR < 0.3){    // chiara
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

  // Kinematic fit
  const reco::TransientTrack muon1TT((*(muoncollection[mcidx_m1].bestTrack())),&bField);  
  const reco::TransientTrack muon2TT((*(muoncollection[mcidx_m2].bestTrack())),&bField);
  KinVtxFitter fitter(
		      { muon1TT, muon2TT },
		      { aux.muon_mass, aux.muon_mass},
		      { aux.muon_sigma, aux.muon_sigma } //some small sigma for the particle mass
		      );
  if ( !fitter.success() ) return;  
  hist->Fill(6);

  KinematicState fitted_cand = fitter.fitted_candidate();
  RefCountedKinematicVertex fitted_vtx = fitter.fitted_refvtx();  
  Jpsi_fit_pt    = fitted_cand.globalMomentum().perp();
  Jpsi_fit_eta   = fitted_cand.globalMomentum().eta();
  Jpsi_fit_phi   = fitted_cand.globalMomentum().phi();
  Jpsi_fit_mass  = fitted_cand.mass();
  Jpsi_fit_vprob = fitter.prob();
  Jpsi_muonsDr   = mu1TV3.DeltaR(mu2TV3);
  Jpsi_muonsDz   = fabs(muoncollection[mcidx_m1].vz() - muoncollection[mcidx_m2].vz() );
  
  // Match L1 / offline
  float bestMatchM1_eta = -99.;
  float bestMatchM1_phi = -99.;
  float bestMatchM1_pt  = -99.;
  float bestMatchM1_dR  = 999.;
  float bestMatchM2_eta = -99.;
  float bestMatchM2_phi = -99.;
  float bestMatchM2_pt  = -99.;
  float bestMatchM2_dR  = 999.;
  const auto &l1MU = iEvent.get(l1MU_); 
  for (l1t::MuonBxCollection::const_iterator it = l1MU.begin(0); it != l1MU.end(0); it++) {
    pat::TriggerObjectStandAlone l1obj(it->p4());
    
    TVector3 l1objTV3;
    l1objTV3.SetPtEtaPhi( it->pt(), it->eta(), it->phi() );
    Float_t deltaRM1 = fabs(mu1TV3.DeltaR(l1objTV3));
    Float_t deltaRM2 = fabs(mu2TV3.DeltaR(l1objTV3));
    
    if (deltaRM1<bestMatchM1_dR){
      bestMatchM1_eta  = it->eta(); 
      bestMatchM1_phi  = it->phi(); 
      bestMatchM1_pt   = it->pt();    
      bestMatchM1_dR   = deltaRM1;
    }
    if (deltaRM2<bestMatchM2_dR){
      bestMatchM2_eta  = it->eta(); 
      bestMatchM2_phi  = it->phi(); 
      bestMatchM2_pt   = it->pt(); 
      bestMatchM2_dR   = deltaRM2;
    }
  }

  // Infos about JPsi candidate
  Jpsi_m1_pt   = muoncollection[mcidx_m1].pt();
  Jpsi_m1_eta  = muoncollection[mcidx_m1].eta();
  Jpsi_m1_abseta = fabs(muoncollection[mcidx_m1].eta());
  Jpsi_m1_phi  = muoncollection[mcidx_m1].phi();
  Jpsi_m1_q    = muoncollection[mcidx_m1].charge();
  Jpsi_m1_looseId  = muoncollection[mcidx_m1].isLooseMuon();
  Jpsi_m1_mediumId = muoncollection[mcidx_m1].isMediumMuon();
  Jpsi_m1_tightId  = muoncollection[mcidx_m1].isTightMuon(PV);
  Jpsi_m1_bestL1pt   = bestMatchM1_pt;
  Jpsi_m1_bestL1eta  = bestMatchM1_eta;
  Jpsi_m1_bestL1phi  = bestMatchM1_phi;
  Jpsi_m1_bestL1dR   = bestMatchM1_dR;  

  Jpsi_m2_pt   = muoncollection[mcidx_m2].pt();
  Jpsi_m2_eta  = muoncollection[mcidx_m2].eta();
  Jpsi_m2_abseta = fabs(muoncollection[mcidx_m2].eta());
  Jpsi_m2_phi  = muoncollection[mcidx_m2].phi();
  Jpsi_m2_alsotag = is_probe_also_tag;
  Jpsi_m2_q    = muoncollection[mcidx_m2].charge();
  Jpsi_m2_looseId  = muoncollection[mcidx_m2].isLooseMuon();
  Jpsi_m2_mediumId = muoncollection[mcidx_m2].isMediumMuon();
  Jpsi_m2_tightId  = muoncollection[mcidx_m2].isTightMuon(PV);
  Jpsi_m2_bestL1pt  = bestMatchM2_pt;
  Jpsi_m2_bestL1eta = bestMatchM2_eta;
  Jpsi_m2_bestL1phi = bestMatchM2_phi;
  Jpsi_m2_bestL1dR  = bestMatchM2_dR;  
  
  if (mcidx_trgobj1>=0) {
    
    Jpsi_m1_trgobj_pt   = trg_obj_collection[mcidx_trgobj1].pt();
    Jpsi_m1_trgobj_eta  = trg_obj_collection[mcidx_trgobj1].eta();
    Jpsi_m1_trgobj_phi  = trg_obj_collection[mcidx_trgobj1].phi();
    Jpsi_m1_trgobj_q    = trg_obj_collection[mcidx_trgobj1].charge();
    TVector3 trgobj1TV3;
    trgobj1TV3.SetPtEtaPhi( Jpsi_m1_trgobj_pt, Jpsi_m1_trgobj_eta, Jpsi_m1_trgobj_phi);
    Jpsi_m1_trgobj_dR = fabs(mu1TV3.DeltaR(trgobj1TV3));
  }

  if (mcidx_trgobj2>=0) {

    Jpsi_m2_trgobj_pt   = trg_obj_collection[mcidx_trgobj2].pt();
    Jpsi_m2_trgobj_eta  = trg_obj_collection[mcidx_trgobj2].eta();
    Jpsi_m2_trgobj_phi  = trg_obj_collection[mcidx_trgobj2].phi();
    Jpsi_m2_trgobj_q    = trg_obj_collection[mcidx_trgobj2].charge();
    TVector3 trgobj2TV3;
    trgobj2TV3.SetPtEtaPhi( Jpsi_m2_trgobj_pt, Jpsi_m2_trgobj_eta, Jpsi_m2_trgobj_phi);
    Jpsi_m2_trgobj_dR = fabs(mu2TV3.DeltaR(trgobj2TV3));
  }

  Jpsi_nonfit_pt   = jpsi_tlv_highest.Pt();
  Jpsi_nonfit_eta  = jpsi_tlv_highest.Eta();
  Jpsi_nonfit_phi  = jpsi_tlv_highest.Phi();
  Jpsi_nonfit_mass = jpsi_tlv_highest.M(); 

  tree_->Fill();

  return;

} //NanoAnalyzer::analyze ends



//**************************************************
//************* additional methods *****************
//**************************************************

void NanoAnalyzerDoubleMu_L1::beginJob(const edm::ParameterSet& iConfig) { }

void NanoAnalyzerDoubleMu_L1::beginRun(const edm::Run &iRun, const edm::EventSetup &iStp) { }

void NanoAnalyzerDoubleMu_L1::fillDescriptions(edm::ConfigurationDescriptions & descriptions) { }

void NanoAnalyzerDoubleMu_L1::endRun(edm::Run const&, edm::EventSetup const&) { }

void NanoAnalyzerDoubleMu_L1::endJob() { }

//define this as a plug-in

// branch title creation
void NanoAnalyzerDoubleMu_L1::createBranch() { 

  tree_->Branch("run", &run, "run/i");
  tree_->Branch("event", &event, "event/l");
  tree_->Branch("nvtx", &nvtx, "nvtx/i");

  // Offline
  tree_->Branch("Jpsi_m1_pt", &Jpsi_m1_pt );
  tree_->Branch("Jpsi_m1_eta", &Jpsi_m1_eta );
  tree_->Branch("Jpsi_m1_abseta", &Jpsi_m1_abseta );
  tree_->Branch("Jpsi_m1_phi", &Jpsi_m1_phi );
  tree_->Branch("Jpsi_m1_q", &Jpsi_m1_q );
  tree_->Branch("Jpsi_m1_looseId"   , &Jpsi_m1_looseId    );
  tree_->Branch("Jpsi_m1_mediumId"  , &Jpsi_m1_mediumId   );
  tree_->Branch("Jpsi_m1_tightId"   , &Jpsi_m1_tightId   );
  tree_->Branch("Jpsi_m1_bestL1pt",   &Jpsi_m1_bestL1pt   );
  tree_->Branch("Jpsi_m1_bestL1eta" , &Jpsi_m1_bestL1eta  );
  tree_->Branch("Jpsi_m1_bestL1phi" , &Jpsi_m1_bestL1phi  );
  tree_->Branch("Jpsi_m1_bestL1dR",   &Jpsi_m1_bestL1dR );

  tree_->Branch("Jpsi_m2_pt", &Jpsi_m2_pt );
  tree_->Branch("Jpsi_m2_eta", &Jpsi_m2_eta );
  tree_->Branch("Jpsi_m2_abseta", &Jpsi_m2_abseta );
  tree_->Branch("Jpsi_m2_phi", &Jpsi_m2_phi );
  tree_->Branch("Jpsi_m2_alsotag", &Jpsi_m2_alsotag );
  tree_->Branch("Jpsi_m2_q", &Jpsi_m2_q );
  tree_->Branch("Jpsi_m2_looseId"   , &Jpsi_m2_looseId  );
  tree_->Branch("Jpsi_m2_mediumId"  , &Jpsi_m2_mediumId   );
  tree_->Branch("Jpsi_m2_tightId"   , &Jpsi_m2_tightId   );
  tree_->Branch("Jpsi_m2_bestL1pt",   &Jpsi_m2_bestL1pt );
  tree_->Branch("Jpsi_m2_bestL1eta" , &Jpsi_m2_bestL1eta );
  tree_->Branch("Jpsi_m2_bestL1phi" , &Jpsi_m2_bestL1phi );
  tree_->Branch("Jpsi_m2_bestL1dR" ,  &Jpsi_m2_bestL1dR );

  // Analysis double-mu trigger
  tree_->Branch("probePathFired", &probePathFired );

  // Reference single-mu trigger
  tree_->Branch("Jpsi_m1_trgobj_pt",   &Jpsi_m1_trgobj_pt );
  tree_->Branch("Jpsi_m1_trgobj_eta",  &Jpsi_m1_trgobj_eta );
  tree_->Branch("Jpsi_m1_trgobj_phi",  &Jpsi_m1_trgobj_phi );
  tree_->Branch("Jpsi_m1_trgobj_q",    &Jpsi_m1_trgobj_q );
  tree_->Branch("Jpsi_m1_trgobj_dR",   &Jpsi_m1_trgobj_dR );

  tree_->Branch("Jpsi_m2_trgobj_pt",   &Jpsi_m2_trgobj_pt );
  tree_->Branch("Jpsi_m2_trgobj_eta",  &Jpsi_m2_trgobj_eta );
  tree_->Branch("Jpsi_m2_trgobj_phi",  &Jpsi_m2_trgobj_phi );
  tree_->Branch("Jpsi_m2_trgobj_q",    &Jpsi_m2_trgobj_q );
  tree_->Branch("Jpsi_m2_trgobj_dR",   &Jpsi_m2_trgobj_dR );

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
  tree_->Branch("Jpsi_muonsDz",     &Jpsi_muonsDz );
}

void NanoAnalyzerDoubleMu_L1::reset(void){

  run = -1;
  event = -1;
  nvtx = -1;

  Jpsi_m1_pt = -99;
  Jpsi_m1_eta = -99;
  Jpsi_m1_abseta = -99;
  Jpsi_m1_phi = -99;
  Jpsi_m1_q = -99;
  Jpsi_m1_looseId = -99;
  Jpsi_m1_mediumId = -99;
  Jpsi_m1_tightId = -99;
  Jpsi_m1_bestL1pt = -99;
  Jpsi_m1_bestL1eta = -99;
  Jpsi_m1_bestL1phi = -99;
  Jpsi_m1_bestL1dR = -99;

  Jpsi_m2_pt = -99;
  Jpsi_m2_eta = -99;
  Jpsi_m2_abseta = -99;
  Jpsi_m2_phi = -99;
  Jpsi_m2_alsotag = -99;
  Jpsi_m2_q = -99;
  Jpsi_m2_looseId = -99;
  Jpsi_m2_mediumId = -99;
  Jpsi_m2_tightId = -99;
  Jpsi_m2_bestL1pt = -99;
  Jpsi_m2_bestL1eta = -99;
  Jpsi_m2_bestL1phi = -99;
  Jpsi_m2_bestL1dR = -99;

  Jpsi_m1_trgobj_pt   = -99;
  Jpsi_m1_trgobj_eta  = -99;
  Jpsi_m1_trgobj_phi  = -99;
  Jpsi_m1_trgobj_q  = -99;
  Jpsi_m1_trgobj_dR = -99;

  Jpsi_m2_trgobj_pt   = -99;
  Jpsi_m2_trgobj_eta  = -99;
  Jpsi_m2_trgobj_phi  = -99;
  Jpsi_m2_trgobj_q  = -99;
  Jpsi_m2_trgobj_dR = -99;

  probePathFired = 0.;

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
  Jpsi_muonsDz     = -99;
}

DEFINE_FWK_MODULE(NanoAnalyzerDoubleMu_L1);
