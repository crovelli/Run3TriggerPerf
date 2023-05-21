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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

// set namespaces
using namespace edm;
using namespace reco;
using namespace std;


//***********************************
// main analyzer class (EDAnalyzer) *
//***********************************

#ifndef CMSSW12plus
// Run 1 and 2
class NanoAnalyzerDoubleMuMC : public edm::EDAnalyzer
#else
// Run 3
class NanoAnalyzerDoubleMuMC : public edm::one::EDAnalyzer<>
#endif
{
public:
  explicit NanoAnalyzerDoubleMuMC(const edm::ParameterSet&);
  ~NanoAnalyzerDoubleMuMC();

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
  edm::EDGetTokenT<BXVector<l1t::Muon>> l1MU_;      
  EDGetTokenT<pat::PackedGenParticleCollection> packedgenparticleToken_;

  Handle<pat::MuonCollection> muons_;
  Handle< reco::VertexCollection > vertices_;
  Handle< edm::TriggerResults> HLTtriggers_;
  Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  Handle<pat::PackedGenParticleCollection> packedgenparticles_   ;

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
  float  Jpsi_m1_genMatchPt;
  float  Jpsi_m1_genMatchEta;
  float  Jpsi_m1_genMatchPhi;
  float  Jpsi_m1_genMatchDR;
  //
  float  Jpsi_m2_pt;
  float  Jpsi_m2_eta;
  float  Jpsi_m2_phi;
  float  Jpsi_m2_mass;
  int    Jpsi_m2_q;   
  float  Jpsi_m2_looseId;   
  float  Jpsi_m2_bestL1pt;
  float  Jpsi_m2_bestL1eta;
  float  Jpsi_m2_bestL1phi;
  float  Jpsi_m2_bestL1dR;  
  float  Jpsi_m2_genMatchPt;
  float  Jpsi_m2_genMatchEta;
  float  Jpsi_m2_genMatchPhi;
  float  Jpsi_m2_genMatchDR;
  //
  int DoubleMu_fired;
  float mu1_matchedDiMu_dR  = 99.;
  float mu2_matchedDiMu_dR  = 99.;
  float mu1_matchedDiMu_pt  = 99.;
  float mu2_matchedDiMu_pt  = 99.;
  float mu1_matchedDiMu_eta = 99.;
  float mu2_matchedDiMu_eta = 99.;
  float mu1_matchedDiMu_phi = 99.;
  float mu2_matchedDiMu_phi = 99.;
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
  float  Jpsi_nonfit_pt;
  float  Jpsi_nonfit_eta;
  float  Jpsi_nonfit_phi;
  float  Jpsi_nonfit_mass;
  float  Jpsi_muonsDr;
  //
  float  JpsiGen_pt;
  float  JpsiGen_eta;
  float  JpsiGen_phi;
  float  JpsiGen_mass;

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

NanoAnalyzerDoubleMuMC::NanoAnalyzerDoubleMuMC(const edm::ParameterSet& iConfig) {
  muonToken_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));

  verticeToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));

  triggerToken_	= consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLT"));
  triggerobjectToken_ = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerobjects"));

  l1MU_ = consumes<BXVector<l1t::Muon>>(iConfig.getParameter<edm::InputTag>("l1MU"));      

  packedgenparticleToken_ = consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedgenparticles"));

  hist = fs->make<TH1F>("cutflow", "cutflow", 10,0,10);
  tree_ = fs->make<TTree>( "tree", "tree" );

  createBranch();

} // end of constructor

NanoAnalyzerDoubleMuMC::~NanoAnalyzerDoubleMuMC() { }

void NanoAnalyzerDoubleMuMC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  reset();
  
  iEvent.getByToken(verticeToken_, vertices_ );   
  iEvent.getByToken(triggerToken_, HLTtriggers_ );
  iEvent.getByToken(triggerobjectToken_, triggerObjects );

  nvtx = vertices_->size();
  
  run = (iEvent.id()).run();
  event = (iEvent.id()).event();
  luminosityBlock = (iEvent.id()).luminosityBlock();
  
  hist->Fill(0);

  // Check if the analysis trigger fired
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*HLTtriggers_);
  std::string finalTriggerName="";  
  for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {    
    if(trigNames.triggerName(i).find("HLT_DoubleMu4_3_LowMass_v")!= std::string::npos){ 
      if(HLTtriggers_->accept(i)){
	DoubleMu_fired = true;
	finalTriggerName=trigNames.triggerName(i); 
      }
    }
  }

  // Take needed collections
  iEvent.getByToken(muonToken_ , muons_ );
  iEvent.getByToken(packedgenparticleToken_ , packedgenparticles_ );

  // Output collections
  std::vector<pat::Muon> muoncollection; 
  muoncollection.clear();
  

  // Offline muons
  for(size_t imuon = 0; imuon < muons_->size(); imuon++){
    const pat::Muon & muon = (*muons_)[imuon];

    // Offline cuts
    if(muon.pt() < 2.5) continue;        
    if(fabs(muon.eta()) > 2.6) continue;     
    if(!(muon.track().isNonnull())) continue;
    
    // This is to further process 
    muoncollection.push_back(muon);  
  }

  if (muoncollection.size() < 2) return; 
  hist->Fill(1);


  // Prepare offline muons pairs
  float jpsi_max_pt = -1;
  int mcidx_m1 = -1;
  int mcidx_m2 = -1;
  TLorentzVector jpsi_tlv_highest;

  for(int im = 0; im < (int)muoncollection.size(); im++){
    for(int jm = im+1; jm < (int)muoncollection.size(); jm++){
      
      const pat::Muon mu1 = muoncollection[im];
      const pat::Muon mu2 = muoncollection[jm];

      TLorentzVector tlv_m1;
      TLorentzVector tlv_m2;
      tlv_m1.SetPtEtaPhiM(mu1.pt(), mu1.eta(), mu1.phi(), aux.mass_muon);
      tlv_m2.SetPtEtaPhiM(mu2.pt(), mu2.eta(), mu2.phi(), aux.mass_muon);

      TLorentzVector tlv_jpsi = (tlv_m1 + tlv_m2);
      float jpsi_mass = tlv_jpsi.M();
      float jpsi_pt = tlv_jpsi.Pt();

      if (mu1.charge() + mu2.charge() !=0) continue;
      if (jpsi_mass < 2.0) continue; 
      if (jpsi_mass > 4.0) continue;

      if(jpsi_max_pt < jpsi_pt){
	jpsi_max_pt = jpsi_pt;
	mcidx_m1 = im;
	mcidx_m2 = jm;
	jpsi_tlv_highest = tlv_jpsi;
      }
    }
  }

  // At least 1 reco J/psi
  if(jpsi_max_pt == -1) return;  
  hist->Fill(2);


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



  // Match with gen muon
  std::vector<float> genMatchPt(2, -99);
  std::vector<float> genMatchEta(2, -99);
  std::vector<float> genMatchPhi(2, -99);
  std::vector<float> genMatchDR(2, -99);

  std::vector<pat::Muon> muonpair;
  muonpair.clear();
  muonpair.push_back(muoncollection[mcidx_m1]);
  muonpair.push_back(muoncollection[mcidx_m2]);

  TVector3 muTV3, genpartTV3;
  for (size_t imu =0; imu < muonpair.size(); imu++){

    float best_match_dR = 999.;
    muTV3.SetPtEtaPhi( muonpair[imu].pt(), muonpair[imu].eta(), muonpair[imu].phi() );

    for (size_t igenp = 0; igenp < packedgenparticles_->size(); ++igenp){

      pat::PackedGenParticle genp = (*packedgenparticles_)[igenp];

      if (abs(genp.pdgId()) != 13) continue;
      genpartTV3.SetPtEtaPhi( genp.pt(), genp.eta(), genp.phi() );

      Float_t deltaR = fabs(muTV3.DeltaR(genpartTV3));
      if (deltaR<0.1 && deltaR < best_match_dR ){  
	best_match_dR = deltaR;
	genMatchPt[imu]  = genpartTV3.Pt(); 
	genMatchEta[imu] = genpartTV3.Eta(); 
	genMatchPhi[imu] = genpartTV3.Phi(); 
	genMatchDR[imu]  = deltaR; 
      }
    } // Loop over gen
  } // Loop over reco muons


  // Reconstruct and store Jpsi gen pt,eta,phi,mass
  TLorentzVector tlv_gen_m1;
  TLorentzVector tlv_gen_m2;
  if ( (genMatchPt[0]>=0.) && (genMatchPt[1]>=0.) && (genMatchPt[0]!=genMatchPt[1]) ){
    tlv_gen_m1.SetPtEtaPhiM(genMatchPt[0], genMatchEta[0], genMatchPhi[0], aux.mass_muon);
    tlv_gen_m2.SetPtEtaPhiM(genMatchPt[1], genMatchEta[1], genMatchPhi[1], aux.mass_muon);

    TLorentzVector tlv_gen_jpsi = (tlv_gen_m1 + tlv_gen_m2);
    JpsiGen_pt   = tlv_gen_jpsi.Pt();
    JpsiGen_eta  = tlv_gen_jpsi.Eta();
    JpsiGen_phi  = tlv_gen_jpsi.Phi();
    JpsiGen_mass = tlv_gen_jpsi.M();
  }

  hist->Fill(3);

  
  // Loop over di-mu paths, when fired
  if (DoubleMu_fired==1) {   

    std::string DoubleMuTrigName = "HLT_DoubleMu4_3_LowMass_v";      
    std::string DoubleMuObjName  = "hltDoubleMu43LowMassL3Filtered";
    
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
	if(pathNamesAll[h].find(DoubleMuTrigName)!= std::string::npos) {
	  isPathExist = true;   
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
  Jpsi_m1_genMatchPt  = genMatchPt[0];
  Jpsi_m1_genMatchEta = genMatchEta[0];
  Jpsi_m1_genMatchPhi = genMatchPhi[0];
  Jpsi_m1_genMatchDR  = genMatchDR[0];

  Jpsi_m2_pt   = muoncollection[mcidx_m2].pt();
  Jpsi_m2_eta  = muoncollection[mcidx_m2].eta();
  Jpsi_m2_phi  = muoncollection[mcidx_m2].phi();
  Jpsi_m2_mass = muoncollection[mcidx_m2].mass();
  Jpsi_m2_q    = muoncollection[mcidx_m2].charge();
  Jpsi_m2_looseId = muoncollection[mcidx_m2].isLooseMuon();
  Jpsi_m2_bestL1pt  = bestMatchM2_pt;
  Jpsi_m2_bestL1eta = bestMatchM2_eta;
  Jpsi_m2_bestL1phi = bestMatchM2_phi;
  Jpsi_m2_bestL1dR  = bestMatchM2_dR;  
  Jpsi_m2_genMatchPt  = genMatchPt[1];
  Jpsi_m2_genMatchEta = genMatchEta[1];
  Jpsi_m2_genMatchPhi = genMatchPhi[1];
  Jpsi_m2_genMatchDR  = genMatchDR[1];

  Jpsi_nonfit_pt   = jpsi_tlv_highest.Pt();
  Jpsi_nonfit_eta  = jpsi_tlv_highest.Eta();
  Jpsi_nonfit_phi  = jpsi_tlv_highest.Phi();
  Jpsi_nonfit_mass = jpsi_tlv_highest.M(); 

  Jpsi_muonsDr = mu1TV3.DeltaR(mu2TV3); 

  tree_->Fill();

  return;

} //NanoAnalyzerDoubleMuMC::analyze ends


//**************************************************
//************* additional methods *****************
//**************************************************

void NanoAnalyzerDoubleMuMC::beginJob(const edm::ParameterSet& iConfig) { }

void NanoAnalyzerDoubleMuMC::beginRun(const edm::Run &iRun, const edm::EventSetup &iStp) { }

void NanoAnalyzerDoubleMuMC::fillDescriptions(edm::ConfigurationDescriptions & descriptions)  { }

void NanoAnalyzerDoubleMuMC::endRun(edm::Run const&, edm::EventSetup const&) { }

void NanoAnalyzerDoubleMuMC::endJob() { }

// branch title creation
void NanoAnalyzerDoubleMuMC::createBranch() { 

  tree_->Branch("run", &run, "run/i");
  tree_->Branch("event", &event, "event/l");
  tree_->Branch("luminosityBlock", &luminosityBlock, "luminosityBlock/i");
  tree_->Branch("nvtx", &nvtx, "nvtx/i");

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
  tree_->Branch("Jpsi_m1_genMatchPt",  &Jpsi_m1_genMatchPt  );
  tree_->Branch("Jpsi_m1_genMatchEta", &Jpsi_m1_genMatchEta );
  tree_->Branch("Jpsi_m1_genMatchPhi", &Jpsi_m1_genMatchPhi );
  tree_->Branch("Jpsi_m1_genMatchDR",  &Jpsi_m1_genMatchDR );

  tree_->Branch("Jpsi_m2_pt", &Jpsi_m2_pt );
  tree_->Branch("Jpsi_m2_eta", &Jpsi_m2_eta );
  tree_->Branch("Jpsi_m2_phi", &Jpsi_m2_phi );
  tree_->Branch("Jpsi_m2_mass", &Jpsi_m2_mass );
  tree_->Branch("Jpsi_m2_q", &Jpsi_m2_q );
  tree_->Branch("Jpsi_m2_looseId"   , &Jpsi_m2_looseId    );
  tree_->Branch("Jpsi_m2_bestL1pt",   &Jpsi_m2_bestL1pt   );
  tree_->Branch("Jpsi_m2_bestL1eta" , &Jpsi_m2_bestL1eta  );
  tree_->Branch("Jpsi_m2_bestL1phi" , &Jpsi_m2_bestL1phi  );
  tree_->Branch("Jpsi_m2_bestL1dR",   &Jpsi_m2_bestL1dR );
  tree_->Branch("Jpsi_m2_genMatchPt",  &Jpsi_m2_genMatchPt );
  tree_->Branch("Jpsi_m2_genMatchEta", &Jpsi_m2_genMatchEta );
  tree_->Branch("Jpsi_m2_genMatchPhi", &Jpsi_m2_genMatchPhi );
  tree_->Branch("Jpsi_m2_genMatchDR",  &Jpsi_m2_genMatchDR );

  tree_->Branch("Jpsi_m1_Dimu_dR",     &mu1_matchedDiMu_dR );
  tree_->Branch("Jpsi_m1_Dimu_pt",     &mu1_matchedDiMu_pt );
  tree_->Branch("Jpsi_m1_Dimu_eta",    &mu1_matchedDiMu_eta );
  tree_->Branch("Jpsi_m1_Dimu_phi",    &mu1_matchedDiMu_phi );
  tree_->Branch("Jpsi_m2_Dimu_dR",     &mu2_matchedDiMu_dR );
  tree_->Branch("Jpsi_m2_Dimu_pt",     &mu2_matchedDiMu_pt );
  tree_->Branch("Jpsi_m2_Dimu_eta",    &mu2_matchedDiMu_eta );
  tree_->Branch("Jpsi_m2_Dimu_phi",    &mu2_matchedDiMu_phi );

  tree_->Branch("DoubleMu_fired", &DoubleMu_fired );

  tree_->Branch("L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4", &L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4 );
  tree_->Branch("L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6", &L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6 );
  tree_->Branch("L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6", &L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6 );
  tree_->Branch("L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5", &L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5 );
  tree_->Branch("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 );
  tree_->Branch("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 );
  tree_->Branch("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2 );
  tree_->Branch("L1_DoubleMu4_SQ_OS_dR_Max1p2", &L1_DoubleMu4_SQ_OS_dR_Max1p2 );
  tree_->Branch("L1_OR", &L1_OR );

  tree_->Branch("JpsiGen_pt"  , &JpsiGen_pt   );
  tree_->Branch("JpsiGen_eta" , &JpsiGen_eta  );
  tree_->Branch("JpsiGen_phi" , &JpsiGen_phi  );
  tree_->Branch("JpsiGen_mass", &JpsiGen_mass );

  tree_->Branch("Jpsi_nonfit_pt",   &Jpsi_nonfit_pt );
  tree_->Branch("Jpsi_nonfit_eta",  &Jpsi_nonfit_eta );
  tree_->Branch("Jpsi_nonfit_phi",  &Jpsi_nonfit_phi);
  tree_->Branch("Jpsi_nonfit_mass", &Jpsi_nonfit_mass );
  tree_->Branch("Jpsi_muonsDr",     &Jpsi_muonsDr );
}

void NanoAnalyzerDoubleMuMC::reset(void){

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
  Jpsi_m1_genMatchPt = -99;
  Jpsi_m1_genMatchEta = -99;
  Jpsi_m1_genMatchPhi = -99;
  Jpsi_m1_genMatchDR = -99;

  Jpsi_m2_pt = -99;
  Jpsi_m2_eta = -99;
  Jpsi_m2_phi = -99;
  Jpsi_m2_mass = -99;
  Jpsi_m2_q = -99;
  Jpsi_m2_looseId = -99;
  Jpsi_m2_bestL1pt = -99;
  Jpsi_m2_bestL1eta = -99;
  Jpsi_m2_bestL1phi = -99;
  Jpsi_m2_bestL1dR = -99;
  Jpsi_m2_genMatchPt = -99;
  Jpsi_m2_genMatchEta = -99;
  Jpsi_m2_genMatchPhi = -99;
  Jpsi_m2_genMatchDR = -99;

  JpsiGen_pt   = -99;
  JpsiGen_eta  = -99;
  JpsiGen_phi  = -99;
  JpsiGen_mass = -99;

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

  Jpsi_nonfit_pt   = -99;
  Jpsi_nonfit_eta  = -99;
  Jpsi_nonfit_phi  = -99;
  Jpsi_nonfit_mass = -99;
  Jpsi_muonsDr     = -99;

  // Not for tree
  dimuobj1_pt  = -999.;
  dimuobj1_eta = -999.;
  dimuobj1_phi = -999.;
  dimuobj2_pt  = -999.;
  dimuobj2_eta = -999.;
  dimuobj2_phi = -999.;
}

DEFINE_FWK_MODULE(NanoAnalyzerDoubleMuMC);
