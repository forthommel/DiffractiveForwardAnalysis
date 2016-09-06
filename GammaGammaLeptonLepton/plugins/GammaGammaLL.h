#ifndef DiffractiveForwardAnalysis_GammaGammaLL_h
#define DiffractiveForwardAnalysis_GammaGammaLL_h

// system include files
#include <fstream>
#include <memory>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/Registry.h"

// HLT information
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

// Generator level collection
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

// Pileup
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "DiffractiveForwardAnalysis/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// Muons collection
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/MuonReco/interface/MuonFwd.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"

// Electrons collection
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"
//#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

// Photons collection
#include "DataFormats/PatCandidates/interface/Photon.h"

// Particle flow collection
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// Vertices collection
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h" 
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

// Jets/MET collection
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

// CT-PPS objects
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

// HPS acceptance
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/AcceptanceTableHelper.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/PrimaryVertexSelector.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/HLTMatcher.h"

// LHC fill information
//#include "DataFormats/Common/interface/ConditionsInEdm.h" // L1 method
//#include "CondFormats/RunInfo/interface/FillInfo.h"
//#include "CondFormats/DataRecord/interface/FillInfoRcd.h" // db method

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>

#define MAX_HLT    10   // Maximum number of HLT to check
/*#define MAX_LL     50   // Maximum number of leptons per event
#define MAX_MUONS  25   // Maximum number of muons per event
#define MAX_ELE    25   // Maximum number of electrons per event
#define MAX_PHO    50   // Maximum number of photons per event
#define MAX_PAIRS  25   // Maximum number of leptons pairs per event
#define MAX_VTX    1000 // Maximum number of primary vertices per event
#define MAX_ET     10000// Maximum number of extra tracks per event
#define MAX_GENMU  25   // Maximum number of generator level muons per event
#define MAX_GENELE 25   // Maximum number of generator level electrons per event
#define MAX_GENPHO 10   // Maximum number of generator level photons per event
#define MAX_GENPRO 8    // Maximum number of generator level protons per event
#define MAX_JETS   30   // Maximum number of jets per event
#define MAX_LOCALPCAND 10 // Maximum number of reconstructed local tracks in RPs
#define MAX_LOCALPPAIRCAND 5 // Maximum number of reconstructed local tracks pairs in RPs*/

#define MASS_MU 0.1057
#define MASS_E  0.000511
#define MASS_P  0.938272029
#define pi 3.14159265359

typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals; 

//
// class declaration
//

class GammaGammaLL : public edm::EDAnalyzer {
   public:
      explicit GammaGammaLL(const edm::ParameterSet&);
      ~GammaGammaLL();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      
      virtual void lookAtTriggers(const edm::Event&, const edm::EventSetup&);
      virtual void analyzeMCEventContent(const edm::Event&);
      virtual void extractElectrons(const edm::Event&, std::map<int, TLorentzVector>*);
      virtual void extractMuons(const edm::Event&, std::map<int, TLorentzVector>*);
      void clearTree();

      // ----------member data ---------------------------

      bool fetchMuons_, fetchElectrons_, fetchProtons_;
      
      unsigned int verb_;

      std::ofstream *logfile_;
      
      // Input tags
      std::string outputFile_;
      std::vector<std::string> leptonsType_;
      std::string hltMenuLabel_;
      std::vector<std::string> triggersList_;

      edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
      edm::EDGetTokenT<reco::VertexCollection> recoVertexToken_;
      edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
      edm::EDGetTokenT< edm::View<pat::Muon> > muonToken_;
      edm::EDGetTokenT< edm::View<pat::Electron> > eleToken_;
      edm::EDGetTokenT< edm::ValueMap<bool> > eleLooseIdMapToken_;
      edm::EDGetTokenT< edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT< edm::ValueMap<bool> > eleTightIdMapToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
      edm::EDGetTokenT< std::vector<PileupSummaryInfo> > pileupToken_;
      edm::EDGetTokenT< edm::View<pat::Jet> > jetToken_;
      edm::EDGetTokenT< edm::View<pat::MET> > metToken_;
      edm::EDGetTokenT< edm::DetSetVector<TotemRPLocalTrack> > totemRPHitToken_;
      edm::EDGetTokenT< edm::View<pat::Photon> > photonToken_;
      edm::EDGetTokenT< edm::View<reco::PFCandidate> > pflowToken_;

      std::vector<edm::InputTag> isoValLabel_;

      bool runOnMC_, printCandidates_;
      double minPtMC_, minEtaMC_;
      double sqrts_;
      unsigned int maxExTrkVtx_;

      // Trigger information
      HLTMatcher* hlts_;
      HLTConfigProvider hltConfig_;
      HLTPrescaleProvider hltPrescale_;

      // Pileup information
      edm::LumiReWeighting *lumiWeights_;
      std::string mcPileupFile_, mcPileupPath_, dataPileupFile_, dataPileupPath_;
      
      // Isolation
      double rhoIso;
      double iso_ch, iso_em, iso_nh; // Electron isolation quantities
      int vtxind; // Primary vertex index (used in loop over vertices)
      int etind; // Extra tracks on vertex index (used in loop over tracks)

      ////// Tree contents //////
      
      // Run/event quantities
      int BX, Run, LumiSection, EventNum;
      //int LHCFillNum, LHCBeamMode;
      //double AvgInstDelLumi, BunchInstLumi[3]; 
      
      // HLT quantities
      int nHLT;
      int HLT_Accept[MAX_HLT], HLT_Prescl[MAX_HLT];
      //char* HLT_Name[MAX_HLT];
      std::vector<std::string>* HLT_Name;
      
      // Generator level quantities
      int nGenMuonCandOutOfAccept, nGenEleCandOutOfAccept, nGenPhotCandOutOfAccept;
      vector<double> GenMuonCand_px, GenMuonCand_py, GenMuonCand_pz, GenMuonCand_p, GenMuonCand_pt, GenMuonCand_eta, GenMuonCand_phi;
      vector<double> GenEleCand_px, GenEleCand_py, GenEleCand_pz, GenEleCand_p, GenEleCand_pt, GenEleCand_eta, GenEleCand_phi;
      double GenPair_p, GenPair_pt, GenPair_mass, GenPair_phi, GenPair_eta, GenPair_dphi, GenPair_dpt, GenPair_3Dangle;
      vector<double> GenPhotCand_p, GenPhotCand_e, GenPhotCand_pt, GenPhotCand_eta, GenPhotCand_phi;
      vector<double> GenProtCand_p, GenProtCand_px, GenProtCand_py, GenProtCand_pz, GenProtCand_pt, GenProtCand_eta, GenProtCand_phi;
      vector<int> GenProtCand_status;

      // HPS quantities
      double xi, t;

      int nLeptonCand, nLeptonsInPrimVertex, nCandidates, nCandidatesInEvent;

      // Pileup reweighting quantities
      double nTruePUafterPUWeight;
      double nTruePUafterPUWeightBXM1, nTruePUafterPUWeightBXP1, nTruePUafterPUWeightBX0;
      double PUWeightTrue;
      int nTruePUforPUWeight;
      int nTruePUforPUWeightBXM1, nTruePUforPUWeightBXP1, nTruePUforPUWeightBX0;
      double Weight;

      // Muon quantities
      vector<double> MuonCand_px, MuonCand_py, MuonCand_pz, MuonCand_p, MuonCand_pt, MuonCand_eta, MuonCand_phi;
      vector<double> MuonCand_vtxx, MuonCand_vtxy, MuonCand_vtxz;
      vector<int> MuonCand_charge;
      vector<double> MuonCand_dxy, MuonCand_dz;
      vector<int> MuonCand_nstatseg, MuonCand_npxlhits, MuonCand_ntrklayers;
      vector<int> MuonCandTrack_nmuchits;
      vector<double> MuonCandTrack_chisq;
      vector<int> MuonCand_isglobal, MuonCand_istracker, MuonCand_isstandalone, MuonCand_ispfmuon, MuonCand_istight;

      // Electron quantities
      int nEleCand;
      vector<double> EleCand_px, EleCand_py, EleCand_pz, EleCand_p, EleCand_e, EleCand_et, EleCand_eta, EleCand_phi;
      vector<double> EleCand_vtxx, EleCand_vtxy, EleCand_vtxz;
      vector<int> EleCand_charge;
      vector<double> EleCandTrack_p, EleCandTrack_pt, EleCandTrack_eta, EleCandTrack_phi;
      vector<double> EleCandTrack_vtxz; 
      vector<double> EleCand_deltaPhi, EleCand_deltaEta;
      vector<double> EleCand_HoverE;
      vector<double> EleCand_trackiso, EleCand_ecaliso, EleCand_hcaliso, EleCand_sigmaIetaIeta;
      vector<double> EleCand_convDist, EleCand_convDcot;
      vector<int> EleCand_ecalDriven;
      vector<int> EleCand_tightID, EleCand_mediumID, EleCand_looseID;
      
      // Photon quantities
      int nPhotonCand;
      double PhotonCand_px[MAX_PHO], PhotonCand_py[MAX_PHO], PhotonCand_pz[MAX_PHO];
      double PhotonCand_p[MAX_PHO], PhotonCand_pt[MAX_PHO];
      double PhotonCand_eta[MAX_PHO], PhotonCand_phi[MAX_PHO], PhotonCand_r9[MAX_PHO];
      double PhotonCand_drtrue[MAX_PHO], PhotonCand_detatrue[MAX_PHO], PhotonCand_dphitrue[MAX_PHO];
      
      // Pair quantities
      vector< pair<int, int> > Pair_candidates;
      vector<int> Pair_lepton1, Pair_lepton2;
      vector<double> Pair_mindist;
      vector<double> Pair_p, Pair_pt, Pair_dpt, Pair_mass, Pair_dphi, Pair_eta, Pair_phi, Pair_3Dangle;
      double PairGamma_mass[MAX_PAIRS][MAX_PHO];
      // Extra tracks
      vector<int> Pair_extratracks1mm, Pair_extratracks2mm, Pair_extratracks3mm,
                  Pair_extratracks4mm, Pair_extratracks5mm, Pair_extratracks1cm,
                  Pair_extratracks2cm, Pair_extratracks3cm, Pair_extratracks4cm,
                  Pair_extratracks5cm, Pair_extratracks10cm;
      
      // Vertex quantities
      vector<int> PrimVertexCand_id, PrimVertexCand_hasdil;
      vector<double> PrimVertexCand_x, PrimVertexCand_y, PrimVertexCand_z;
      vector<int> PrimVertexCand_tracks, PrimVertexCand_matchedtracks, PrimVertexCand_unmatchedtracks;
      vector<double> PrimVertexCand_chi2;
      vector<int> PrimVertexCand_ndof;
      
      // Extra tracks on vertex quantities
      vector<int> ExtraTrack_purity, ExtraTrack_nhits, ExtraTrack_charge, ExtraTrack_ndof, ExtraTrack_vtxId;
      vector<double> ExtraTrack_p, ExtraTrack_pt, ExtraTrack_px, ExtraTrack_py, ExtraTrack_pz, ExtraTrack_eta, ExtraTrack_phi;
      vector<double> ExtraTrack_chi2;
      vector<double> ExtraTrack_vtxdxyz, ExtraTrack_vtxT, ExtraTrack_vtxZ;
      vector<double> ExtraTrack_x, ExtraTrack_y, ExtraTrack_z;
      vector<double> ClosestExtraTrack_vtxdxyz, ClosestHighPurityExtraTrack_vtxdxyz;
      vector<int> ClosestExtraTrack_id, ClosestHighPurityExtraTrack_id;
      int nQualityExtraTrack;

      // Jets/MET quantities
      vector<double> JetCand_px, JetCand_py, JetCand_pz, JetCand_e, JetCand_eta, JetCand_phi;
      double HighestJet_e, HighestJet_eta, HighestJet_phi;
      double SumJet_e;
      double Etmiss, Etmiss_phi, Etmiss_x, Etmiss_y, Etmiss_z, Etmiss_significance;

      // CTPPS quantities
      vector<double> LocalProtCand_x, LocalProtCand_y, LocalProtCand_z, LocalProtCand_xSigma, LocalProtCand_ySigma;
      vector<double> LocalProtCand_Tx, LocalProtCand_Ty, LocalProtCand_TxSigma, LocalProtCand_TySigma;
      vector<int> LocalProtCand_arm, LocalProtCand_side;

      int nLocalProtPairCand;
      vector<double> LocalProtPairCand_mass, LocalProtPairCand_pt, LocalProtPairCand_y;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
TFile* file_;
TTree* tree_;

#endif
