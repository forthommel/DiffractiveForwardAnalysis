// -*- C++ -*-
//
// Package:    DiffractiveForwardAnalysis/GammaGammaLeptonLepton
// Class:      GammaGammaGammaGamma
// 
/**\class GammaGammaGammaGamma GammaGammaGammaGamma.cc DiffractiveForwardAnalysis/GammaGammaLeptonLepton/plugins/GammaGammaGammaGamma.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Forthomme
//         Created:  Sun, 29 May 2016 19:05:16 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/VertexFinder.h"
//
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "TTree.h"
#include "TLorentzVector.h"

#define MAX_PHOTONS 100
#define MAX_PHOTON_PAIRS 100

class GammaGammaGammaGamma : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit GammaGammaGammaGamma(const edm::ParameterSet&);
      ~GammaGammaGammaGamma();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      bool passPhotonId(const edm::Ptr< pat::Photon >& photon_ref) const;

      edm::EDGetTokenT< edm::View<pat::Photon> > fPhotonToken;
      edm::EDGetTokenT< edm::View<reco::Conversion> > fConversionToken;
      edm::EDGetTokenT< edm::View<pat::PackedCandidate> > fPFCandidateToken;
      edm::EDGetTokenT< edm::View<reco::Vertex> > fVertexToken;
      edm::EDGetTokenT< edm::ValueMap<bool> > fPhotonMediumIdBoolMapToken;
      edm::EDGetTokenT< edm::ValueMap<vid::CutFlowResult> > fPhotonMediumIdFullInfoMapToken;
      
      // The first map simply has pass/fail for each particle
      edm::Handle< edm::ValueMap<bool> > fPhotonMediumIdDecisions;
      // The second map has the full info about the cut flow
      edm::Handle< edm::ValueMap<vid::CutFlowResult> > fPhotonMediumIdCutflowData;

      //edm::EDGetTokenT< reco::BeamSpot > fBeamSpotToken;
      
      TTree* fTree;

      //VertexFinder fVertexFinder;
      //
      int aRunId;
      int aLSId;
      int aEventId;

      int aNumPhotons;
      double aPhotonPt[MAX_PHOTONS];
      double aPhotonEta[MAX_PHOTONS];
      double aPhotonPhi[MAX_PHOTONS];
      double aPhotonVtxX[MAX_PHOTONS];
      double aPhotonVtxY[MAX_PHOTONS];
      double aPhotonVtxZ[MAX_PHOTONS];

      int aNumPhotonPairs;
      double aPhotonPairVertexDist[MAX_PHOTON_PAIRS];
      double aPhotonPairMass[MAX_PHOTON_PAIRS];
      double aPhotonPairPt[MAX_PHOTON_PAIRS];
      double aPhotonPairDpt[MAX_PHOTON_PAIRS];
      double aPhotonPairDphi[MAX_PHOTON_PAIRS];
};

//define this as a plug-in
DEFINE_FWK_MODULE(GammaGammaGammaGamma);
