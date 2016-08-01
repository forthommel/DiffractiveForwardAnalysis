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
//         Created:  Sun, 22 May 2016 19:05:16 GMT
//
//

#include "GammaGammaGammaGamma.h"

//
// constructors and destructor
//
GammaGammaGammaGamma::GammaGammaGammaGamma(const edm::ParameterSet& iConfig) :
  fPhotonToken(consumes< edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photonTag"))),
  fConversionToken(consumes< edm::View<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("conversionTag"))),
  fPFCandidateToken(consumes< edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("pfCandidateTag"))),
  fVertexToken(consumes< edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertexTag"))),
  //fBeamSpotToken(consumes< edm::View<reco::BeamSpot> >(iConfig.getParameter<edm::InputTag>("beamspotTag"))),
  fPhotonMediumIdBoolMapToken(consumes< edm::ValueMap<bool> >(
    iConfig.getUntrackedParameter<edm::InputTag>("photonMedIdBoolMapTag", edm::InputTag("egmPhotonIDs:mvaPhoID-PHYS14-PU20bx25-nonTrig-V1-wp90"))
  )),
  fPhotonMediumIdFullInfoMapToken(consumes< edm::ValueMap< vid::CutFlowResult > >(
    iConfig.getUntrackedParameter<edm::InputTag>("photonMedIdInfoMapTag", edm::InputTag("egmPhotonIDs:mvaPhoID-PHYS14-PU20bx25-nonTrig-V1-wp90"))
  )),
  fTree(0)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  //TFileDirectory dir = fs->mkdir("DiPhotonAnalyzer");
  TFileDirectory dir = fs->mkdir(".");

  fTree = fs->make<TTree>("ntp1", "diphoton candidates tree");
}


GammaGammaGammaGamma::~GammaGammaGammaGamma()
{}


// ------------ method called for each event  ------------
void
GammaGammaGammaGamma::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle< edm::View<pat::Photon> > photons;
  iEvent.getByToken(fPhotonToken, photons);

  edm::Handle< edm::View<reco::Conversion> > conversions;
  iEvent.getByToken(fConversionToken, conversions);

  edm::Handle< edm::View<pat::PackedCandidate> > pfcandidates;
  iEvent.getByToken(fPFCandidateToken, pfcandidates);

  edm::Handle< edm::View<reco::Vertex> > vertices;
  iEvent.getByToken(fVertexToken, vertices);

  // The first map simply has pass/fail for each particle
  iEvent.getByToken(fPhotonMediumIdBoolMapToken, fPhotonMediumIdDecisions);
  // The second map has the full info about the cut flow
  iEvent.getByToken(fPhotonMediumIdFullInfoMapToken, fPhotonMediumIdCutflowData);

  // Get MVA values and categories (optional)
  /*edm::Handle<edm::ValueMap<float> > mvaValues;
  iEvent.getByToken(mvaValuesMapToken_,mvaValues);
  edm::Handle<edm::ValueMap<int> > mvaCategories;
  iEvent.getByToken(mvaCategoriesMapToken_,mvaCategories);*/

  //fVertexFinder.SetEventHandles(conversions, pfcandidates, vertices);

  TLorentzVector ph1, ph2;
 
  aRunId = iEvent.id().run();
  aLSId = iEvent.id().luminosityBlock();
  aEventId = iEvent.id().event();

  unsigned int i=0, j=0;
  for (std::vector< edm::Ptr<pat::Photon> >::const_iterator photon1_ptr=photons->ptrs().begin(); photon1_ptr!=photons->ptrs().end(); photon1_ptr++) {
    const edm::Ptr<pat::Photon> photon1 = *photon1_ptr;
    if (!passPhotonId(photon1)) continue;
    
    aPhotonPt[i] = photon1->pt();
    aPhotonEta[i] = photon1->eta();
    aPhotonPhi[i] = photon1->phi();
    aPhotonVtxX[i] = photon1->vx();
    aPhotonVtxY[i] = photon1->vy();
    aPhotonVtxZ[i] = photon1->vz();

    ph1.SetPxPyPzE(photon1->p4().px(), photon1->p4().py(), photon1->p4().pz(), photon1->p4().E());

    for (std::vector< edm::Ptr<pat::Photon> >::const_iterator photon2_ptr=photon1_ptr+1; photon2_ptr!=photons->ptrs().end(); photon2_ptr++) {
      const edm::Ptr<pat::Photon> photon2 = *photon2_ptr;
      if (!passPhotonId(photon2)) continue;

      //reco::Vertex vtx_matched = fVertexFinder.FindVertex(*photon, *photon2);

      ph2.SetPxPyPzE(photon2->p4().px(), photon2->p4().py(), photon2->p4().pz(), photon2->p4().E());
      const TLorentzVector photon_pair = ph1+ph2;
      double dphi = fabs(ph1.Phi()-ph2.Phi());

      //aPhotonPairVertexDist[j] = std::sqrt((photon->vertex()-photon2->vertex()).mag2());
      aPhotonPairVertexDist[j] = std::sqrt(pow(photon1->vx()-photon2->vx(), 2)+
                                           pow(photon1->vy()-photon2->vy(), 2)+
                                           pow(photon1->vz()-photon2->vz(), 2));
      aPhotonPairDphi[j] = (dphi<TMath::Pi()) ? dphi : 2.*TMath::Pi()-dphi; // dphi lies in [-pi, pi]
      aPhotonPairDpt[j] = fabs(ph1.Pt()-ph2.Pt());
      aPhotonPairMass[j] = photon_pair.M();
      aPhotonPairPt[j] = photon_pair.Pt();
      j++;
    }

    i++;
  }
  aNumPhotons = i;
  aNumPhotonPairs = j;

  /*for (edm::View<reco::Vertex>::const_iterator vtx=vertices->begin(); vtx!=vertices->end(); vtx++) {
    if (passVertexCriteria(*vtx, *photon1, *photon2)) matched_vertex = true;
  }*/

  fTree->Fill();

}

bool
GammaGammaGammaGamma::passPhotonId(const edm::Ptr< pat::Photon >& photon_ref) const
{
  const bool pass_medium = (*fPhotonMediumIdDecisions)[photon_ref];

  return pass_medium;
}


// ------------ method called once each job just before starting event loop  ------------
void 
GammaGammaGammaGamma::beginJob()
{
  fTree->Branch("run_id", &aRunId, "run_id/I");
  fTree->Branch("lumi_section", &aLSId, "lumi_section/I");
  fTree->Branch("event_id", &aEventId, "event_id/I");
  fTree->Branch("num_photons", &aNumPhotons, "num_photons/I");
  fTree->Branch("photon_pt", aPhotonPt, "photon_pt[num_photons]/D");
  fTree->Branch("photon_eta", aPhotonEta, "photon_eta[num_photons]/D");
  fTree->Branch("photon_phi", aPhotonPhi, "photon_phi[num_photons]/D");
  fTree->Branch("photon_vtxx", aPhotonVtxX, "photon_vtxx[num_photons]/D");
  fTree->Branch("photon_vtxy", aPhotonVtxY, "photon_vtxy[num_photons]/D");
  fTree->Branch("photon_vtxz", aPhotonVtxZ, "photon_vtxz[num_photons]/D");
  fTree->Branch("num_photonpairs", &aNumPhotonPairs, "num_photonpairs/I");
  fTree->Branch("photonpair_vertex_dist", aPhotonPairVertexDist, "photonpair_vertex_dist[num_photonpairs]/D");
  fTree->Branch("photonpair_mass", aPhotonPairMass, "photonpair_mass[num_photonpairs]/D");
  fTree->Branch("photonpair_pt", aPhotonPairPt, "photonpair_pt[num_photonpairs]/D");
  fTree->Branch("photonpair_dphi", aPhotonPairDphi, "photonpair_dphi[num_photonpairs]/D");
  fTree->Branch("photonpair_dpt", aPhotonPairDpt, "photonpair_dpt[num_photonpairs]/D");
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GammaGammaGammaGamma::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GammaGammaGammaGamma::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
