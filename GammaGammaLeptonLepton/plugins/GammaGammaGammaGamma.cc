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
  fPhotonToken                      (consumes< edm::View<pat::Photon> >               (iConfig.getParameter<edm::InputTag>("photonTag"))),
  fConversionToken                  (consumes< std::vector<reco::Conversion> >        (iConfig.getParameter<edm::InputTag>("conversionTag"))),
  fElectronToken                    (consumes< edm::View<pat::Electron> >             (iConfig.getParameter<edm::InputTag>("electronTag"))),
  fRecoElectronToken                (consumes< std::vector<reco::GsfElectron> >       (iConfig.getParameter<edm::InputTag>("recoElectronTag"))),
  fPFCandidateToken                 (consumes< edm::View<pat::PackedCandidate> >      (iConfig.getParameter<edm::InputTag>("pfCandidateTag"))),
  fVertexToken                      (consumes< edm::View<reco::Vertex> >              (iConfig.getParameter<edm::InputTag>("vertexTag"))),
  fBeamSpotToken                    (consumes< reco::BeamSpot >                       (iConfig.getParameter<edm::InputTag>("beamspotTag"))),
  fPhotonMediumIdBoolMapToken       (consumes< edm::ValueMap<bool> >                  (iConfig.getParameter<edm::InputTag>("photonMedIdBoolMapTag"))),
  //fPhotonMediumIdFullInfoMapToken   (consumes< edm::ValueMap< vid::CutFlowResult > >  (iConfig.getParameter<edm::InputTag>("photonMedIdInfoMapTag"))),
  fPhotonMVAIdValuesMapToken        (consumes< edm::ValueMap<float> >                 (iConfig.getParameter<edm::InputTag>("photonMVAIdValuesMapTag"))),
  fPhotonMVAIdCategoriesMapToken    (consumes< edm::ValueMap<int> >                   (iConfig.getParameter<edm::InputTag>("photonMVAIdCategoriesMapTag"))),
  fPhotonMinPt    (iConfig.getParameter<double>("photonMinPt")),
  fFetchProtons   (iConfig.getParameter<bool>  ("fetchProtons")),
  fTree(0)
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  //TFileDirectory dir = fs->mkdir("DiPhotonAnalyzer");
  TFileDirectory dir = fs->mkdir(".");

  fTree = fs->make<TTree>("ntp1", "diphoton candidates tree");

  if (fFetchProtons) {
    fProtonToken = consumes< edm::DetSetVector<TotemRPLocalTrack> >(iConfig.getParameter<edm::InputTag>("protonTag"));
  }
}


GammaGammaGammaGamma::~GammaGammaGammaGamma()
{}

void
GammaGammaGammaGamma::clearTree()
{
  aNumPhotons = 0;
  for (unsigned int i=0; i<MAX_PHOTONS; i++) {
    aPhotonPt[i] = aPhotonEta[i] = aPhotonPhi[i] = 0.;
    aPhotonVtxX[i] = aPhotonVtxY[i] = aPhotonVtxZ[i] = 0.;
    aPhotonIsoTrack[i] = aPhotonIsoECAL[i] = aPhotonIsoHCAL[i] = aPhotonIsoCalo[i] = 0.;
    aPhotonR9[i] = 0.;
    aPhotonElectronVeto[i] = aPhotonConverted[i] = 0;
  }
  aNumPhotonPairs = 0;
  for (unsigned int i=0; i<MAX_PHOTON_PAIRS; i++) {
    aPhotonPairVertexDist[i] = aPhotonPairMass[i] = aPhotonPairPt[i] = aPhotonPairDpt[i] = aPhotonPairDphi[i] = 0.;
  }
  aNumProtons = 0;
  for (unsigned int i=0; i<MAX_PROTONS; i++) {
    aProtonX[i] = aProtonZ[i] = aProtonY[i] = 0.;
    aProtonXsigma[i] = aProtonYsigma[i] = 0.;
    aProtonArm[i] = aProtonSide[i] = 0;
  }
}

// ------------ method called for each event  ------------
void
GammaGammaGammaGamma::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  clearTree();

  // The first map simply has pass/fail for each particle
  iEvent.getByToken(fPhotonMediumIdBoolMapToken, fPhotonMediumIdDecisions);
  // The second map has the full info about the cut flow
  //iEvent.getByToken(fPhotonMediumIdFullInfoMapToken, fPhotonMediumIdCutflowData);
  // Get MVA values and categories (optional)
  edm::Handle<edm::ValueMap<float> > photonIdMVAvalues;
  iEvent.getByToken(fPhotonMVAIdValuesMapToken, photonIdMVAvalues);
  edm::Handle<edm::ValueMap<int> > photonIdMVAcategories;
  iEvent.getByToken(fPhotonMVAIdCategoriesMapToken, photonIdMVAcategories);

  /*edm::Handle< edm::View<pat::PackedCandidate> > pfcandidates;
  iEvent.getByToken(fPFCandidateToken, pfcandidates);

  fVertexFinder.SetEventHandles(conversions, pfcandidates, vertices);*/

  TLorentzVector ph1, ph2;
 
  aRunId = iEvent.id().run();
  aLumiSection = iEvent.id().luminosityBlock();
  aEventNum = iEvent.id().event();

  edm::Handle< edm::View<pat::Photon> > photons;
  iEvent.getByToken(fPhotonToken, photons); 

  edm::Handle< std::vector<reco::Conversion> > conversions;
  iEvent.getByToken(fConversionToken, conversions);

  edm::Handle< std::vector<reco::GsfElectron> > reco_electrons;
  iEvent.getByToken(fRecoElectronToken, reco_electrons);
  /*edm::Handle< edm::View<pat::Electron> > electrons;
  iEvent.getByToken(fElectronToken, electrons);*/

  edm::Handle< reco::BeamSpot > beamspots;
  iEvent.getByToken(fBeamSpotToken, beamspots);
  const reco::BeamSpot& beamspot = *beamspots.product();

  unsigned int i=0, j=0;
  for (std::vector< edm::Ptr<pat::Photon> >::const_iterator photon1_ptr=photons->ptrs().begin(); photon1_ptr!=photons->ptrs().end(); photon1_ptr++) {
    const edm::Ptr<pat::Photon> photon1 = *photon1_ptr;
std::cerr << "photon 1 retrieved: pt=" << photon1->pt() << std::endl;
    if (photon1->pt()<fPhotonMinPt) continue;
    if (!passPhotonId(photon1)) continue;
    
    aPhotonPt[i] = photon1->pt();
    aPhotonEta[i] = photon1->eta();
    aPhotonPhi[i] = photon1->phi();
    aPhotonVtxX[i] = photon1->vx();
    aPhotonVtxY[i] = photon1->vy();
    aPhotonVtxZ[i] = photon1->vz();
    aPhotonIsoTrack[i] = photon1->trackIso();
    aPhotonIsoECAL[i] = photon1->ecalIso();
    aPhotonIsoHCAL[i] = photon1->hcalIso();
    aPhotonIsoCalo[i] = photon1->caloIso();
    aPhotonR9[i] = photon1->r9();

    const float lxyMin = 2.0,
                probMin = 1e-6;
    const int nHitsBeforeVtxMax = 0;

    aPhotonElectronVeto[i] = !ConversionTools::hasMatchedPromptElectron(photon1->superCluster(), reco_electrons, conversions, beamspot.position(), lxyMin, probMin, nHitsBeforeVtxMax);
    aPhotonConverted[i] = ConversionTools::hasMatchedConversion(*photon1->superCluster(), conversions, beamspot.position(), lxyMin, probMin, 1);
    if (aPhotonConverted[i]) {
      reco::ConversionRef conv = ConversionTools::matchedConversion(*photon1->superCluster(), conversions, beamspot.position());
std::cerr << "conversion found for photon 1 --> number of tracks: " <<  conv->nTracks() << std::endl;
      if (conv->nTracks()<2) continue;
      aPhotonVtxX[i] = conv->conversionVertex().x();
      aPhotonVtxY[i] = conv->conversionVertex().y();
      aPhotonVtxZ[i] = conv->conversionVertex().z();
    }

    ph1.SetPxPyPzE(photon1->p4().px(), photon1->p4().py(), photon1->p4().pz(), photon1->p4().E());

    for (std::vector< edm::Ptr<pat::Photon> >::const_iterator photon2_ptr=photon1_ptr+1; photon2_ptr!=photons->ptrs().end(); photon2_ptr++) {
      const edm::Ptr<pat::Photon> photon2 = *photon2_ptr;
std::cerr << "photon 2 retrieved: pt=" << photon2->pt() << std::endl;
      if (photon2->pt()<fPhotonMinPt) continue;
      if (!passPhotonId(photon2)) continue;

      /*aPhotonIdMVAvalue[j] = (*photonIdMVAvalues)[photon2];
      aPhotonIdMVAcategory[j] = (*photonIdMVAcategories)[photon2];*/

      //reco::Vertex vtx_matched = fVertexFinder.FindVertex(*photon, *photon2);

      ph2.SetPxPyPzE(photon2->p4().px(), photon2->p4().py(), photon2->p4().pz(), photon2->p4().E());
      const TLorentzVector photon_pair = ph1+ph2;
      double dphi = fabs(ph1.Phi()-ph2.Phi());
      // ensure dphi lies in [-pi, pi]
      while (dphi> TMath::Pi()) dphi -= 2.*TMath::Pi();
      while (dphi<-TMath::Pi()) dphi += 2.*TMath::Pi();

      double vx = photon2->vx(), vy = photon2->vy(), vz = photon2->vz();
      if (ConversionTools::hasMatchedConversion(*photon2->superCluster(), conversions, beamspot.position(), lxyMin, probMin, 1)) {
        reco::ConversionRef conv2 = ConversionTools::matchedConversion(*photon2->superCluster(), conversions, beamspot.position());
        if (conv2->nTracks()>=2) {
          vx = conv2->conversionVertex().x();
          vy = conv2->conversionVertex().y();
          vz = conv2->conversionVertex().z();
        }
      }

      //aPhotonPairVertexDist[j] = std::sqrt((photon->vertex()-photon2->vertex()).mag2());
      aPhotonPairVertexDist[j] = std::sqrt(pow(aPhotonVtxX[i]-vx, 2)+
                                           pow(aPhotonVtxY[i]-vy, 2)+
                                           pow(aPhotonVtxZ[i]-vz, 2));
      aPhotonPairDphi[j] = dphi;
      aPhotonPairDpt[j] = fabs(ph1.Pt()-ph2.Pt());
      aPhotonPairMass[j] = photon_pair.M();
      aPhotonPairPt[j] = photon_pair.Pt();
std::cerr << "photon pair with mass=" << photon_pair.M() << std::endl;
      j++;
    }

    i++;
  }
  aNumPhotons = i;
  aNumPhotonPairs = j;

  /*edm::Handle< edm::View<reco::Vertex> > vertices;
  iEvent.getByToken(fVertexToken, vertices);

  for (edm::View<reco::Vertex>::const_iterator vtx=vertices->begin(); vtx!=vertices->end(); vtx++) {
    if (passVertexCriteria(*vtx, *photon1, *photon2)) matched_vertex = true;
  }*/

  if (fFetchProtons) {
    edm::Handle< edm::DetSetVector<TotemRPLocalTrack> > rp_coll;
    iEvent.getByToken(fProtonToken, rp_coll);

    for (edm::DetSetVector<TotemRPLocalTrack>::const_iterator rp=rp_coll->begin(); rp!=rp_coll->end(); rp++) {
      const unsigned int det_id = rp->detId();
      const unsigned short arm  = (det_id%100==2), // 0->F, 1->N (003/103->F, 002/102->N)
                           side = (det_id/100);    // 0->L, 1->R (002/003->L, 102/103->R)
      aNumProtons = 0;
      for (edm::DetSet<TotemRPLocalTrack>::const_iterator proton=rp->begin(); proton!=rp->end(); proton++) {
        if (!proton->isValid()) continue;
        aProtonX[aNumProtons] = (proton->getX0())/1.e3;
        aProtonY[aNumProtons] = (proton->getY0())/1.e3;
        aProtonZ[aNumProtons] = (proton->getZ0())/1.e3;
        aProtonXsigma[aNumProtons] = (proton->getX0Sigma())/1.e3;
        aProtonYsigma[aNumProtons] = (proton->getY0Sigma())/1.e3;
        aProtonArm[aNumProtons] = arm;
        aProtonSide[aNumProtons] = side;
      }
    }
  }

  fTree->Fill();

}

bool
GammaGammaGammaGamma::passPhotonId(const edm::Ptr< pat::Photon >& photon_ref) const
{
  //const bool pass_medium = (*fPhotonMediumIdDecisions)[photon_ref->photonCore()];
  const bool pass_medium = true;

  return pass_medium;
}


// ------------ method called once each job just before starting event loop  ------------
void 
GammaGammaGammaGamma::beginJob()
{
  fTree->Branch("Run", &aRunId, "Run/I");
  fTree->Branch("LumiSection", &aLumiSection, "LumiSection/I");
  fTree->Branch("EventNum", &aEventNum, "EventNum/I");
  fTree->Branch("nPhotonCand", &aNumPhotons, "nPhotonCand/I");
  fTree->Branch("PhotonCand_pt", aPhotonPt, "PhotonCand_pt[nPhotonCand]/D");
  fTree->Branch("PhotonCand_eta", aPhotonEta, "PhotonCand_eta[nPhotonCand]/D");
  fTree->Branch("PhotonCand_phi", aPhotonPhi, "PhotonCand_phi[nPhotonCand]/D");
  fTree->Branch("PhotonCand_vtxx", aPhotonVtxX, "PhotonCand_vtxx[nPhotonCand]/D");
  fTree->Branch("PhotonCand_vtxy", aPhotonVtxY, "PhotonCand_vtxy[nPhotonCand]/D");
  fTree->Branch("PhotonCand_vtxz", aPhotonVtxZ, "PhotonCand_vtxz[nPhotonCand]/D");
  fTree->Branch("PhotonCand_MVAId_value", aPhotonIdMVAvalue, "PhotonCand_MVAId_value[nPhotonCand]/D");
  fTree->Branch("PhotonCand_MVAId_category", aPhotonIdMVAcategory, "PhotonCand_MVAId_category[nPhotonCand]/I");
  fTree->Branch("PhotonCand_trackIso", aPhotonIsoTrack, "PhotonCand_trackIso[nPhotonCand]/D");
  fTree->Branch("PhotonCand_ecalIso", aPhotonIsoECAL, "PhotonCand_ecalIso[nPhotonCand]/D");
  fTree->Branch("PhotonCand_hcalIso", aPhotonIsoHCAL, "PhotonCand_hcalIso[nPhotonCand]/D");
  fTree->Branch("PhotonCand_caloIso", aPhotonIsoCalo, "PhotonCand_caloIso[nPhotonCand]/D");
  fTree->Branch("PhotonCand_r9", aPhotonR9, "PhotonCand_r9[nPhotonCand]/D");
  fTree->Branch("PhotonCand_electronVeto", aPhotonElectronVeto, "PhotonCand_electronVeto[nPhotonCand]/D");
  fTree->Branch("PhotonCand_converted", aPhotonConverted, "PhotonCand_converted[nPhotonCand]/D");

  fTree->Branch("nPhotonPairCand", &aNumPhotonPairs, "nPhotonPairCand/I");
  fTree->Branch("PhotonPairCand_vertex_dist", aPhotonPairVertexDist, "PhotonPairCand_vertex_dist[nPhotonPairCand]/D");
  fTree->Branch("PhotonPairCand_mass", aPhotonPairMass, "PhotonPairCand_mass[nPhotonPairCand]/D");
  fTree->Branch("PhotonPairCand_pt", aPhotonPairPt, "PhotonPairCand_pt[nPhotonPairCand]/D");
  fTree->Branch("PhotonPairCand_dphi", aPhotonPairDphi, "PhotonPairCand_dphi[nPhotonPairCand]/D");
  fTree->Branch("PhotonPairCand_dpt", aPhotonPairDpt, "PhotonPairCand_dpt[nPhotonPairCand]/D");

  if (fFetchProtons) {
    fTree->Branch("nLocalProtCand", &aNumProtons, "nLocalProtCand/I");
    fTree->Branch("LocalProtCand_x", aProtonX, "LocalProtCand_x[nLocalProtCand]/D");
    fTree->Branch("LocalProtCand_y", aProtonY, "LocalProtCand_y[nLocalProtCand]/D");
    fTree->Branch("LocalProtCand_z", aProtonZ, "LocalProtCand_z[nLocalProtCand]/D");    
    fTree->Branch("LocalProtCand_xSigma", aProtonXsigma, "LocalProtCand_xSigma[nLocalProtCand]/D");
    fTree->Branch("LocalProtCand_ySigma", aProtonYsigma, "LocalProtCand_ySigma[nLocalProtCand]/D");
    fTree->Branch("LocalProtCand_arm", aProtonArm, "LocalProtCand_arm[nLocalProtCand]/I");
    fTree->Branch("LocalProtCand_side", aProtonSide, "LocalProtCand_side[nLocalProtCand]/I");    
  }
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
