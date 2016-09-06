#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h" // RECO
//#include "DataFormats/PatCandidates/interface/PFParticle.h" // PAT
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // miniAOD
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DiffractiveForwardAnalysis/Utilities/interface/VertexCandidateMap.h"

// Shamefully stolen from flashgg's VertexMapFromCandidateProducer
// (https://github.com/cms-analysis/flashgg/blob/master/MicroAOD/plugins/VertexMapFromCandidateProducer.cc)


// use PackedCandidate::fromPV to produce a vertex map
// Require fromPV(int) > FromPVCut
// From https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#PV_Assignment
//   "The definition normally used for isolation calculations is fromPV() > 1;
//    the definition used for CHS subtraction in jets is fromPV() > 0."

class VertexMapFromCandidateProducer : public edm::EDProducer
{
 public:
  VertexMapFromCandidateProducer(const edm::ParameterSet&);

 private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  edm::EDGetTokenT< edm::View<reco::Vertex> > vertexToken_;
  edm::EDGetTokenT< edm::View<reco::PFCandidate> > pfcandidateToken_;
  unsigned int fromPVgt_, fromPVgtIfDz_;
  double dzCut_;
};

VertexMapFromCandidateProducer::VertexMapFromCandidateProducer(const edm::ParameterSet& iConfig) :
  vertexToken_     (consumes< edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("VertexTag"))),
  pfcandidateToken_(consumes< edm::View<reco::PFCandidate> >(iConfig.getParameter<edm::InputTag>("PFCandidatesTag"))), // RECO
  //pfcandidateToken_(consumes< edm::View<pat::PFParticle> >(iConfig.getParameter<edm::InputTag>("PFCandidatesTag"))), // PAT
  //pfcandidateToken_(consumes< edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("PFCandidatesTag"))), // miniAOD
  fromPVgt_        (iConfig.getParameter<unsigned int>("FromPVCut")), // 0 for standard CHS, 2 for PUPPI
  fromPVgtIfDz_    (iConfig.getParameter<unsigned int>("FromPVCutIfPassDz")), // 0 for CHS or PUPPI
  dzCut_           (iConfig.getParameter<double>("DzCut"))
{
  produces<VertexCandidateMap>();
}

void
VertexMapFromCandidateProducer::produce(edm::Event& evt, const edm::EventSetup&)
{
  edm::Handle< edm::View<reco::Vertex> > primaryVertices;
  evt.getByToken( vertexToken_, primaryVertices );

  edm::Handle< edm::View<reco::PFCandidate> > pfCandidates; // RECO
  //edm::Handle< edm::View<pat::PFParticle> > pfCandidates; // PAT
  //edm::Handle< edm::View<pat::PackedCandidate> > pfCandidates; // miniAOD
  evt.getByToken( pfcandidateToken_, pfCandidates );

  std::auto_ptr<VertexCandidateMap> assoc( new VertexCandidateMap );

  for (unsigned int i=0; i<pfCandidates->size(); i++) {
    //edm::Ptr<pat::PackedCandidate> cand = pfCandidates->ptrAt(i);
    edm::Ptr<reco::PFCandidate> cand = pfCandidates->ptrAt(i);
    if (cand->charge()==0) { continue; } // skip neutrals

    for (unsigned int j=0; j<primaryVertices->size(); j++) {
      if (cand->fromPV(j)>fromPVgt_) {
        assoc->emplace_back(primaryVertices->ptrAt(j), cand);
      }
      else if (cand->fromPV(j)>fromPVgtIfDz_) {
        // This section to support extra cut on dZ in puppi code
        const double absdz = fabs(cand->dz(primaryVertices->ptrAt(j)->position()));
        if (absdz<dzCut_) {
          assoc->emplace_back(primaryVertices->ptrAt(j), cand);
        }
      }
    }
  }

  std::stable_sort(assoc->begin(), assoc->end(), compare_by_vtx());

  evt.put(assoc);
}

DEFINE_FWK_MODULE(VertexMapFromCandidateProducer);
