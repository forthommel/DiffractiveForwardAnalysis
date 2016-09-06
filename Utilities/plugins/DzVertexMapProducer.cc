#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h" // PAT

#include "DiffractiveForwardAnalysis/Utilities/interface/VertexCandidateMap.h"

class DzVertexMapProducer : public edm::EDProducer
{
 public:
  DzVertexMapProducer(const edm::ParameterSet&);

 private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  edm::EDGetTokenT< edm::View<reco::Vertex> > vertexToken_;
  edm::EDGetTokenT< edm::View<pat::PFParticle> > pfcandidateToken_;
  double maxAllowedDz_;
  bool useEachTrackOnce_;
};

DzVertexMapProducer::DzVertexMapProducer(const edm::ParameterSet& iConfig) :
  vertexToken_     (consumes< edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("VertexTag"))),
  pfcandidateToken_(consumes< edm::View<pat::PFParticle> >(iConfig.getParameter<edm::InputTag>("PFCandidatesTag"))),
  maxAllowedDz_    (iConfig.getParameter<double>("MaxAllowedDz")), // in cm
  useEachTrackOnce_(iConfig.getParameter<bool>("UseEachTrackOnce"))
{
  produces<VertexCandidateMap>();
}

void
DzVertexMapProducer::produce(edm::Event& evt, const edm::EventSetup&)
{
  edm::Handle< edm::View<reco::Vertex> > primaryVertices;
  evt.getByToken( vertexToken_, primaryVertices );

  edm::Handle< edm::View<pat::PFParticle> > pfCandidates;
  evt.getByToken( pfcandidateToken_, pfCandidates );

  std::auto_ptr<VertexCandidateMap> assoc(new VertexCandidateMap);

  if (useEachTrackOnce_) {
    // Associate a track to the closest vertex only, and only if dz < maxAllowedDz_
    for (unsigned int i=0; i<pfCandidates->size(); i++) {
      edm::Ptr<pat::PFParticle> cand = pfCandidates->ptrAt(i);
      if (cand->charge()==0) { continue; } // skip neutrals

      double closestDz = maxAllowedDz_;
      unsigned int closestDzIndex = -1;

      for (unsigned int j=0; j<primaryVertices->size(); j++) {
        edm::Ptr<reco::Vertex> vtx = primaryVertices->ptrAt(j);

        const double dz = fabs(cand->originalObject()->dz(vtx->position()));
        if (dz<closestDz) {
          closestDz = dz;
          closestDzIndex = j;
        }
      }
      if (closestDz<maxAllowedDz_) {
        edm::Ptr<reco::Vertex> vtx = primaryVertices->ptrAt(closestDzIndex);
        assoc->emplace_back(vtx, cand);
      }
    }
  }
  else { /* i.e. if !useEachTrackOnce_ */
    // Allow a track to be associated to multiple vertices if it's close to each of them
    for (unsigned int i=0; i<pfCandidates->size(); i++) {
      edm::Ptr<pat::PFParticle> cand = pfCandidates->ptrAt(i);

      if (cand->charge()==0) { continue; } // skip neutrals

      for (unsigned int j=0; j<primaryVertices->size(); j++) {
        edm::Ptr<reco::Vertex> vtx = primaryVertices->ptrAt(j);

        const double dz = fabs( cand->originalObject()->dz(vtx->position()));
        if (dz<maxAllowedDz_) {
          assoc->emplace_back(vtx, cand);
        }
      }
    }
  }

  std::stable_sort( assoc->begin(), assoc->end(), compare_by_vtx() );

  evt.put( assoc );
}

DEFINE_FWK_MODULE( DzVertexMapProducer );
