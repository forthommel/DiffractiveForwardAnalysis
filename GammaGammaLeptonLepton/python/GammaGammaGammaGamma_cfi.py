import FWCore.ParameterSet.Config as cms

gggg_aod = cms.EDAnalyzer('GammaGammaGammaGamma',
    photonTag = cms.InputTag('selectedPatPhotons'),
    beamspotTag = cms.InputTag('offlineBeamSpot'),
    conversionTag = cms.InputTag('particleFlowEGamma', 'conversions'),
    pfCandidateTag = cms.InputTag('PFCandidates'),
    vertexTag = cms.InputTag('offlinePrimaryVertices'),
    protonTag = cms.InputTag('totemRPLocalTrackFitter'),
    photonMedIdBoolMapTag = cms.InputTag('egmPhotonIDs:mvaPhoID-PHYS14-PU20bx25-nonTrig-V1-wp90'),
    photonMedIdInfoMapTag = cms.InputTag('egmPhotonIDs:mvaPhoID-PHYS14-PU20bx25-nonTrig-V1-wp90'),
    photonMinPt = cms.double(75.),
)

gggg_miniaod = gggg_aod.clone(
    photonTag = cms.InputTag('slimmedPhotons'),
    conversionTag = cms.InputTag('reducedEgamma', 'reducedConversions'),
    pfCandidateTag = cms.InputTag('packedPFCandidates'),
    vertexTag = cms.InputTag('offlineSlimmedPrimaryVertices'),
)

#from RecoEgamma.PhotonIdentification.PhotonIDValueMapProducer_cfi import photonIDValueMapProducer
#from RecoEgamma.PhotonIdentification.PhotonMVAValueMapProducer_cfi import photonMVAValueMapProducer

#gggg = cms.Sequence(photonIDValueMapProducer * photonMVAValueMapProducer * gggg_aod)
