import FWCore.ParameterSet.Config as cms

gggg_aod = cms.EDAnalyzer('GammaGammaGammaGamma',
    photonTag = cms.InputTag('patPhotons'),
    beamspotTag = cms.InputTag('offlineBeamSpot'),
    #conversionTag = cms.InputTag('particleFlowEGamma', 'conversions'),
    conversionTag = cms.InputTag('conversions'),
    electronTag = cms.InputTag('patElectrons'),
    recoElectronTag = cms.InputTag('gedGsfElectrons'),
    pfCandidateTag = cms.InputTag('PFCandidates'),
    vertexTag = cms.InputTag('offlinePrimaryVertices'),
    protonTag = cms.InputTag('totemRPLocalTrackFitter'),
    #photonMedIdBoolMapTag = cms.InputTag('egmPhotonIDs:mvaPhoID-PHYS14-PU20bx25-nonTrig-V1-wp90'),
    #photonMedIdInfoMapTag = cms.InputTag('egmPhotonIDs:mvaPhoID-PHYS14-PU20bx25-nonTrig-V1-wp90'),
    photonMedIdBoolMapTag = cms.InputTag('egmPhotonIDs:mvaPhoID-Spring15-25ns-nonTrig-V2p1-wp90'),
    photonMedIdInfoMapTag = cms.InputTag('egmPhotonIDs:mvaPhoID-Spring15-25ns-nonTrig-V2p1-wp90'),
    photonMVAIdValuesMapTag = cms.InputTag('photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values'),
    photonMVAIdCategoriesMapTag = cms.InputTag('photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Categories'),
    photonMinPt = cms.double(75.),
)

gggg_miniaod = gggg_aod.clone(
    photonTag = cms.InputTag('slimmedPhotons'),
    conversionTag = cms.InputTag('reducedEgamma', 'reducedConversions'),
    electronTag = cms.InputTag('slimmedElectrons'),
    pfCandidateTag = cms.InputTag('packedPFCandidates'),
    vertexTag = cms.InputTag('offlineSlimmedPrimaryVertices'),
)

#from RecoEgamma.PhotonIdentification.PhotonIDValueMapProducer_cfi import photonIDValueMapProducer
#from RecoEgamma.PhotonIdentification.PhotonMVAValueMapProducer_cfi import photonMVAValueMapProducer

#gggg = cms.Sequence(photonIDValueMapProducer * photonMVAValueMapProducer * gggg_aod)
