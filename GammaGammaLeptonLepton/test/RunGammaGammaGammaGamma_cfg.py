import FWCore.ParameterSet.Config as cms

process = cms.Process("gggg")

runOnMC = False
useAOD = True # AOD or MiniAOD?

#########################
#    General options    #
#########################

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options   = cms.untracked.PSet(
    #wantSummary = cms.untracked.bool(True),
    #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#########################
#      Input files      #
#########################

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#      '/store/relval/CMSSW_8_0_14/DoubleEG/RAW-RECO/ZElectron-80X_dataRun2_relval_v15_RelVal_doubEG2016B-v1/10000/1ADA6A0A-EB4A-E611-AF74-0025905A612C.root',
#        '/store/relval/CMSSW_8_1_0_pre9/RelValProdMinBias_13/AODSIM/81X_mcRun2_asymptotic_v2-v1/10000/0E584269-FE52-E611-993E-0CC47A4C8E0E.root',
      #'/store/mc/Run2015D/DoubleEG/AOD/04Dec2015-v1/10000/04D11E1B-BB9E-E511-AC1A-047D7B881D62.root',
      #'/store/data/Run2016B/DoubleEG/AOD/PromptReco-v2/000/273/158/00000/006772B7-E019-E611-AEBE-02163E014583.root',
      #'/store/data/Run2016C/DoubleEG/AOD/PromptReco-v2/000/275/603/00000/0C964C97-893A-E611-B3D2-02163E0128F1.root',
      #'/store/data/Run2016C/DoubleEG/AOD/PromptReco-v2/000/275/601/00000/D211B105-753A-E611-A7B0-02163E013668.root',
      #'/store/data/Run2016C/DoubleEG/AOD/PromptReco-v2/000/275/657/00000/0AE1FDE0-673B-E611-9ECF-02163E011A27.root',
#'/store/data/Run2016C/DoubleEG/AOD/PromptReco-v2/000/275/657/00000/0AE1FDE0-673B-E611-9ECF-02163E011A27.root',
'/store/data/Run2016B/DoubleEG/AOD/01Jul2016-v1/00000/04D8B592-B745-E611-80F9-02163E0148EE.root'
    ),
    #firstEvent = cms.untracked.uint32(0)
)

#########################
#        Triggers       #
#########################

process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.HLTPaths = cms.vstring(
    'HLT_DoublePhoton85_v*', 'HLT_DoublePhoton60_v*',
)

#########################
#      Preskimming      #
#########################
process.load("Configuration.StandardSequences.GeometryDB_cff") ## FIXME need to ensure that this is the good one
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4),
    maxAbsZ = cms.double(15),
    maxd0 = cms.double(2)
)

#########################
#     PAT-ification     #
#########################
## Look at https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#Core_Tools for more information

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from Configuration.EventContent.EventContent_cff import *

#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning, patExtraAodEventContent
#from PhysicsTools.PatAlgos.tools.coreTools import *

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_offline*PrimaryVertices*_*_*',
        'keep *_gsfElectrons*_*_*',
        'keep *_selectedPatPhotons*_*_*',
        'keep *_*hoton*_*_*',
        'keep *_selectedPatJets*_*_*',
        'keep *_*MET*_*_*',
        'keep *_*particleFlow*_*_*',
        #*patEventContentNoCleaning
    ),
)

from DiffractiveForwardAnalysis.GammaGammaLeptonLepton.RemovePATMCMatching_cfi import removePATMCMatching

if not runOnMC:
    #names = ['Photons', 'Electrons', 'Muons', 'Jets', 'METs']
    removePATMCMatching(process)

#########################
#       Photon ID       #
#########################

from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDPhotonIdProducer, setupVIDPhotonSelection, setupAllVIDIdsInModule, DataFormat

switchOnVIDPhotonIdProducer(process, DataFormat.AOD)
#setupAllVIDIdsInModule(process, 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_PHYS14_PU20bx25_nonTrig_V1_cff', setupVIDPhotonSelection)
#setupAllVIDIdsInModule(process, 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2_cff', setupVIDPhotonSelection)
setupAllVIDIdsInModule(process, 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2p1_cff', setupVIDPhotonSelection)

#########################
#       Analysis        #
#########################

process.load("DiffractiveForwardAnalysis.GammaGammaLeptonLepton.GammaGammaGammaGamma_cfi")

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("output.root"),
    closeFileFast = cms.untracked.bool(True)
)

process.gggg_aod.TriggersList = process.hltFilter.HLTPaths
process.gggg_aod.RunOnMC = cms.untracked.bool(runOnMC)
process.gggg_aod.fetchProtons = cms.bool(True)

process.p = cms.Path(
    process.hltFilter*
    process.egmPhotonIDSequence*
    process.patDefaultSequence*
    process.gggg_aod
)
