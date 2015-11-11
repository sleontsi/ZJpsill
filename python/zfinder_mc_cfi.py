import FWCore.ParameterSet.Config as cms

process = cms.Process("ZFinder")

# Set up message output and logging
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'START53_V29B::All' ##mc

#process.MessageLogger.cerr.FwkReport.reportEvery = 100  # Report status ever 100 events
process.MessageLogger.cerr.FwkReport.reportEvery = 1000  # Report status ever 100 events

# Number of events from each file to process. It should be -1 (all) when
# running for an analysis
N_EVENTS_TO_PROCESS = -1
if N_EVENTS_TO_PROCESS != -1:
    print "NOT RUNNING ON ALL EVENTS IN THE FILE!"
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(N_EVENTS_TO_PROCESS)
        )

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 'file:/hdfs/cms/user/turkewitz/tmp/jpsi_inclusive_mc_test_file.root')
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string("test9mc.root")
        )

#
# rho value for isolation
#

from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets  # the 4 references the rParam = 0.4
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

#
# particle flow isolation
#

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.pfiso = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence)

#
# ZFinder
#

# Import ZDefinitions
from ZFinder.Event.ZDefinitions_cfi import zdefs

process.ZFinder = cms.EDAnalyzer('ZFinder',
        # General tags
        ecalElectronsInputTag  = cms.InputTag("gsfElectrons"),
        muonsInputTag          = cms.InputTag("muons"),
        hfElectronsInputTag    = cms.InputTag("hfRecoEcalCandidate"),
        hfClustersInputTag     = cms.InputTag("hfEMClusters"),
        conversionsInputTag    = cms.InputTag("allConversions"),
        beamSpotInputTag       = cms.InputTag("offlineBeamSpot"),
        rhoIsoInputTag         = cms.InputTag("kt6PFJetsForIsolation", "rho"),
        primaryVertexInputTag  = cms.InputTag("offlinePrimaryVertices"),
        ntElectronsInputTag    = cms.InputTag("photons"),
        ak5PFJetsInputTag      = cms.InputTag("ak5PFJets"),
        isoValInputTags        = cms.VInputTag(
            cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
            cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
            cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')
            ),
        # MC, but still required to be something for data
        pileupInputTag = cms.InputTag("addPileupInfo"),
        generatorInputTag = cms.InputTag("genParticles"),
        # ZDefinitions from ZFinder.ZFinder.ZDefinitions_cfi
        ZDefinitions = zdefs,
        ##pileup_era = cms.string("ABCD") # defaults to ABCD
        pileup_era = cms.string("B") # defaults to ABCD
        )

# RUN
process.p = cms.Path(process.kt6PFJetsForIsolation * process.pfiso * process.ZFinder)
process.schedule = cms.Schedule(process.p)
