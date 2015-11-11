import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

# Set up message output and logging
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
##process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'FT53_V21A_AN6::All' ##data
#process.GlobalTag.globaltag = 'START53_V29B::All' ##mc

process.MessageLogger.cerr.FwkReport.reportEvery = 1000  # Report status ever 1000 events

# Number of events from each file to process. It should be -1 (all) when
# running for an analysis
N_EVENTS_TO_PROCESS = 1000
if N_EVENTS_TO_PROCESS != -1:
    print "NOT RUNNING ON ALL EVENTS IN THE FILE!"
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(N_EVENTS_TO_PROCESS)
        )

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
# di-muon AOD
#    fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/data/Run2012C/DoubleMuParked/AOD/22Jan2013-v1/10000/0002ACB4-C96C-E211-A96F-20CF3027A628.root')
# di-electron AOD
  fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/data/Run2012C/DoubleElectron/AOD/22Jan2013-v1/20000/00B50FA7-E867-E211-9F1B-003048FFD754.root')
#  fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/data/Run2012C/DoubleElectron/AOD/22Jan2013-v1/20000/1E64EF26-0268-E211-9B23-002618943882.root')
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string("zfinder.root")
        )

# Run only on lumis specified in the lumi file
# Recipe from:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePythonTips#Use_a_JSON_file_of_good_lumi_sec
from FWCore.ParameterSet.Types import untracked, VLuminosityBlockRange
from FWCore.PythonUtilities.LumiList import LumiList
##json_file for electrons
json_file = "/afs/cern.ch/user/s/sleontsi/Files/Cert_Run2012ABCD_22Jan2013ReReco.json" # File location
run_2012abcd_lumis = LumiList(filename = json_file).getCMSSWString().split(',')
process.source.lumisToProcess = untracked(VLuminosityBlockRange(run_2012abcd_lumis))

#
# electron regression
#

from ZFinder.Event.electron_regression_cfi import CalibratedElectrons, RandomNumberGeneratorService, ElectronEnergyRegressions
process.RandomNumberGeneratorService = RandomNumberGeneratorService
process.CalibratedElectrons = CalibratedElectrons
process.eleRegressionEnergy = ElectronEnergyRegressions

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
##process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.eleIsoSequence = setupPFElectronIso(process, 'CalibratedElectrons:calibratedGsfElectrons')
process.pfiso = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence)

##TESTING muon pfIso
##process.muonIsoSequence = setupPFMuonIso(process, 'muons')
##process.pfiso = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence + process.muonIsoSequence)

#
# ZFinder
#

# Import ZDefinitions
#from ZFinder.Event.ZDefinitions_cfi import zdefs

process.ZFinder = cms.EDAnalyzer('ZFinder',
        # General tags
        # Use the calibrated electrons we make with process.CalibratedElectrons
        ecalElectronsInputTag = cms.InputTag("CalibratedElectrons", "calibratedGsfElectrons"),
        ##ecalElectronsInputTag  = cms.InputTag("gsfElectrons"),
        muonsInputTag          = cms.InputTag("muons"),
        conversionsInputTag    = cms.InputTag("allConversions"),
        beamSpotInputTag       = cms.InputTag("offlineBeamSpot"),
        rhoIsoInputTag         = cms.InputTag("kt6PFJetsForIsolation", "rho"),
        primaryVertexInputTag  = cms.InputTag("offlinePrimaryVertices"),
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
        #ZDefinitions = zdefs,
        pileup_era = cms.string("ABCD") # defaults to ABCD
        ##pileup_era = cms.string("B") # defaults to ABCD
        )

# RUN
process.p = cms.Path(process.kt6PFJetsForIsolation * process.eleRegressionEnergy * process.CalibratedElectrons * process.pfiso * process.ZFinder)
process.schedule = cms.Schedule(process.p)
