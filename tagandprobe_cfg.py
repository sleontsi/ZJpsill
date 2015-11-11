import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

# Set up message output and logging
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

##TODO make this an option, maybe make a different file for data/mc to avoid forgetting
##about to switch out the glbal tag
process.GlobalTag.globaltag = 'FT53_V21A_AN6::All' ##data
#process.GlobalTag.globaltag = 'START53_V29B::All' ##mc

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
    ##fileNames = cms.untracked.vstring( 'file:/local/cms/user/turkewitz/ZPhysics/JpsiSkim/jpsiSkimUpdated.root')
    ##fileNames = cms.untracked.vstring( 'file:/hdfs/cms/user/turkewitz/ZPhysics/JPsiSkim/DoubleElectron/2012B/jpsiSkimUpdated/jpsiSkimUpdated_999-pool.root')
    ##fileNames = cms.untracked.vstring( 'file:/hdfs/cms/user/turkewitz/ZPhysics/JPsiSkim/MuOnia/Run2012B/jpsiTriggerSkim_198_1_U0I.root')
    ##fileNames = cms.untracked.vstring( 'file:/local/cms/user/turkewitz/ZPhysics/PromptJpsi/JpsiMM_8TeV_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_RECO.root')
    ##fileNames = cms.untracked.vstring( 'file:/local/cms/user/turkewitz/ZPhysics/tmp/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_00037C53-AAD1-E111-B1BE-003048D45F38.root')
    ##fileNames = cms.untracked.vstring( 'file:/hdfs/cms/user/turkewitz/tmp/jpsi_inclusive_mc_test_file.root')
    ###fileNames = cms.untracked.vstring( 'file:/local/cms/user/turkewitz/ZPhysics/tmp/0012F37A-CE09-E211-ABDA-00261894396F.root')
    ##fileNames = cms.untracked.vstring( 'file:/local/cms/user/turkewitz/MuonTriggerSkimTest2012B/MuonTriggerSkimTest2012B_200-pool.root')
    fileNames = cms.untracked.vstring( 'file:/local/cms/user/turkewitz/ZPhysics/tmp/6642D763-9872-E211-BF3C-00259074AE5C.root')
    ##fileNames = cms.untracked.vstring( 'file:/home/user1/turkewitz/Work/CMSSW_5_3_13_ZJPsi/src/jpsiMuMu_Zee_Skim.root')
    ##fileNames = cms.untracked.vstring( 'file:/home/user1/turkewitz/Work/CMSSW_5_3_13_ZJPsi/src/jpsiMuMu_Zmumu_Skim.root')
    #fileNames = cms.untracked.vstring( 'file:/home/user1/turkewitz/Work/CMSSW_5_3_13_ZJPsi/src/jpsiSkimMuonsUpdated.root')
    ##fileNames = cms.untracked.vstring( 'root://xrootd.unl.edu//store/mc/Summer12_DR53X/JPsiToMuMu_2MuPtEtaFilter_tuneD6T_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7A-v2/00000/0012F37A-CE09-E211-ABDA-00261894396F.root')
    ## fileNames = cms.untracked.vstring( 'file:/hdfs/cms/user/turkewitz/ZPhysics/JPsiSkim/DoubleElectron/jpsiSkim/jpsiSkim_000-pool.root')
    ##fileNames = cms.untracked.vstring( 'file:/home/user1/turkewitz/Work/CMSSW_5_3_13_ZJPsi/src/JPsiFilter/jpsiSkim.root')
    ##fileNames = cms.untracked.vstring( 'file:/local/cms/phedex/store/data/Run2012A/DoubleElectron/AOD/22Jan2013-v1/20000/003EC246-5E67-E211-B103-00259059642E.root')
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string("test9e2c.root")
        )

# Run only on lumis specified in the lumi file
# Recipe from:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePythonTips#Use_a_JSON_file_of_good_lumi_sec
from FWCore.ParameterSet.Types import untracked, VLuminosityBlockRange
from FWCore.PythonUtilities.LumiList import LumiList
##json_file for electrons
json_file = "/home/user1/turkewitz/Work/CMSSW_5_3_13_ZJPsi/src/Metadata/lumi_json/Run2012ABCD.json" # File location
run_2012abcd_lumis = LumiList(filename = json_file).getCMSSWString().split(',')
process.source.lumisToProcess = untracked(VLuminosityBlockRange(run_2012abcd_lumis))

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
#from ZFinder.Event.ZDefinitions_cfi import zdefs

process.TagAndProbe = cms.EDAnalyzer('TagAndProbe',
        # General tags
        ecalElectronsInputTag  = cms.InputTag("gsfElectrons"),
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
        ##pileup_era = cms.string("ABCD") # defaults to ABCD
        pileup_era = cms.string("B") # defaults to ABCD
        )

# RUN
process.p = cms.Path(process.kt6PFJetsForIsolation * process.pfiso * process.TagAndProbe)
process.schedule = cms.Schedule(process.p)
