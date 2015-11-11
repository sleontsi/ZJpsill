// -*- C++ -*-
// TODO update this description
//
// Package:    TagAndProbe
// Class:      TagAndProbe
//
/**\class TagAndProbe TagAndProbe.cc ZShape/TagAndProbe/src/TagAndProbe.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Alexander Gude
//         Created:  Thu Aug  8 15:19:00 CDT 2013
// $Id$
//
//


// system include files
#include <memory>

// standard library files
#include <map>  // std::map
#include <string>  // std::string
#include <utility>  // std::pair
#include <algorithm>  // std::sort, std::swap
#include <iostream>  // std::cout, std::endl

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// CMSSW
#include "FWCore/ServiceRegistry/interface/Service.h" // edm::Service
#include "CommonTools/UtilAlgos/interface/TFileService.h" // TFileService
#include "DataFormats/Common/interface/Handle.h"  // edm::Handle
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"  // reco::PhotonCollection
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"  // reco::SuperClusterCollection, reco::SuperClusterRef
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"  // reco::RecoEcalCandidateCollection
#include "DataFormats/MuonReco/interface/MuonFwd.h" // reco::MuonCollection
#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"  // EgammaCutBasedEleId::PassWP, EgammaCutBasedEleId::*
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h" // GenEventInfoProduct
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"  // PileupSummaryInfo
#include "DataFormats/HLTReco/interface/TriggerEvent.h" // trigger::TriggerEvent
#include "DataFormats/TrackReco/interface/Track.h" //reco::Track
#include "DataFormats/TrackReco/interface/TrackFwd.h" // reco::TrackCollection
#include "DataFormats/JetReco/interface/PFJetCollection.h" //
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

// ZFinder 
#include "ZFinder/Event/interface/ZDefinition.h"  // ZDefinition
#include "ZFinder/Event/interface/ZFinderEvent.h"  // ZFinderEvent
#include "ZFinder/Event/interface/ZFinderPlotter.h"  // ZFinderPlotter
#include "ZFinder/Event/interface/ZFinderCuts.h"  // ZFinderCuts

// CMSSW

//
// class declaration
//

class TagAndProbe : public edm::EDAnalyzer {
  public:
    explicit TagAndProbe(const edm::ParameterSet&);
    ~TagAndProbe();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

    // ----------member data ---------------------------
    const edm::ParameterSet& iConfig_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TagAndProbe::TagAndProbe(const edm::ParameterSet& iConfig) : iConfig_(iConfig) {
  //now do what ever initialization is needed

  edm::Service<TFileService> fs;

}

TagAndProbe::~TagAndProbe() {
}


//
// member functions
//

// ------------ method called for each event  ------------
void TagAndProbe::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace zf;

  zf::ZFinderEvent zfe(iEvent, iSetup, iConfig_);

  //edm::Handle<reco::Track> generalTracks;
  edm::Handle<reco::TrackCollection> generalTracks;
  if ( !iEvent.getByLabel("generalTracks", generalTracks) ) {
    std::cout << "ERROR: Could not find generalTracks" << std::endl ;
    return ;
  }

  for (unsigned int i = 0; i < zfe.reco_jpsi.m.size() ; ++i ) {
    for (size_t j = 0; j < generalTracks->size(); ++j) {
      const reco::Track & track = generalTracks->at(j);
      if (track.pt() > 2) {
        std::cout << track.pt() << std::endl;
      }
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void TagAndProbe::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void TagAndProbe::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void TagAndProbe::beginRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a run  ------------
void TagAndProbe::endRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block  ------------
void TagAndProbe::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void TagAndProbe::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TagAndProbe::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TagAndProbe);
