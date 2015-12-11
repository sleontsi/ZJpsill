// -*- C++ -*-
//
// Package:    ZFinder
// Class:      ZFinder
//
/**\class ZFinder ZFinder.cc ZShape/ZFinder/src/ZFinder.cc

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

// ZFinder
#include "ZFinder/Event/interface/ZDefinition.h"  // ZDefinition
#include "ZFinder/Event/interface/ZDefinitionPlotter.h"  // ZDefinitionPlotter
#include "ZFinder/Event/interface/ZDefinitionWorkspace.h"  // ZDefinitionWorkspace
#include "ZFinder/Event/interface/ZFinderEvent.h"  // ZFinderEvent
#include "ZFinder/Event/interface/ZFinderPlotter.h"  // ZFinderPlotter
#include "ZFinder/Event/interface/ZFinderTree.h"  // ZFinderTree
#include "ZFinder/Event/interface/ZFinderCuts.h"  // ZFinderCuts

// HLT information
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

//
// class declaration
//

class ZFinder : public edm::EDAnalyzer {
  public:
    explicit ZFinder(const edm::ParameterSet&);
    ~ZFinder();

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
    zf::ZFinderTree *z_tree;

    zf::ZFinderPlotter *zfp_all;
    zf::ZFinderPlotter *zfp_dimuon_jpsi, *zfp_dimuon_jpsi_soft, *zfp_dimuon_jpsi_vtx_compatible, *zfp_dimuon_jpsi_primary_vertex, *zfp_jpsi, *zfp_prompt_jpsi;
    zf::ZFinderPlotter *zfp_dielectron_z, *zfp_dielectron_z_good, *zfp_dielectron_z_good_compatible_vertex, *zfp_z_to_electrons, 
      *zfp_z_to_electrons_and_good_dimuon_jpsi, *zfp_z_to_electrons_and_jpsi, *zfp_z_to_electrons_and_prompt_jpsi;
    zf::ZFinderPlotter *zfp_dimuon_z, *zfp_dimuon_z_good, *zfp_dimuon_z_good_compatible_vertex, *zfp_z_to_muons, 
      *zfp_z_to_muons_and_good_dimuon_jpsi, *zfp_z_to_muons_and_jpsi, *zfp_z_to_muons_and_prompt_jpsi;
    zf::ZFinderPlotter *zfp_all_mc, *zfp_jpsi_mc;
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
ZFinder::ZFinder(const edm::ParameterSet& iConfig) : iConfig_(iConfig) {

  // Setup plotters
  edm::Service<TFileService> fs;

  TFileDirectory tdir_tree(fs->mkdir("Tree"));

  z_tree = new zf::ZFinderTree(tdir_tree, true);

}

ZFinder::~ZFinder() {
  //for (auto& i_zdeft : z_tuples_) {
  //  delete i_zdeft;
  //}
}


//
// member functions
//

// ------------ method called for each event  ------------
void ZFinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace zf;

//  if (iEvent.id().event() < 1081870)
//    return;
//
//  if (iEvent.id().event() == 1081878)
//   std::printf("------------------------> entered the bad event\n"); 
//
//  std::printf("%d\n", iEvent.id().event());

  zf::ZFinderEvent zfe(iEvent, iSetup, iConfig_);

  if (is_Jpsimumu) {
    if (!zfe.found_jpsi)                                                                        
      return;
  }   
  if (is_Zmumu) {                                                                                                                        
    if (!zfe.found_z_to_muons_mass)
      return;
  }
  if (is_Jpsiee) {
    if (!zfe.found_jpsi_from_electrons)
      return;
  }
  if (is_Zee) {
    if (!zfe.found_z_to_electrons_mass)
     return;
  }

  //Fill Tree
  z_tree->Fill(zfe);
}

// ------------ method called once each job just before starting event loop  ------------
void ZFinder::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void ZFinder::endJob() {
}

// ------------ method called when starting to processes a run  ------------
//void ZFinder::beginRun(edm::Run const&, edm::EventSetup const&) {
void ZFinder::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
}

// ------------ method called when ending the processing of a run  ------------
void ZFinder::endRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block  ------------
void ZFinder::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void ZFinder::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ZFinder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZFinder);
