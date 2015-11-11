#ifndef ZFINDER_ZFINDEREVENT_H_
#define ZFINDER_ZFINDEREVENT_H_

// Standard Library
#include <map>  // std::map
#include <string>  // std::string
#include <utility>  // std::pair
#include <vector>  // std::vector

// CMSSW
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"  // reco::GsfElectron
#include "DataFormats/EgammaCandidates/interface/Photon.h"  // reco::Photon
#include "DataFormats/MuonReco/interface/Muon.h" // reco::Muon
#include "DataFormats/JetReco/interface/PFJet.h" //
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"  // reco::GenParticle
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"  // reco::RecoEcalCandidate
#include "FWCore/Framework/interface/Event.h"  // edm::Event, edm::EventSetup
#include "FWCore/ParameterSet/interface/ParameterSet.h"  // edm::ParameterSet
#include "FWCore/Utilities/interface/InputTag.h"  // edm::InputTag
#include "DataFormats/HLTReco/interface/TriggerObject.h"  // trigger::TriggerObject
#include "DataFormats/Common/interface/TriggerResults.h"  // trigger::TriggerResults
#include "FWCore/Common/interface/TriggerNames.h"  // trigger::TriggerNames
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h" // edm::LumiReWeighting

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

// ZFinder
#include "ZFinder/Event/interface/ZFinderElectron.h"  // ZFinderElectron, ZFinderElectron
#include "ZFinder/Event/interface/ZFinderCuts.h"  // ZFinderMuon, ZFinderMuon

//Math
#include <math.h>
#include <TLorentzVector.h>

namespace zf {

  // Cut level struct
  struct CutLevel{
    // Constructor sets all values to false
    CutLevel() {
      pass = false;
      t0p1_pass = false;
      t0p1_eff = 1.;
      t1p0_pass = false;
      t1p0_eff = 1.;
      event_weight = 1.;
    }
    bool pass;
    bool t0p1_pass;
    bool t1p0_pass;
    double t0p1_eff;
    double t1p0_eff;
    double event_weight;
  };


  // Used to match cut levels to names
  typedef std::pair<std::string, CutLevel> cutlevel_pair;
  // Used to pass around cut levels
  typedef std::vector<cutlevel_pair> cutlevel_vector;
  // Used to pass around trigger objects for matching
  typedef std::pair<const trigger::TriggerObject*, double> trig_dr_pair;
  typedef std::vector<trig_dr_pair> trig_dr_vec;


  class ZFinderEvent{
    public:
      // Constructor. Although iEvent, iSetup, and iConfig violate our naming
      // convention, they are almost ubiquitous in CMSSW code
      ZFinderEvent() {}
      ZFinderEvent(
          const edm::Event& iEvent,
          const edm::EventSetup& iSetup,
          const edm::ParameterSet& iConfig
          );
      // Destructor
      ~ZFinderEvent();

      // Data or MC
      bool is_real_data;

      // Beam Spot
      struct Beamspot{
        double x;
        double y;
        double z;
      } reco_bs;

      // Primary vertexes
      struct Vertexes{
        unsigned int num;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
        double primary_x;
        double primary_y;
        double primary_z;
        reco::Vertex primary_vert;
      } truth_vert, reco_vert;

      // Event ID
      struct EventID{
        unsigned int run_num;
        unsigned int lumi_num;
        unsigned int event_num;
      } id;

      // Z Data
      struct ZData{
        double m;
        double pt;
        double y;
        double phi;
        double phistar;
        double eta;
        std::vector<double> muon0_pT;
        std::vector<double> muon1_pT;
        std::vector<double> muon0_eta;
        std::vector<double> muon1_eta;
        std::vector<double> muon0_phi;
        std::vector<double> muon1_phi;
        std::vector<double> vtx_prob;
        double vtx_x;
        double vtx_y;
        double vtx_z;
        math::PtEtaPhiMLorentzVector zlv;
        TransientVertex vtx;
      } reco_z, reco_z_from_muons, truth_z_muons, truth_z_electrons;

      // JPsi Data
      struct JPsiData{
        std::vector<double> m;
        std::vector<double> pt;
        std::vector<double> y;
        std::vector<double> phistar;
        std::vector<double> eta;
        std::vector<double> phi;
        std::vector<double> tau_xy;
        std::vector<double> tau_z;
        std::vector<double> distance_x;
        std::vector<double> distance_y;
        std::vector<double> distance_z;
        std::vector<double> distance;
        std::vector<double> dist_err;
        std::vector<double> chi2;
        std::vector<double> distance_xy;
        std::vector<double> dist_err_xy;
        std::vector<double> chi2_xy;
        std::vector<double> vtx_x;
        std::vector<double> vtx_y;
        std::vector<double> vtx_z;
        std::vector<double> vtx_prob;
        std::vector<double> jpsi_efficiency;
        std::vector<double> jpsi_scale_factor;
        std::vector<double> cos_jpsi_mu_plus;
        std::vector<double> cos_jpsi_mu_minus;
        std::vector<double> muons_delta_phi;
        std::vector<double> muons_delta_eta;
        std::vector<double> muons_deltaR;
        std::vector<double> z_delta_phi;
        std::vector<double> four_lepton_mass;
        std::vector<reco::Muon> muon0;
        std::vector<reco::Muon> muon1;
        std::vector<double> muon0_pT;
        std::vector<double> muon1_pT;
        std::vector<double> muon0_eta;
        std::vector<double> muon1_eta;
        std::vector<double> muon0_phi;
        std::vector<double> muon1_phi;
        std::vector<double> jpsi_acc_eff;
        std::vector<double> muon0_efficiency;
        std::vector<double> muon1_efficiency;
        std::vector<double> muon0_scale_factor;
        std::vector<double> muon1_scale_factor;
        std::vector<double> muon0_deltaR_to_z_muons;
        std::vector<double> muon1_deltaR_to_z_muons;
        std::vector<double> muon0_deltaR_to_truth_muons;
        std::vector<double> muon1_deltaR_to_truth_muons;
        std::vector<double> iso_mu0;
        std::vector<double> iso_sum_charged_hadron_pt_mu0;
        std::vector<double> iso_sum_charged_particle_pt_mu0;
        std::vector<double> iso_sum_neutral_hadron_et_mu0;
        std::vector<double> iso_sum_photon_et_mu0;
        std::vector<double> iso_sum_pileup_pt_mu0;
        std::vector<double> iso_mu1;
        std::vector<double> iso_sum_charged_hadron_pt_mu1;
        std::vector<double> iso_sum_charged_particle_pt_mu1;
        std::vector<double> iso_sum_neutral_hadron_et_mu1;
        std::vector<double> iso_sum_photon_et_mu1;
        std::vector<double> iso_sum_pileup_pt_mu1;
        std::vector<double> trigger_object_mu0_pt;
        std::vector<double> trigger_object_mu1_pt;
        std::vector<bool> has_muons_in_eta_window;
        std::vector<bool> has_high_pt_muons;
        std::vector<bool> has_soft_id_muons;
        std::vector<bool> has_muons_with_compatible_vertex;
        std::vector<bool> has_dimuon_vertex_compatible_with_primary_vertex;
        std::vector<bool> is_high_pt;
        std::vector<bool> is_in_rap_window;
        std::vector<bool> is_within_jpsi_mass_window;
        std::vector<bool> is_prompt;
      } reco_jpsi, reco_jpsi_from_electrons, truth_jpsi;

      struct Jets{
        std::vector<double> pt;
        std::vector<double> phi;
        std::vector<double> eta;
        std::vector<double> btag_discriminator;
      } reco_jets, reco_muon_jets;
      // Event weight, used for things like pileup reweighting. Most
      // other weights are cut dependent (efficiencies for example) and
      // so are store with the cuts in each electron, with a combined
      // efficiency calculated in ZDefinition.
      double event_weight;


      // These are the special, selected electrons used to make the Z
      ZFinderElectron* e0;
      ZFinderElectron* e1;

      reco::Muon z_muon0;
      reco::Muon z_muon1;

      void set_e0(ZFinderElectron* electron) { e0 = electron; }
      void set_e1(ZFinderElectron* electron) { e1 = electron; }
      void set_both_e(ZFinderElectron* electron0, ZFinderElectron* electron1) { e0 = electron0; e1 = electron1; }
      ZFinderElectron* e0_truth;
      ZFinderElectron* e1_truth;
      void set_e0_truth(ZFinderElectron* electron) { e0_truth = electron; }
      void set_e1_truth(ZFinderElectron* electron) { e1_truth = electron; }
      void set_both_e_truth(ZFinderElectron* electron0, ZFinderElectron* electron1) { e0_truth = electron0; e1_truth = electron1; }
      ZFinderElectron* e0_trig;
      ZFinderElectron* e1_trig;
      void set_e0_trig(ZFinderElectron* electron) { e0_trig = electron; }
      void set_e1_trig(ZFinderElectron* electron) { e1_trig = electron; }
      void set_both_e_trig(ZFinderElectron* electron0, ZFinderElectron* electron1) { e0_trig = electron0; e1_trig = electron1; }

      std::vector<const reco::Candidate*> jpsi_muon0;
      std::vector<const reco::Candidate*> jpsi_muon1;

      const reco::Candidate* z_truth_muon0;
      const reco::Candidate* z_truth_muon1;

      const reco::Candidate* z_truth_electron0;
      const reco::Candidate* z_truth_electron1;
      //void set_mu0(ZFinderMuon* muon) { mu0 = muon; }
      //void set_mu1(ZFinderMuon* muon) { mu1 = muon; }
      //void set_both_mu(ZFinderMuon* muon0, ZFinderMuon* muon1) { mu0 = muon0; mu1 = muon1; }
      //void set_both_mu(const reco::Muon &muon0, const reco::Muon &muon1) { mu0 = muon0; mu1 = muon1; }

      // Access pruned lists of the internal electrons
      std::vector<ZFinderElectron*>* FilteredElectrons();
      std::vector<ZFinderElectron*>* AllElectrons() { return FilteredElectrons(); }
      std::vector<ZFinderElectron*>* FilteredElectrons(const std::string& cut_name);

      //std::vector<reco::PFJet> jets;

      // Number of Electrons
      int n_reco_electrons;

      // Number of Anti Electrons
      int n_reco_anti_electrons;

      // Number of Muons
      int n_reco_muons;
      
      //TODO testing jpsi->ee
      int n_reco_jpsi_from_electrons;

      // Number of Jets
      int n_reco_jets;
      int n_reco_muon_jets;

      bool found_four_muons;

      bool found_high_pt_muons_from_z;
      bool found_good_muons_from_z;
      bool found_dimuon_z_compatible_vertex; 
      bool found_z_to_muons_mass;

      bool found_high_pt_electrons_from_z;
      bool found_good_electrons_from_z;
      bool found_dielectron_z_compatible_vertex; 
      bool found_z_to_electrons_mass;

      std::vector<bool> found_dimuon_jpsi_with_soft_id_and_high_pt_muons;
      std::vector<bool> found_dimuon_jpsi_with_good_muons_and_compatible_muon_vertex;
      std::vector<bool> found_good_dimuon_jpsi_compatible_with_primary_vertex;
      bool found_jpsi;


      //TODO jpsi->ee
      bool found_dimuon_jpsi_from_electrons_with_soft_id_and_high_pt_muons;
      bool found_dimuon_jpsi_from_electrons_with_good_muons_and_compatible_muon_vertex;
      bool found_good_dimuon_jpsi_from_electrons_compatible_with_primary_vertex;
      bool found_jpsi_from_electrons;
      //--------------

      bool found_truth_jpsi_with_high_pt_muons;

      //reco::TrackRef GetElectronTrackRef(const reco::GsfElectron & e);
      reco::TrackRef GetMuonTrackRef(const reco::Muon & mu);


      //TODO testing for jpsi->ee
      //reco::TrackRef GetElectronTrackRef(const reco::GsfElectron & el);

      // Output
      void PrintElectrons(const int TYPE = 0, const bool PRINT_CUTS = false);  // 0 is reco, 1 is truth, 2 is trig
      void PrintTruthElectrons(const bool PRINT_CUTS = false) { PrintElectrons(1, PRINT_CUTS); }
      void PrintRecoElectrons(const bool PRINT_CUTS = false) { PrintElectrons(0, PRINT_CUTS); }
      void PrintTrigElectrons(const bool PRINT_CUTS = false) { PrintElectrons(2, PRINT_CUTS); }

      //TODO is this needed? Can it be safely cleaned up for J/Psi part of the code?
      // Access ZDefinition information
      void AddZDef(const std::string NAME, cutlevel_vector PASS_OBJ) { zdef_map_[NAME] = PASS_OBJ; }
      const cutlevel_vector* GetZDef(const std::string& NAME) const;
      bool ZDefPassed(const std::string& NAME) const;
      void PrintZDefs(const bool VERBOSE = false) const;

    protected:
      // These variables are defined at the top of ZFinderEvent.cc to
      // avoid compilation issues
      static const double TRIG_DR_;

      void SetLumiEventWeight(const edm::Event& iEvent);
      void SetMCEventWeight(const edm::Event& iEvent);

      // Called by the constructor to handle MC and Data separately
      void InitReco(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      void InitTruth(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      void InitTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup);

      void InitGSFElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup);

      // Update the Z Info from e0, e1
      void InitZFromElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      void InitZFromMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup);

      // Update the JPsi Info from two muons
      void InitJPsi(const reco::Muon& mu0, const reco::Muon& mu1, const TransientVertex &dimuon_vertex);
      //TODO testing
      void InitJPsiFromElectrons(const reco::GsfElectron& e0, const reco::GsfElectron& e1, const TransientVertex &dielectron_vertex);

      void InitJets(const edm::Event& iEvent, const edm::EventSetup& iSetup) ;

      // Initialize all variables to safe values
      void InitVariables();

      // Input tags
      struct InputTags{
        edm::InputTag ecal_electron;
        edm::InputTag nt_electron;
        edm::InputTag muon;
        edm::InputTag jet;
        edm::InputTag conversion;
        edm::InputTag beamspot;
        edm::InputTag rho_iso;
        edm::InputTag vertex;
        edm::InputTag pileup;
        edm::InputTag generator;
        std::vector<edm::InputTag> iso_vals;
      } inputtags_;

      // Find matching trigger objects
      const trig_dr_vec* GetMatchedTriggerObjects(
          const edm::Event& iEvent,
          const std::vector<std::string>& trig_names,
          const double ETA, const double PHI, const double DR_CUT
          );
      const trigger::TriggerObject* GetBestMatchedTriggerObject(
          const edm::Event& iEvent,
          const std::vector<std::string>& trig_names,
          const double ETA, const double PHI
          );
      bool TriggerMatch(
          const edm::Event& iEvent,
          const std::vector<std::string>& trig_names,
          const double ETA, const double PHI, const double DR_CUT
          );

      double JpsiMuonTruthMatch(const reco::Muon&);
      double JpsiMuonZMuonMatch(const reco::Muon&);
      double sumPtSquared(const reco::Vertex&);

      // A list of all electrons, split into reco and gen
      std::vector<ZFinderElectron*> reco_electrons_;
      std::vector<ZFinderElectron*> reco_anti_electrons_;
      ZFinderElectron* AddRecoElectron(reco::GsfElectron electron);
      void AddRecoElectron(zf::ZFinderElectron zf_electron);
      ZFinderElectron* AddRecoElectron(reco::RecoEcalCandidate electron);
      ZFinderElectron* AddRecoElectron(reco::Photon electron);

      std::vector<ZFinderElectron*> truth_electrons_;
      ZFinderElectron* AddTruthElectron(reco::GenParticle electron);

      std::vector<ZFinderElectron*> hlt_electrons_;
      ZFinderElectron* AddHLTElectron(trigger::TriggerObject electron);

      // Calculate phistar
      static double ReturnPhistar(const double& eta0, const double& phi0, const double& eta1, const double& phi1);

      // Sorting functions
      static bool SortByPTHighLowElectron(const ZFinderElectron* e0, const ZFinderElectron* e1) { return (e0->pt > e1->pt); }

      // Print cuts
      void PrintCuts(ZFinderElectron* zf_elec);

      // Store ZDefinition Information
      std::map<std::string, cutlevel_vector> zdef_map_;
      // Pileup reweighting
      static edm::LumiReWeighting* lumi_weights_;

  };
}  // namespace zf
#endif  // ZFINDER_ZFINDEREVENT_H_
