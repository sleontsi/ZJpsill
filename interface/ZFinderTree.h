#ifndef ZFINDER_ZFINDERTREE_H_
#define ZFINDER_ZFINDERTREE_H_

// Standard Library
#include <string>  // string
#include <utility>  // pair

// ROOT
#include "TBranch.h"  // TBranch
#include "TDirectory.h"  // TDirectory
#include "TTree.h"  // TTree
#include <TLorentzVector.h>

// CMSSW
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// ZFinder Code
#include "ZFinderEvent.h"  // ZFinderEvent

namespace zf {
  class ZFinderTree {
    public:
      // Constructor
      ZFinderTree(
          TFileDirectory& tdir,
          const bool IS_MC = false
          );

      // destructor
      ~ZFinderTree();

      // Add event
      void Fill(const ZFinderEvent& zf_event);

      // Wrapper around TTree::GetCurrentFile()
      TFile* GetCurrentFile();

    protected:
      double ZMassPDG = 91.1876;

      // Structs that map to the branches
      struct jpsi_branch {
        void clear_values() {
          jpsi_m        .clear();
          jpsi_pt       .clear();
          jpsi_y        .clear();
          jpsi_phi      .clear();
          jpsi_eta      .clear();
          jpsi_vtx_prob .clear();
          jpsi_vtx_x    .clear();
          jpsi_vtx_y    .clear();
          jpsi_vtx_z    .clear();

          jpsi_tau_xy      .clear();
          jpsi_tau_z       .clear();
          jpsi_distance_xy .clear();
          jpsi_distance_z  .clear();

          jpsi_eff         .clear();
          jpsi_acc_eff     .clear();
          jpsi_scale_factor.clear();

          muon0_pt         .clear();
          muon0_eta        .clear();
          muon0_phi        .clear();
          muon0_d0         .clear();
          muon0_dxy        .clear();
          muon0_dz         .clear();
          muon0_d0err      .clear();
          muon0_dxyerr     .clear();
          muon0_dzerr      .clear();
          muon0_trkKink    .clear(); 
          muon0_glbKink    .clear(); 

          muon1_pt         .clear();
          muon1_eta        .clear();
          muon1_phi        .clear();
          muon1_d0         .clear();
          muon1_dxy        .clear();
          muon1_dz         .clear();
          muon1_d0err      .clear();
          muon1_dxyerr     .clear();
          muon1_dzerr      .clear();
          muon1_trkKink    .clear(); 
          muon1_glbKink    .clear(); 

          muon0_charge     .clear();
          muon1_charge     .clear();

          reco_jpsi_soft_id.clear();
        }
        // Constructor
        jpsi_branch() {
          clear_values();
        }
        std::vector<double> jpsi_m;
        std::vector<double> jpsi_pt;
        std::vector<double> jpsi_y;
        std::vector<double> jpsi_phi;
        std::vector<double> jpsi_eta;
        std::vector<double> jpsi_vtx_prob;
        std::vector<double> jpsi_vtx_x;
        std::vector<double> jpsi_vtx_y;
        std::vector<double> jpsi_vtx_z;
        std::vector<double> jpsi_tau_xy;
        std::vector<double> jpsi_tau_z;
        std::vector<double> jpsi_distance_xy;
        std::vector<double> jpsi_distance_z;
        std::vector<double> jpsi_eff;
        std::vector<double> jpsi_acc_eff;
        std::vector<double> jpsi_scale_factor;
        std::vector<double> muon0_pt;
        std::vector<double> muon0_eta;
        std::vector<double> muon0_phi;
        std::vector<double> muon0_d0;
        std::vector<double> muon0_dxy;
        std::vector<double> muon0_dz;
        std::vector<double> muon0_d0err;
        std::vector<double> muon0_dxyerr;
        std::vector<double> muon0_dzerr;
        std::vector<double> muon0_trkKink;
        std::vector<double> muon0_glbKink;
        std::vector<double> muon1_pt;
        std::vector<double> muon1_eta;
        std::vector<double> muon1_phi;
        std::vector<double> muon1_d0;
        std::vector<double> muon1_dxy;
        std::vector<double> muon1_dz;
        std::vector<double> muon1_d0err;
        std::vector<double> muon1_dxyerr;
        std::vector<double> muon1_dzerr;
        std::vector<double> muon1_trkKink;
        std::vector<double> muon1_glbKink;
        std::vector<double> reco_jpsi_soft_id;
        std::vector<double> muon0_iso_sum_charged_hadron_pt;                                     
        std::vector<double> muon0_iso_sum_charged_particle_pt;
        std::vector<double> muon0_iso_sum_neutral_hadron_et;
        std::vector<double> muon0_iso_sum_photon_et;
        std::vector<double> muon0_iso_sum_pileup_pt;
        std::vector<double> muon0_iso;
        std::vector<double> muon1_iso_sum_charged_hadron_pt;                                     
        std::vector<double> muon1_iso_sum_charged_particle_pt;
        std::vector<double> muon1_iso_sum_neutral_hadron_et;
        std::vector<double> muon1_iso_sum_photon_et;
        std::vector<double> muon1_iso_sum_pileup_pt;
        std::vector<double> muon1_iso;
        std::vector<int>    muon0_charge;
        std::vector<int>    muon1_charge;
        std::vector<int>    has_muons_in_eta_window;
        std::vector<int>    has_high_pt_muons;
      } reco_jpsi_, truth_jpsi_, reco_jpsi_from_electrons_;

      struct four_lepton_branch {
        void clear_values() {
          muon0_pt .clear();
          muon1_pt .clear();
          muon2_pt .clear();
          muon3_pt .clear();
          muon0_eta.clear();
          muon1_eta.clear();
          muon2_eta.clear();
          muon3_eta.clear();
          muon0_phi.clear();
          muon1_phi.clear();
          muon2_phi.clear();
          muon3_phi.clear();
          vtx_chi2.clear();
          vtx_ndf.clear();
          vtx_prob.clear();
        }
        // Constructor
        four_lepton_branch() {
          clear_values();
        }
        std::vector<double> muon0_pt; 
        std::vector<double> muon1_pt;
        std::vector<double> muon2_pt; 
        std::vector<double> muon3_pt;
        std::vector<double> muon0_eta;
        std::vector<double> muon1_eta;
        std::vector<double> muon2_eta;
        std::vector<double> muon3_eta;
        std::vector<double> muon0_phi;
        std::vector<double> muon1_phi;
        std::vector<double> muon2_phi;
        std::vector<double> muon3_phi;
        std::vector<double> vtx_chi2;
        std::vector<double> vtx_ndf;
        std::vector<double> vtx_prob;
      } four_lepton_vertex_;

      struct z_branch {
        void clear_values() {
          z_m        .clear(); //= -1000;
          z_pt       .clear(); //= -1000;
          z_y        .clear(); //= -1000;
          z_phi      .clear(); //= -1000;
          z_phistar  .clear(); //= -1000;
          z_eta      .clear(); //= -1000;
          z_vtx_prob .clear(); //= -1000;
          z_vtx_x    .clear(); //= -1000;
          z_vtx_y    .clear(); //= -1000;
          z_vtx_z    .clear(); //= -1000;

          daughter0_pt     .clear();//= -1000;
          daughter0_eta    .clear();//= -1000;
          daughter0_phi    .clear();//= -1000;
          daughter0_d0     .clear();
          daughter0_dxy    .clear();
          daughter0_dz     .clear();
          daughter0_d0err  .clear();
          daughter0_dxyerr .clear();
          daughter0_dzerr  .clear();

          daughter0_trkKink.clear(); 
          daughter0_glbKink.clear(); 

          daughter1_pt     .clear();//= -1000;
          daughter1_eta    .clear();//= -1000;
          daughter1_phi    .clear();//= -1000;
          daughter1_d0     .clear();
          daughter1_dxy    .clear();
          daughter1_dz     .clear();
          daughter1_d0err  .clear();
          daughter1_dxyerr .clear();
          daughter1_dzerr  .clear();
          daughter1_trkKink.clear(); 
          daughter1_glbKink.clear(); 

          daughter0_charge .clear();//= 0;
          daughter1_charge .clear();//= 0;

          passed_triggers.clear();
          truth_pdg_id.clear();
        }
        // Constructor
        z_branch() {
          clear_values();
        }
        std::vector<double> z_m;
        std::vector<double> z_pt;
        std::vector<double> z_y;
        std::vector<double> z_phi;
        std::vector<double> z_phistar;
        std::vector<double> z_eta;
        std::vector<double> z_vtx_prob;
        std::vector<double> z_vtx_x;
        std::vector<double> z_vtx_y;
        std::vector<double> z_vtx_z;
        std::vector<double> daughter0_pt;
        std::vector<double> daughter0_eta;
        std::vector<double> daughter0_phi;
        std::vector<double> daughter0_d0;
        std::vector<double> daughter0_dxy;
        std::vector<double> daughter0_dz;
        std::vector<double> daughter0_d0err;
        std::vector<double> daughter0_dxyerr;
        std::vector<double> daughter0_dzerr;
        std::vector<double> daughter0_trkKink;
        std::vector<double> daughter0_glbKink;
        std::vector<double> daughter1_pt;
        std::vector<double> daughter1_eta;
        std::vector<double> daughter1_phi;
        std::vector<double> daughter1_d0;
        std::vector<double> daughter1_dxy;
        std::vector<double> daughter1_dz;
        std::vector<double> daughter1_d0err;
        std::vector<double> daughter1_dxyerr;
        std::vector<double> daughter1_dzerr;
        std::vector<double> daughter1_trkKink;
        std::vector<double> daughter1_glbKink;
        std::vector<double> daughter0_charge;
        std::vector<double> daughter1_charge;
        std::vector<int>  truth_pdg_id;
        std::vector<std::string> passed_triggers;
      } reco_z_, reco_z_from_muons_, truth_z_muons_,truth_z_electrons_;

      struct event_branch {
        void clear_values() {
          event_weight                                                                .clear();//= -1.0;
          event_number                                                                .clear();//= 0;
          run_number                                                                  .clear();//= 0;
          n_verts                                                                     .clear();//= 0;
          found_good_muons_from_z                                                     .clear();//= false;
          found_good_electrons_from_z                                                 .clear();//= false;
        }

        // Constructor
        event_branch() {
          clear_values();
        }
        std::vector<double>        event_weight;
        std::vector<unsigned int>  event_number;
        std::vector<unsigned int>  run_number;
        std::vector<int>           n_verts;
        std::vector<bool>          found_good_muons_from_z;
        std::vector<bool>          found_good_electrons_from_z;

      } event_;

      // File Directory to write to
      TDirectory* tdir_;

      // Use the MC or reco data
      const bool IS_MC_;

      // The tuples
      TTree* tree_;

      //// Set up a variable size branch for the weights
      //int weight_size_;
      //double weight_fsr_;
      //int weight_cteq_size_;
      //int weight_mstw_size_;
      //int weight_nnpdf_size_;
      //static constexpr int MAX_SIZE_ = 100;
      //static constexpr int MAX_SIZE_PDF_ = 110;
      //
      //// Although vectors seem like the right solution, since TTrees need
      //// the memory used for the array to be static, an array is
      //// (unfortunately) the best choice
      //double weights_[MAX_SIZE_];
      //int weight_ids_[MAX_SIZE_];
      //double weights_cteq_[MAX_SIZE_PDF_];
      //double weights_mstw_[MAX_SIZE_PDF_];
      //double weights_nnpdf_[MAX_SIZE_PDF_];
      //
      //// We insert the weights and the IDs into this vector, and then
      //// read it out into the array before filling the tree
      //std::vector<std::pair<int, double>> weight_id_vector_;
      //
      //// Get the weight of the cuts
      //void FillCutWeights(cutlevel_vector const * const CUT_LEVEL_VECTOR);
      double GetTotalWeight(cutlevel_vector const * const CUT_LEVEL_VECTOR);
  };
}  // namespace zf
#endif  // ZFINDER_ZFINDERTREE_H_
