#ifndef ZFINDER_ZFINDERPLOTTER_H_
#define ZFINDER_ZFINDERPLOTTER_H_

// Root
#include <TH1D.h>  // TH1D
#include <TH2.h>  // TH2D

// CMSSW
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/MuonReco/interface/Muon.h" // reco::Muon
#include "DataFormats/MuonReco/interface/MuonSelectors.h" //muon::isSoftMuon

// ZFinder Code
#include "ZFinderEvent.h"  // ZFinderEvent


namespace zf {
    class ZFinderPlotter{
        public:
            // Constructor
            ZFinderPlotter(TFileDirectory& tdir, const bool USE_MC = false, const bool APPLY_MUON_MIN_PT = false, const bool APPLY_SOFT_MUONS = false,
                const bool APPLY_DIMUON_VTX_COMPATIBILITY = false, const bool APPLY_JPSI_MASS_WINDOW = false, const bool APPLY_VERTEX_Z_POS_WINDOW = false,
                const bool APPLY_PROMPT_JPSI_WINDOW = false);

            // Add events
            void Fill(
                    const ZFinderEvent& zf_event,
                    const double EVENT_WEIGHT = 1.
                    );
            //void CalculateJpsiLifetime( const ZFinderEvent::JPsiData& , const ZFinderEvent::ZFromMuonsData& ); 

        protected:
            // Histograms
            TH1D* z_mass_all_;
            TH1D* z_mass_coarse_;
            TH1D* z_mass_fine_;
            TH1D* z_rapidity_;
            TH1D* z_pt_;
            TH1D* z_vtx_prob_;
            TH1D* phistar_;

            TH1D* z_from_muons_mass_all_;
            TH1D* z_from_muons_mass_coarse_;
            TH1D* z_from_muons_mass_fine_;
            TH1D* z_from_muons_rapidity_;
            TH1D* z_from_muons_pt_;
            TH1D* z_from_muons_vtx_prob_;
            TH1D* z_from_muons_phistar_;

            TH1D* muon0_from_z_pt_;
            TH1D* muon1_from_z_pt_;
            TH1D* muon0_from_z_eta_;
            TH1D* muon1_from_z_eta_;
            TH1D* muon0_from_z_phi_;
            TH1D* muon1_from_z_phi_;
            TH1D* muon0_from_z_charge_;
            TH1D* muon1_from_z_charge_;

            TH1D* e0_pt_;
            TH1D* e1_pt_;
            TH1D* e0_eta_;
            TH1D* e1_eta_;
            TH1D* e0_phi_;
            TH1D* e1_phi_;
            TH1D* e0_charge_;
            TH1D* e1_charge_;

            TH1D* baseweights_;
            TH1D* pileup_;
            TH1D* nelectrons_;

            TH1D* jpsi_mass_all_;


            TH1D* jpsi_mass_coarse_;

            TH1D* jpsi_mass_fine_;

            //TDOO jpsi->ee testing
            //---------------------
            TH1D* jpsi_from_electrons_mass_fine_;
            TH1D* jpsi_from_electrons_pt_;
            TH1D* jpsi_from_electrons_tau_xy_very_fine_;
            //------------------------

            TH1D* jpsi_mass_fine_ptUnder10_;
            TH1D* jpsi_mass_fine_pt10to15_;
            TH1D* jpsi_mass_fine_pt15to20_;
            TH1D* jpsi_mass_fine_pt20to25_;
            TH1D* jpsi_mass_fine_pt25to30_;
            TH1D* jpsi_mass_fine_ptAbove30_;

            TH1D* jpsi_four_lepton_mass_;

            TH1D* jpsi_rapidity_;
            TH1D* jpsi_pt_;


            TH2D* jpsi_pt_vs_rap_;
            TH2D* jpsi_pt_vs_rap_fine_;
            TH2D* jpsi_pt_vs_rap_polarization_long;
            TH2D* jpsi_pt_vs_rap_polarization_TPlusZero;
            TH1D* jpsi_efficiency_;
            TH1D* jpsi_acc_eff_;
            TH1D* jpsi_scale_factor_;
            TH2D* jpsi_reco_pt_vs_jpsi_truth_pt_;
            TH1D* jpsi_truth_pt_minus_jpsi_reco_pt_;
            TH1D* jpsi_trigger_obj_mu0_pt_;
            TH1D* jpsi_trigger_obj_mu1_pt_;
            TH1D* jpsi_trigger_obj_mu0_pt_minus_reco_mu0_pt_;
            TH1D* jpsi_trigger_obj_mu1_pt_minus_reco_mu1_pt_;

            TH1D* jpsi_cos_mu_plus_;
            TH1D* jpsi_cos_mu_plus_jpsi_pt_8to8p5_;
            TH1D* jpsi_cos_mu_plus_jpsi_pt_8p5to9_;
            TH1D* jpsi_cos_mu_plus_jpsi_pt_9to10_;
            TH1D* jpsi_cos_mu_plus_jpsi_pt_10to15_;
            TH1D* jpsi_cos_mu_plus_jpsi_pt_15to20_;
            TH1D* jpsi_cos_mu_plus_lambdaNeg1_;
            TH1D* jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8to8p5_;
            TH1D* jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8p5to9_;
            TH1D* jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_9to10_;
            TH1D* jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_10to15_;
            TH1D* jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_15to20_;
            TH1D* jpsi_cos_mu_plus_lambda1_;
            TH1D* jpsi_cos_mu_plus_lambda1_jpsi_pt_8to8p5_;
            TH1D* jpsi_cos_mu_plus_lambda1_jpsi_pt_8p5to9_;
            TH1D* jpsi_cos_mu_plus_lambda1_jpsi_pt_9to10_;
            TH1D* jpsi_cos_mu_plus_lambda1_jpsi_pt_10to15_;
            TH1D* jpsi_cos_mu_plus_lambda1_jpsi_pt_15to20_;
            TH1D* jpsi_cos_mu_minus_;

            TH1D* mu0_pt_;
            TH1D* mu1_pt_;
            TH1D* mu0_eta_;
            TH1D* mu1_eta_;
            TH1D* mu0_efficiency_;
            TH1D* mu1_efficiency_;
            TH1D* mu0_scale_factor_;
            TH1D* mu1_scale_factor_;
            TH1D* mu0_phi_;
            TH1D* mu1_phi_;
            TH1D* mu0_charge_;
            TH1D* mu1_charge_;
            TH1D* mu0_tracker_layers_;
            TH1D* mu1_tracker_layers_;
            TH1D* mu0_deltaR_truth_;
            TH1D* mu1_deltaR_truth_;
            TH1D* mu0_deltaR_z_muon_;
            TH1D* mu1_deltaR_z_muon_;
            TH1D* n_truth_matched_jpsi_muons_;

            TH1D* jpsi_distance_;
            TH1D* jpsi_dist_err_;
            TH1D* jpsi_chi2_;
            TH1D* jpsi_distance_xy_;
            TH1D* jpsi_dist_err_xy_;
            TH1D* jpsi_chi2_xy_;
            TH1D* jpsi_tau_xy_;
            TH1D* jpsi_tau_xy_fine_;
            TH1D* jpsi_tau_xy_very_fine_;

            TH1D* jpsi_tau_xy_very_fine_ptUnder10_;
            TH1D* jpsi_tau_xy_very_fine_pt10to15_;
            TH1D* jpsi_tau_xy_very_fine_pt15to20_;
            TH1D* jpsi_tau_xy_very_fine_pt20to25_;
            TH1D* jpsi_tau_xy_very_fine_pt25to30_;
            TH1D* jpsi_tau_xy_very_fine_ptAbove30_;

            TH1D* jpsi_tau_xy_very_fine_rap0_0to0_3_;
            TH1D* jpsi_tau_xy_very_fine_rap0_3to0_6_;
            TH1D* jpsi_tau_xy_very_fine_rap0_6to0_9_;
            TH1D* jpsi_tau_xy_very_fine_rap0_9to1_2_;
            TH1D* jpsi_tau_xy_very_fine_rap1_2to1_5_;
            TH1D* jpsi_tau_xy_very_fine_rap1_5to1_8_;
            TH1D* jpsi_tau_xy_very_fine_rap1_8to2_1_;
            TH1D* jpsi_tau_xy_very_fine_rap2_1to2_4_;

            TH1D* jpsi_tau_xy_very_fine_rap0_0to0_9pt10to15_;
            TH1D* jpsi_tau_xy_very_fine_rap0_9to1_2pt10to15_;
            TH1D* jpsi_tau_xy_very_fine_rap1_2to2_1pt10to15_;

            TH1D* jpsi_tau_xy_very_fine_above_12_tracker_layers_;

            TH1D* jpsi_tau_xy_very_fine_similar_pt_muons_;

            TH1D* jpsi_tau_xy_dimuon_continuum_bg_;

            TH1D* jpsi_tau_z_;
            TH1D* jpsi_tau_z_fine_;
            TH1D* jpsi_tau_z_very_fine_;
            TH1D* jpsi_zpt_difference_;
            TH1D* jpsi_zmumupt_difference_;

            TH1D* dimuon_vtx_prob_;
            TH1D* dimuon_delta_phi_;
            TH1D* dimuon_delta_eta_;
            TH1D* dimuon_deltaR_;

            TH1D* z_jpsi_delta_phi_;

            TH1D* jet_pt_;
            TH1D* jet_eta_;
            TH1D* jet_btag_discriminator_;
            TH1D* muon_jet_pt_;
            TH1D* muon_jet_pt_diff_z_pt_;
            TH1D* muon_jet_pt_diff_dimuon_pt_;
            TH1D* muon_jet_eta_;
            TH1D* muon_jet_btag_discriminator_;

            TH2D* muon_jet_pt_z_pt_;
            TH2D* muon_jet_pt_dimuon_pt_;
            TH2D* muon_jet_phi_z_phi_;
            TH2D* muon_jet_phi_dimuon_phi_;

            TH1D* jpsi_vtx_distance_z_vtx_x_;
            TH1D* jpsi_vtx_distance_z_vtx_y_;
            TH1D* jpsi_vtx_distance_z_vtx_z_;

            TH1D* jpsi_iso_mu0_;
            TH1D* jpsi_iso_sum_charged_hadron_pt_mu0_;
            TH1D* jpsi_iso_sum_charged_particle_pt_mu0_;
            TH1D* jpsi_iso_sum_neutral_hadron_et_mu0_;
            TH1D* jpsi_iso_sum_photon_et_mu0_;
            TH1D* jpsi_iso_sum_pileup_pt_mu0_;

            TH1D* jpsi_iso_mu1_;
            TH1D* jpsi_iso_sum_charged_hadron_pt_mu1_;
            TH1D* jpsi_iso_sum_charged_particle_pt_mu1_;
            TH1D* jpsi_iso_sum_neutral_hadron_et_mu1_;
            TH1D* jpsi_iso_sum_photon_et_mu1_;
            TH1D* jpsi_iso_sum_pileup_pt_mu1_;

            TH2D* jpsi_mass_vs_chi2_;
            TH2D* jpsi_tau_xy_vs_tau_z_;
            TH2D* jpsi_tau_xy_vs_distance_z_;
            TH2D* jpsi_tau_z_vs_distance_z_;
            TH2D* dimuon_mass_vs_dimuon_tau_xy_;

            TH2D* dimuon_pt_vs_zee_pt_;
            TH2D* dimuon_pt_vs_zmumu_pt_;

            TH2D* dimuon_mass_vs_dimuon_tau_xy_fine_;
            TH2D* dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_;
            TH2D* dimuon_mass_vs_dimuon_tau_xy_10to14_;
            TH2D* dimuon_mass_vs_dimuon_tau_xy_14to18_;
            TH2D* dimuon_mass_vs_dimuon_tau_xy_18to30_;
            TH2D* dimuon_mass_vs_dimuon_tau_xy_30to100_;

            TH2D* low_rap_dimuon_mass_vs_dimuon_tau_xy_fine_;
            TH2D* low_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_;
            TH2D* low_rap_dimuon_mass_vs_dimuon_tau_xy_10to14_;
            TH2D* low_rap_dimuon_mass_vs_dimuon_tau_xy_14to18_;
            TH2D* low_rap_dimuon_mass_vs_dimuon_tau_xy_18to30_;
            TH2D* low_rap_dimuon_mass_vs_dimuon_tau_xy_30to100_;

            TH2D* high_rap_dimuon_mass_vs_dimuon_tau_xy_fine_;
            TH2D* high_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_;
            TH2D* high_rap_dimuon_mass_vs_dimuon_tau_xy_10to14_;
            TH2D* high_rap_dimuon_mass_vs_dimuon_tau_xy_14to18_;
            TH2D* high_rap_dimuon_mass_vs_dimuon_tau_xy_18to30_;
            TH2D* high_rap_dimuon_mass_vs_dimuon_tau_xy_30to100_;



            TH2D* dimuon_mass_vs_dimuon_tau_xy_fine_weighted_;

            TH1D* vtx_x_;
            TH1D* vtx_y_;
            TH1D* vtx_z_;

            TH1D* primary_vtx_x_;
            TH1D* primary_vtx_y_;
            TH1D* primary_vtx_z_;

            TH1D* z_vtx_x_;
            TH1D* z_vtx_y_;
            TH1D* z_vtx_z_;

            TH1D* z_from_muons_vtx_x_;
            TH1D* z_from_muons_vtx_y_;
            TH1D* z_from_muons_vtx_z_;

            TH1D* dimuon_vtx_x_;
            TH1D* dimuon_vtx_y_;
            TH1D* dimuon_vtx_z_;

            TH1D* jpsi_truth_vtx_x_minus_jpsi_reco_vtx_x_;
            TH1D* jpsi_truth_vtx_y_minus_jpsi_reco_vtx_y_;
            TH1D* jpsi_truth_vtx_z_minus_jpsi_reco_vtx_z_;

            TH1D* primary_vtx_x_minus_z_vtx_x_;
            TH1D* primary_vtx_y_minus_z_vtx_y_;
            TH1D* primary_vtx_z_minus_z_vtx_z_;

            TH1D* primary_vtx_x_minus_zmumu_vtx_x_;
            TH1D* primary_vtx_y_minus_zmumu_vtx_y_;
            TH1D* primary_vtx_z_minus_zmumu_vtx_z_;

            TH2D* primary_vtx_x_vs_z_vtx_x_;
            TH2D* primary_vtx_y_vs_z_vtx_y_;
            TH2D* primary_vtx_z_vs_z_vtx_z_;

            TH2D* primary_vtx_x_vs_zmuons_vtx_x_;
            TH2D* primary_vtx_y_vs_zmuons_vtx_y_;
            TH2D* primary_vtx_z_vs_zmuons_vtx_z_;

            TH1D* nmuons_;
            TH1D* njets_;
            TH1D* n_muonjets_;
            TH1D* njpsis_;

            // Use the MC or reco data
            const bool USE_MC_;
            const bool APPLY_MUON_MIN_PT_;
            const bool APPLY_SOFT_MUONS_;
            const bool APPLY_DIMUON_VTX_COMPATIBILITY_;
            const bool APPLY_JPSI_MASS_WINDOW_;
            const bool APPLY_VERTEX_Z_POS_WINDOW_;
            const bool APPLY_PROMPT_JPSI_WINDOW_;
            // Plotting variables
            static const int X_SIZE = 1280;
            static const int Y_SIZE = 640;
    };
}  // namespace zf
#endif  // ZFINDER_ZFINDERPLOTTER_H_
