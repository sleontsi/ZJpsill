#include "ZFinder/Event/interface/ZFinderPlotter.h"

// Standard Library
#include <string>  // string

// Root
#include <TCanvas.h>  // TCanvas

// ZFinder Code
#include "ZFinder/Event/interface/ZFinderElectron.h"  // ZFinderElectron
#include "ZFinder/Event/interface/ZFinderCuts.h"  // ZFinderCuts


namespace zf {
  // Constructor
  ZFinderPlotter::ZFinderPlotter(TFileDirectory& tdir, const bool USE_MC, const bool APPLY_MUON_MIN_PT, const bool APPLY_SOFT_MUONS, const bool APPLY_DIMUON_VTX_COMPATIBILITY, 
      const bool APPLY_JPSI_MASS_WINDOW , const bool APPLY_VERTEX_Z_POS_WINDOW, const bool APPLY_PROMPT_JPSI_WINDOW)
    : USE_MC_(USE_MC), APPLY_MUON_MIN_PT_(APPLY_MUON_MIN_PT), APPLY_SOFT_MUONS_(APPLY_SOFT_MUONS), APPLY_DIMUON_VTX_COMPATIBILITY_(APPLY_DIMUON_VTX_COMPATIBILITY),
    APPLY_JPSI_MASS_WINDOW_(APPLY_JPSI_MASS_WINDOW), APPLY_VERTEX_Z_POS_WINDOW_(APPLY_VERTEX_Z_POS_WINDOW), APPLY_PROMPT_JPSI_WINDOW_(APPLY_PROMPT_JPSI_WINDOW) {
      /*
       * Initialize a set of histograms and associate them with a given TDirectory.
       */

    }

  void ZFinderPlotter::Fill(
      const ZFinderEvent& zfe,
      const double EVENT_WEIGHT
      ) {
    /*
     * Given a zfe, fills all the histograms.
     */
    // Z Info
    double event_weight = zfe.event_weight;
    if (!USE_MC_) {
      for (unsigned int i = 0; i < zfe.reco_jpsi.m.size() ; ++i ) {
        if (!zfe.is_real_data && APPLY_SOFT_MUONS_ && zfe.reco_jpsi.has_soft_id_muons.at(i) &&
            zfe.reco_jpsi.has_high_pt_muons.at(i) && zfe.reco_jpsi.has_muons_in_eta_window.at(i) ) {
          event_weight = event_weight * zfe.reco_jpsi.jpsi_scale_factor.at(i);
          //TODO decide how to handle multiple j/psi case!!
          break;
        }
      }


    }
  }
  //void ZFinderPlotter::CalculateJpsiLifetime(const ZFinderEvent::JPsiData &jpsi_data, const ZFinderEvent::ZFromMuonsData &z_from_muon ) {
  //  for (unsigned int i = 0; i < jpsi_data.m.size() ; ++i ) {
  //  }
  //}
}  // namespace zf
