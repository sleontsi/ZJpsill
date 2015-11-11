#include "ZFinder/Event/interface/MuonEfficiency.h"
#include "ZFinder/Event/interface/JpsiEfficiencyTables.h"
#include <math.h>

namespace zf {
  double GetEfficiency (const double EFF_TABLE[][7], 
      const double ETA,
      const double PT ) {
    //TODO to modify this, need GetEfficiency to take number of rows of the table perhaps?
    //alternatively, may want to use vectors
    int eta_min = 0;
    int eta_max = 1;
    int pt_min = 2;
    int pt_max = 3;
    int eff = 4;
    //int eff_err_low = 5;
    //int eff_err_hi = 6;

    double efficiency = -1000;
    for ( int i=0; i < zf::SOFT_MUON_DATA_EFF_TABLE_ROWS; ++i) {
      if (fabs(ETA) >= EFF_TABLE[i][eta_min] && fabs(ETA) < EFF_TABLE[i][eta_max]  && 
        (PT >= EFF_TABLE[i][pt_min] && PT < EFF_TABLE[i][pt_max] )) {
          efficiency =  EFF_TABLE[i][eff];
        }
    }
    return efficiency;
  }
  //TODO comine this with previous function, kind of clunky
  double GetAccEff (const double EFF_TABLE[][7], 
      const double ETA,
      const double PT ) {
    //TODO to modify this, need GetEfficiency to take number of rows of the table perhaps?
    //alternatively, may want to use vectors
    int eta_min = 0;
    int eta_max = 1;
    int pt_min = 2;
    int pt_max = 3;
    int eff = 4;
    //int eff_err_low = 5;
    //int eff_err_hi = 6;

    double acc_eff = -1000;
    for ( int i=0; i < zf::SOFT_MUON_DATA_ACC_EFF_TABLE_ROWS; ++i) {
      //note should not use fabs here, perhaps should combine bins or something though
      if (ETA >= EFF_TABLE[i][eta_min] && ETA < EFF_TABLE[i][eta_max]  && 
        (PT >= EFF_TABLE[i][pt_min] && PT < EFF_TABLE[i][pt_max] )) {
          acc_eff =  EFF_TABLE[i][eff];
        }
    }
    return acc_eff;
  }
} //namespace zf
