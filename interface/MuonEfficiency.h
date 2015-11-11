#ifndef ZFINDER_MUONEFFICIENCY_H_
#define ZFINDER_MUONEFFICIENCY_H_

namespace zf {
      double GetEfficiency (const double EFF_TABLE[][7], 
          const double ETA,
          const double PT ) ;
      double GetAccEff (const double EFF_TABLE[][7], 
          const double ETA,
          const double PT ) ;
}
#endif // ZFINDER_MUONEFFICIENCY_H_
