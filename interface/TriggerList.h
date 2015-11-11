#ifndef ZFINDER_TRIGGERLIST_H_
#define ZFINDER_TRIGGERLIST_H_

/* These vectors contain all of the trigger filters that make up the triggers
 * we use. */

namespace zf {
  static const std::vector<std::string> ET_ET_TIGHT = {
    "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter"  // ET-ET (ET == ECAL TOTAL) Tight Leg
  };
  static const std::vector<std::string> ET_ET_DZ = {
    "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ"  // ET-ET Vertex distance
  };
  static const std::vector<std::string> ET_ET_LOOSE = {
    "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter"  // ET-ET Loose Leg
  };
  static const std::vector<std::string> ET_NT_ET_TIGHT = {
    "hltEle27TightIdLooseIsoTrackIsoFilter"  // ET-NT Tight ECAL Leg AND ET-HF TightLoose ECAL Leg
  };

  static const std::vector<std::string> ET_HF_ET_TIGHT = ET_NT_ET_TIGHT;

  static const std::vector<std::string> ET_HF_ET_LOOSE = {
    "hltEle23TightIdLooseIsoTrackIsoFilter"  // ET-HF LooseTight ECAL Leg
  };

  static const std::vector<std::string> ET_HF_HF_TIGHT = {
    "hltHFEMPt30TightFilter" // ET-HF LooseTight HF Leg (30pt)
  };

  static const std::vector<std::string> ET_HF_HF_LOOSE = {
    "hltHFEMTightFilter"  // ET-HF TightLoose HF Leg (15pt)
  };

  static const std::vector<std::string> SINGLE_ELECTRON_TRIGGER = {
    "hltEle27WP80TrackIsoFilter"  // Single Electron Trigger
  };

  // The giant combined list
  static const std::vector<std::string> ALL_TRIGGERS = {
    "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ",
    "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter",
    "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",
    "hltEle23TightIdLooseIsoTrackIsoFilter",
    "hltEle27TightIdLooseIsoTrackIsoFilter",
    "hltHFEMPt30TightFilter",
    "hltHFEMTightFilter",
    "hltEle27WP80TrackIsoFilter"
  };
  static const std::vector<std::string> DOUBLE_MUON_TIGHT_LEG_TRIGGER = {
    //"hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8"
    //"hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered17"
    //"hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17"
    "hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17",
    "hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17"
  };
  //hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8
  static const std::vector<std::string> DOUBLE_MUON_LOOSE_LEG_TRIGGER = {
    //"hltL3pfL1DoubleMu10MuOpenL1f0L2f10L3Filtered17"
    //"hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered8"
    //"hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered13"
    //"hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8"
    "hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8",
    "hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8"
  };

  //TODO testing Jpsi trigger
  static const std::vector<std::string> JPSI_TRIGGER = {
    "hltVertexmumuFilterJpsi" //HLT_DIMUON0_JPSI ???
  };

  static const std::vector<std::string> DOUBLE_MUON_TRIGGER = {
    "hltDiMuonGlb17Glb8DzFiltered0p2" //HLT_DIMUON0_JPSI ???
  };

}  // namespace zf
#endif  // ZFINDER_TRIGGERLIST_H_
