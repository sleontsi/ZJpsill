#include "ZFinder/Event/interface/ZFinderTree.h"

// Standard Library
#include <algorithm>
#include <vector>  // std::min


#include <typeinfo>


// ZFinder Code


namespace zf {
  // Constructor
  //ZFinderTree::ZFinderTree(const ZDefinition& zdef, TFileDirectory& tdir, const bool IS_MC) : IS_MC_(IS_MC) {
  ZFinderTree::ZFinderTree(TFileDirectory& tdir, const bool IS_MC) : IS_MC_(IS_MC) {

    // Make the Tree to write to
    //tree_ = new TTree(zdef.NAME.c_str(), zdef.NAME.c_str());
    tree_ = new TTree("zfinder_tree", "zfinder_tree");

    tree_->Branch("reco_z_z_m",                 &reco_z_.z_m);
    tree_->Branch("reco_z_z_pt",                &reco_z_.z_pt);
    tree_->Branch("reco_z_z_y",                 &reco_z_.z_y);
    tree_->Branch("reco_z_z_phi",               &reco_z_.z_phi);
    tree_->Branch("reco_z_z_phistar",           &reco_z_.z_phistar);           
    tree_->Branch("reco_z_z_eta",               &reco_z_.z_eta);               
    tree_->Branch("reco_z_z_vtx_prob",          &reco_z_.z_vtx_prob);            
    tree_->Branch("reco_z_z_vtx_x",             &reco_z_.z_vtx_x);             
    tree_->Branch("reco_z_z_vtx_y",             &reco_z_.z_vtx_y);             
    tree_->Branch("reco_z_z_vtx_z",             &reco_z_.z_vtx_z);          
    tree_->Branch("reco_z_daughter0_pt",        &reco_z_.daughter0_pt);        
    tree_->Branch("reco_z_daughter0_eta",       &reco_z_.daughter0_eta);      
    tree_->Branch("reco_z_daughter0_phi",       &reco_z_.daughter0_phi);    
    tree_->Branch("reco_z_daughter1_pt",        &reco_z_.daughter1_pt);          
    tree_->Branch("reco_z_daughter1_eta",       &reco_z_.daughter1_eta);        
    tree_->Branch("reco_z_daughter1_phi",       &reco_z_.daughter1_phi);        
    tree_->Branch("reco_z_daughter0_charge",    &reco_z_.daughter0_charge);      
    tree_->Branch("reco_z_daughter1_charge",    &reco_z_.daughter1_charge);     

    tree_->Branch("reco_z_from_muons_z_m",                 &reco_z_from_muons_.z_m);
    tree_->Branch("reco_z_from_muons_z_pt",                &reco_z_from_muons_.z_pt);
    tree_->Branch("reco_z_from_muons_z_y",                 &reco_z_from_muons_.z_y);
    tree_->Branch("reco_z_from_muons_z_phi",               &reco_z_from_muons_.z_phi);
    tree_->Branch("reco_z_from_muons_z_phistar",           &reco_z_from_muons_.z_phistar);           
    tree_->Branch("reco_z_from_muons_z_eta",               &reco_z_from_muons_.z_eta);               
    tree_->Branch("reco_z_from_muons_z_vtx_prob",          &reco_z_from_muons_.z_vtx_prob);            
    tree_->Branch("reco_z_from_muons_z_vtx_x",             &reco_z_from_muons_.z_vtx_x);             
    tree_->Branch("reco_z_from_muons_z_vtx_y",             &reco_z_from_muons_.z_vtx_y);             
    tree_->Branch("reco_z_from_muons_z_vtx_z",             &reco_z_from_muons_.z_vtx_z);          
    tree_->Branch("reco_z_from_muons_daughter0_pt",        &reco_z_from_muons_.daughter0_pt);        
    tree_->Branch("reco_z_from_muons_daughter0_eta",       &reco_z_from_muons_.daughter0_eta);      
    tree_->Branch("reco_z_from_muons_daughter0_phi",       &reco_z_from_muons_.daughter0_phi);    
    tree_->Branch("reco_z_from_muons_daughter1_pt",        &reco_z_from_muons_.daughter1_pt);          
    tree_->Branch("reco_z_from_muons_daughter1_eta",       &reco_z_from_muons_.daughter1_eta);        
    tree_->Branch("reco_z_from_muons_daughter1_phi",       &reco_z_from_muons_.daughter1_phi);        
    tree_->Branch("reco_z_from_muons_daughter0_charge",    &reco_z_from_muons_.daughter0_charge);      
    tree_->Branch("reco_z_from_muons_daughter1_charge",    &reco_z_from_muons_.daughter1_charge);     

    tree_->Branch("reco_jpsi_jpsi_m",                  &reco_jpsi_.jpsi_m);                
    tree_->Branch("reco_jpsi_jpsi_pt",                 &reco_jpsi_.jpsi_pt);             
    tree_->Branch("reco_jpsi_jpsi_y",                  &reco_jpsi_.jpsi_y);            
    tree_->Branch("reco_jpsi_jpsi_phi",                &reco_jpsi_.jpsi_phi);           
    tree_->Branch("reco_jpsi_jpsi_eta",                &reco_jpsi_.jpsi_eta);
    tree_->Branch("reco_jpsi_jpsi_vtx_prob",           &reco_jpsi_.jpsi_vtx_prob);
    tree_->Branch("reco_jpsi_jpsi_vtx_x",              &reco_jpsi_.jpsi_vtx_x);
    tree_->Branch("reco_jpsi_jpsi_vtx_y",              &reco_jpsi_.jpsi_vtx_y);          
    tree_->Branch("reco_jpsi_jpsi_vtx_z",              &reco_jpsi_.jpsi_vtx_z);
    tree_->Branch("reco_jpsi_jpsi_tau_xy",             &reco_jpsi_.jpsi_tau_xy);
    tree_->Branch("reco_jpsi_jpsi_tau_z",              &reco_jpsi_.jpsi_tau_z);
    tree_->Branch("reco_jpsi_jpsi_distance_xy",        &reco_jpsi_.jpsi_distance_xy);
    tree_->Branch("reco_jpsi_jpsi_distance_z",         &reco_jpsi_.jpsi_distance_z);
    tree_->Branch("reco_jpsi_jpsi_eff",                &reco_jpsi_.jpsi_eff);
    tree_->Branch("reco_jpsi_jpsi_acc_eff",            &reco_jpsi_.jpsi_acc_eff);
    tree_->Branch("reco_jpsi_jpsi_scale_factor",       &reco_jpsi_.jpsi_scale_factor);
    tree_->Branch("reco_jpsi_muon0_pt",                &reco_jpsi_.muon0_pt);
    tree_->Branch("reco_jpsi_muon0_eta",               &reco_jpsi_.muon0_eta);
    tree_->Branch("reco_jpsi_muon0_phi",               &reco_jpsi_.muon0_phi);
    tree_->Branch("reco_jpsi_muon1_pt",                &reco_jpsi_.muon1_pt);
    tree_->Branch("reco_jpsi_muon1_eta",               &reco_jpsi_.muon1_eta);
    tree_->Branch("reco_jpsi_muon1_phi",               &reco_jpsi_.muon1_phi);
    tree_->Branch("reco_jpsi_muon0_charge",            &reco_jpsi_.muon0_charge);
    tree_->Branch("reco_jpsi_muon1_charge",            &reco_jpsi_.muon1_charge);
    tree_->Branch("reco_jpsi_has_muons_in_eta_window", &reco_jpsi_.has_muons_in_eta_window);
    tree_->Branch("reco_jpsi_has_high_pt_muons",       &reco_jpsi_.has_high_pt_muons);

    tree_->Branch("reco_jpsi_from_electrons_jpsi_m",                  &reco_jpsi_from_electrons_.jpsi_m);   
    tree_->Branch("reco_jpsi_from_electrons_jpsi_pt",                 &reco_jpsi_from_electrons_.jpsi_pt);  
    tree_->Branch("reco_jpsi_from_electrons_jpsi_y",                  &reco_jpsi_from_electrons_.jpsi_y);   
    tree_->Branch("reco_jpsi_from_electrons_jpsi_phi",                &reco_jpsi_from_electrons_.jpsi_phi); 
    tree_->Branch("reco_jpsi_from_electrons_jpsi_eta",                &reco_jpsi_from_electrons_.jpsi_eta); 
    tree_->Branch("reco_jpsi_from_electrons_jpsi_vtx_prob",           &reco_jpsi_from_electrons_.jpsi_vtx_prob);      
    tree_->Branch("reco_jpsi_from_electrons_jpsi_vtx_x",              &reco_jpsi_from_electrons_.jpsi_vtx_x); 
    tree_->Branch("reco_jpsi_from_electrons_jpsi_vtx_y",              &reco_jpsi_from_electrons_.jpsi_vtx_y); 
    tree_->Branch("reco_jpsi_from_electrons_jpsi_vtx_z",              &reco_jpsi_from_electrons_.jpsi_vtx_z); 
    tree_->Branch("reco_jpsi_from_electrons_jpsi_tau_xy",             &reco_jpsi_from_electrons_.jpsi_tau_xy);
    tree_->Branch("reco_jpsi_from_electrons_jpsi_tau_z",              &reco_jpsi_from_electrons_.jpsi_tau_z); 
    tree_->Branch("reco_jpsi_from_electrons_jpsi_distance_xy",        &reco_jpsi_from_electrons_.jpsi_distance_xy);   
    tree_->Branch("reco_jpsi_from_electrons_jpsi_distance_z",         &reco_jpsi_from_electrons_.jpsi_distance_z);    
    tree_->Branch("reco_jpsi_from_electrons_jpsi_eff",                &reco_jpsi_from_electrons_.jpsi_eff);
    tree_->Branch("reco_jpsi_from_electrons_jpsi_acc_eff",            &reco_jpsi_from_electrons_.jpsi_acc_eff);
    tree_->Branch("reco_jpsi_from_electrons_jpsi_scale_factor",       &reco_jpsi_from_electrons_.jpsi_scale_factor);  
    tree_->Branch("reco_jpsi_from_electrons_muon0_pt",                &reco_jpsi_from_electrons_.muon0_pt); 
    tree_->Branch("reco_jpsi_from_electrons_muon0_eta",               &reco_jpsi_from_electrons_.muon0_eta);
    tree_->Branch("reco_jpsi_from_electrons_muon0_phi",               &reco_jpsi_from_electrons_.muon0_phi);
    tree_->Branch("reco_jpsi_from_electrons_muon1_pt",                &reco_jpsi_from_electrons_.muon1_pt); 
    tree_->Branch("reco_jpsi_from_electrons_muon1_eta",               &reco_jpsi_from_electrons_.muon1_eta);
    tree_->Branch("reco_jpsi_from_electrons_muon1_phi",               &reco_jpsi_from_electrons_.muon1_phi);
    tree_->Branch("reco_jpsi_from_electrons_muon0_charge",            &reco_jpsi_from_electrons_.muon0_charge); 
    tree_->Branch("reco_jpsi_from_electrons_muon1_charge",            &reco_jpsi_from_electrons_.muon1_charge);
    tree_->Branch("reco_jpsi_from_electrons_has_muons_in_eta_window", &reco_jpsi_from_electrons_.has_muons_in_eta_window); 
    tree_->Branch("reco_jpsi_from_electrons_has_high_pt_muons",       &reco_jpsi_from_electrons_.has_high_pt_muons); 

    if (IS_MC_) {
      tree_->Branch("truth_z_z_m",                 &truth_z_electrons_.z_m);
      tree_->Branch("truth_z_z_pt",                &truth_z_electrons_.z_pt);
      tree_->Branch("truth_z_z_y",                 &truth_z_electrons_.z_y);
      tree_->Branch("truth_z_z_phi",               &truth_z_electrons_.z_phi);
      tree_->Branch("truth_z_z_phistar",           &truth_z_electrons_.z_phistar);           
      tree_->Branch("truth_z_z_eta",               &truth_z_electrons_.z_eta);               
      tree_->Branch("truth_z_z_vtx_prob",          &truth_z_electrons_.z_vtx_prob);            
      tree_->Branch("truth_z_z_vtx_x",             &truth_z_electrons_.z_vtx_x);             
      tree_->Branch("truth_z_z_vtx_y",             &truth_z_electrons_.z_vtx_y);             
      tree_->Branch("truth_z_z_vtx_z",             &truth_z_electrons_.z_vtx_z);          
      tree_->Branch("truth_z_daughter0_pt",        &truth_z_electrons_.daughter0_pt);        
      tree_->Branch("truth_z_daughter0_eta",       &truth_z_electrons_.daughter0_eta);      
      tree_->Branch("truth_z_daughter0_phi",       &truth_z_electrons_.daughter0_phi);    
      tree_->Branch("truth_z_daughter1_pt",        &truth_z_electrons_.daughter1_pt);          
      tree_->Branch("truth_z_daughter1_eta",       &truth_z_electrons_.daughter1_eta);        
      tree_->Branch("truth_z_daughter1_phi",       &truth_z_electrons_.daughter1_phi);        
      tree_->Branch("truth_z_daughter0_charge",    &truth_z_electrons_.daughter0_charge);      
      tree_->Branch("truth_z_daughter1_charge",    &truth_z_electrons_.daughter1_charge);     

      tree_->Branch("truth_z_from_muons_z_m",                 &truth_z_muons_.z_m);
      tree_->Branch("truth_z_from_muons_z_pt",                &truth_z_muons_.z_pt);
      tree_->Branch("truth_z_from_muons_z_y",                 &truth_z_muons_.z_y);
      tree_->Branch("truth_z_from_muons_z_phi",               &truth_z_muons_.z_phi);
      tree_->Branch("truth_z_from_muons_z_phistar",           &truth_z_muons_.z_phistar);           
      tree_->Branch("truth_z_from_muons_z_eta",               &truth_z_muons_.z_eta);               
      tree_->Branch("truth_z_from_muons_z_vtx_prob",          &truth_z_muons_.z_vtx_prob);            
      tree_->Branch("truth_z_from_muons_z_vtx_x",             &truth_z_muons_.z_vtx_x);             
      tree_->Branch("truth_z_from_muons_z_vtx_y",             &truth_z_muons_.z_vtx_y);             
      tree_->Branch("truth_z_from_muons_z_vtx_z",             &truth_z_muons_.z_vtx_z);          
      tree_->Branch("truth_z_from_muons_daughter0_pt",        &truth_z_muons_.daughter0_pt);        
      tree_->Branch("truth_z_from_muons_daughter0_eta",       &truth_z_muons_.daughter0_eta);      
      tree_->Branch("truth_z_from_muons_daughter0_phi",       &truth_z_muons_.daughter0_phi);    
      tree_->Branch("truth_z_from_muons_daughter1_pt",        &truth_z_muons_.daughter1_pt);          
      tree_->Branch("truth_z_from_muons_daughter1_eta",       &truth_z_muons_.daughter1_eta);        
      tree_->Branch("truth_z_from_muons_daughter1_phi",       &truth_z_muons_.daughter1_phi);        
      tree_->Branch("truth_z_from_muons_daughter0_charge",    &truth_z_muons_.daughter0_charge);      
      tree_->Branch("truth_z_from_muons_daughter1_charge",    &truth_z_muons_.daughter1_charge);     

      tree_->Branch("truth_jpsi_jpsi_m",                  &truth_jpsi_.jpsi_m);   
      tree_->Branch("truth_jpsi_jpsi_pt",                 &truth_jpsi_.jpsi_pt);  
      tree_->Branch("truth_jpsi_jpsi_y",                  &truth_jpsi_.jpsi_y);   
      tree_->Branch("truth_jpsi_jpsi_phi",                &truth_jpsi_.jpsi_phi); 
      tree_->Branch("truth_jpsi_jpsi_eta",                &truth_jpsi_.jpsi_eta); 
      tree_->Branch("truth_jpsi_jpsi_vtx_prob",           &truth_jpsi_.jpsi_vtx_prob);      
      tree_->Branch("truth_jpsi_jpsi_vtx_x",              &truth_jpsi_.jpsi_vtx_x); 
      tree_->Branch("truth_jpsi_jpsi_vtx_y",              &truth_jpsi_.jpsi_vtx_y); 
      tree_->Branch("truth_jpsi_jpsi_vtx_z",              &truth_jpsi_.jpsi_vtx_z); 
      tree_->Branch("truth_jpsi_jpsi_tau_xy",             &truth_jpsi_.jpsi_tau_xy);
      tree_->Branch("truth_jpsi_jpsi_tau_z",              &truth_jpsi_.jpsi_tau_z); 
      tree_->Branch("truth_jpsi_jpsi_distance_xy",        &truth_jpsi_.jpsi_distance_xy);   
      tree_->Branch("truth_jpsi_jpsi_distance_z",         &truth_jpsi_.jpsi_distance_z);    
      tree_->Branch("truth_jpsi_jpsi_eff",                &truth_jpsi_.jpsi_eff);
      tree_->Branch("truth_jpsi_jpsi_acc_eff",            &truth_jpsi_.jpsi_acc_eff);
      tree_->Branch("truth_jpsi_jpsi_scale_factor",       &truth_jpsi_.jpsi_scale_factor);  
      tree_->Branch("truth_jpsi_muon0_pt",                &truth_jpsi_.muon0_pt); 
      tree_->Branch("truth_jpsi_muon0_eta",               &truth_jpsi_.muon0_eta);
      tree_->Branch("truth_jpsi_muon0_phi",               &truth_jpsi_.muon0_phi);
      tree_->Branch("truth_jpsi_muon1_pt",                &truth_jpsi_.muon1_pt); 
      tree_->Branch("truth_jpsi_muon1_eta",               &truth_jpsi_.muon1_eta);
      tree_->Branch("truth_jpsi_muon1_phi",               &truth_jpsi_.muon1_phi);
      tree_->Branch("truth_jpsi_muon0_charge",            &truth_jpsi_.muon0_charge); 
      tree_->Branch("truth_jpsi_muon1_charge",            &truth_jpsi_.muon1_charge);
      tree_->Branch("truth_jpsi_has_muons_in_eta_window", &truth_jpsi_.has_muons_in_eta_window); 
      tree_->Branch("truth_jpsi_has_high_pt_muons",       &truth_jpsi_.has_high_pt_muons); 

    }

    //-----------------
    tree_->Branch("event_weight",                                                                &event_.event_weight);
    tree_->Branch("event_number",                                                                &event_.event_number);
    tree_->Branch("run_number",                                                                  &event_.run_number);  
    tree_->Branch("n_verts",                                                                     &event_.n_verts);     
    tree_->Branch("found_good_muons_from_z",                                                     &event_.found_good_muons_from_z);   
    tree_->Branch("found_dimuon_z_compatible_vertex",                                            &event_.found_dimuon_z_compatible_vertex);
    tree_->Branch("found_good_electrons_from_z",                                                 &event_.found_good_electrons_from_z);    tree_->Branch("found_dielectron_z_compatible_vertex",                                            &event_.found_dielectron_z_compatible_vertex);       
    tree_->Branch("found_dimuon_jpsi_with_soft_id_and_high_pt_muons",                            &event_.found_dimuon_jpsi_with_soft_id_and_high_pt_muons);                           
    tree_->Branch("found_dimuon_jpsi_with_good_muons_and_compatible_muon_vertex",                &event_.found_dimuon_jpsi_with_good_muons_and_compatible_muon_vertex);               
    tree_->Branch("found_good_dimuon_jpsi_compatible_with_primary_vertex",                       &event_.found_good_dimuon_jpsi_compatible_with_primary_vertex);                      
    tree_->Branch("found_dimuon_jpsi_from_electrons_with_muons_in_eta_window",                   &event_.found_dimuon_jpsi_from_electrons_with_muons_in_eta_window);
    tree_->Branch("found_dimuon_jpsi_from_electrons_with_high_pt_muons",                         &event_.found_dimuon_jpsi_from_electrons_with_high_pt_muons);      
    tree_->Branch("found_dimuon_jpsi_from_electrons_with_soft_id_and_high_pt_muons",             &event_.found_dimuon_jpsi_from_electrons_with_soft_id_and_high_pt_muons);
    tree_->Branch("found_dimuon_jpsi_from_electrons_with_good_muons_and_compatible_muon_vertex", &event_.found_dimuon_jpsi_from_electrons_with_good_muons_and_compatible_muon_vertex);
    tree_->Branch("found_good_dimuon_jpsi_from_electrons_compatible_with_primary_vertex",        &event_.found_good_dimuon_jpsi_from_electrons_compatible_with_primary_vertex);       
    //-------------------

    //if (IS_MC_) {
    //  tree_->Branch("weight_size", &weight_size_, "weight_size/I");
    //  tree_->Branch("weights", weights_, "weights[weight_size]/D");
    //  tree_->Branch("weight_ids", weight_ids_, "weight_ids[weight_size]/I");
    //  tree_->Branch("weight_cteq_size", &weight_cteq_size_, "weight_cteq_size/I");
    //  tree_->Branch("weights_cteq", weights_cteq_, "weights_cteq[weight_cteq_size]/D");
    //  tree_->Branch("weight_mstw_size", &weight_mstw_size_, "weight_mstw_size/I");
    //  tree_->Branch("weights_mstw", weights_mstw_, "weights_mstw[weight_mstw_size]/D");
    //  tree_->Branch("weight_nnpdf_size", &weight_nnpdf_size_, "weight_nnpdf_size/I");
    //  tree_->Branch("weights_nnpdf", weights_nnpdf_, "weights_nnpdf[weight_nnpdf_size]/D");
    //  tree_->Branch("weight_fsr", &weight_fsr_, "weight_fsr/D");
    //}
  }

  ZFinderTree::~ZFinderTree() {
    tree_->Write();
    // Clean up our pointer
    delete tree_;
  }

  void ZFinderTree::Fill(const ZFinderEvent& zfe) {
    
    // Clear our branches
    reco_z_.clear_values();
    reco_z_from_muons_.clear_values();
    truth_z_electrons_.clear_values();
    truth_z_muons_.clear_values();
    reco_jpsi_.clear_values();
    reco_jpsi_from_electrons_.clear_values();
    truth_jpsi_.clear_values();
    event_.clear_values();

    // Reco
    
    if (true) {
      int n_jpsi = 0;
// sleontsi from here down use for Z->ee
      for (unsigned int iZ = 0; iZ < zfe.reco_z.muon0_pT.size() ; ++iZ ) {
        math::PtEtaPhiMLorentzVector Z_l1(zfe.reco_z.muon0_pT.at(iZ), 
                                          zfe.reco_z.muon0_eta.at(iZ),
                                          zfe.reco_z.muon0_phi.at(iZ),
                                          0.510998928/1000.);
        math::PtEtaPhiMLorentzVector Z_l2(zfe.reco_z.muon1_pT.at(iZ), 
                                          zfe.reco_z.muon1_eta.at(iZ),
                                          zfe.reco_z.muon1_phi.at(iZ),
                                          0.510998928/1000.);

// sleontsi from here down use for Z->mumu
 //     for (unsigned int iZ = 0; iZ < zfe.reco_z_from_muons.muon0_pT.size() ; ++iZ ) {
 //       math::PtEtaPhiMLorentzVector Z_l1(zfe.reco_z_from_muons.muon0_pT.at(iZ), 
 //                                         zfe.reco_z_from_muons.muon0_eta.at(iZ),
 //                                         zfe.reco_z_from_muons.muon0_phi.at(iZ),
 //                                         0.1056583715);
 //       math::PtEtaPhiMLorentzVector Z_l2(zfe.reco_z_from_muons.muon1_pT.at(iZ), 
 //                                         zfe.reco_z_from_muons.muon1_eta.at(iZ),
 //                                         zfe.reco_z_from_muons.muon1_phi.at(iZ),
 //                                         0.1056583715);
        double fourMass = 0.;
        for (unsigned int i = 0; i < zfe.reco_jpsi.m.size() ; ++i ) {
          n_jpsi++;

          math::PtEtaPhiMLorentzVector Jpsi_l1(zfe.reco_jpsi.muon0.at(i).pt(), 
                                               zfe.reco_jpsi.muon0.at(i).eta(),
                                               zfe.reco_jpsi.muon0.at(i).phi(),
                                               0.1056583715);
          math::PtEtaPhiMLorentzVector Jpsi_l2(zfe.reco_jpsi.muon1.at(i).pt(), 
                                               zfe.reco_jpsi.muon1.at(i).eta(),
                                               zfe.reco_jpsi.muon1.at(i).phi(),
                                               0.1056583715);
          fourMass = (Z_l1+Z_l2+Jpsi_l1+Jpsi_l2).mass();
          if (fourMass<66. || fourMass>116.)
            continue;

          if (zfe.reco_jpsi.m.at(i)<2.6 || zfe.reco_jpsi.m.at(i)>3.6)
            continue;

          if (Z_l1.pt()<MIN_Z_MUON_PT || Z_l2.pt()<MIN_Z_MUON_PT || Jpsi_l1.pt()<MIN_JPSI_LEADING_MUON_PT || Jpsi_l2.pt()<MIN_JPSI_LEADING_MUON_PT)
            continue;

//          if (Z_l1.pt()==Jpsi_l1.pt() || Z_l1.pt()==Jpsi_l2.pt() || Z_l2.pt()==Jpsi_l1.pt() || Z_l2.pt()==Jpsi_l2.pt())
//            continue;

          // Z->ee
          reco_z_.z_m        .push_back((Z_l1+Z_l2).mass());
          reco_z_.z_pt       .push_back((Z_l1+Z_l2).pt());
          reco_z_.z_y        .push_back(zfe.reco_z.y);
          reco_z_.z_phi      .push_back(zfe.reco_z.phi);
          reco_z_.z_phistar  .push_back(zfe.reco_z.phistar);
          reco_z_.z_eta      .push_back((Z_l1+Z_l2).eta());
          reco_z_.z_vtx_prob .push_back(zfe.reco_z.vtx_prob.at(iZ));
          reco_z_.z_vtx_x    .push_back(zfe.reco_z.vtx_x);
          reco_z_.z_vtx_y    .push_back(zfe.reco_z.vtx_y);
          reco_z_.z_vtx_z    .push_back(zfe.reco_z.vtx_z);
          reco_z_.daughter0_pt     .push_back(zfe.reco_z.muon0_pT.at(iZ));
          reco_z_.daughter0_eta    .push_back(zfe.reco_z.muon0_eta.at(iZ));
          reco_z_.daughter0_phi    .push_back(zfe.reco_z.muon0_phi.at(iZ));
          reco_z_.daughter1_pt     .push_back(zfe.reco_z.muon1_pT.at(iZ));
          reco_z_.daughter1_eta    .push_back(zfe.reco_z.muon1_eta.at(iZ));
          reco_z_.daughter1_phi    .push_back(zfe.reco_z.muon1_phi.at(iZ));
          reco_z_.daughter0_charge .push_back(zfe.e0->charge);
          reco_z_.daughter1_charge .push_back(zfe.e1->charge);

          // z->mumu
 //         reco_z_from_muons_.z_m        .push_back((Z_l1+Z_l2).mass());
 //         reco_z_from_muons_.z_pt       .push_back((Z_l1+Z_l2).pt());
 //         reco_z_from_muons_.z_y        .push_back(zfe.reco_z_from_muons.y);
 //         reco_z_from_muons_.z_phi      .push_back(zfe.reco_z_from_muons.phi);
 //         reco_z_from_muons_.z_phistar  .push_back(zfe.reco_z_from_muons.phistar);
 //         reco_z_from_muons_.z_eta      .push_back((Z_l1+Z_l2).eta());
 //         reco_z_from_muons_.z_vtx_prob .push_back(zfe.reco_z_from_muons.vtx_prob.at(iZ));
 //         reco_z_from_muons_.z_vtx_x    .push_back(zfe.reco_z_from_muons.vtx_x);
 //         reco_z_from_muons_.z_vtx_y    .push_back(zfe.reco_z_from_muons.vtx_y);
 //         reco_z_from_muons_.z_vtx_z    .push_back(zfe.reco_z_from_muons.vtx_z);
 //         reco_z_from_muons_.daughter0_pt     .push_back(zfe.reco_z_from_muons.muon0_pT.at(iZ));
 //         reco_z_from_muons_.daughter0_eta    .push_back(zfe.reco_z_from_muons.muon0_eta.at(iZ));
 //         reco_z_from_muons_.daughter0_phi    .push_back(zfe.reco_z_from_muons.muon0_phi.at(iZ));
 //         reco_z_from_muons_.daughter1_pt     .push_back(zfe.reco_z_from_muons.muon1_pT.at(iZ));
 //         reco_z_from_muons_.daughter1_eta    .push_back(zfe.reco_z_from_muons.muon1_eta.at(iZ));
 //         reco_z_from_muons_.daughter1_phi    .push_back(zfe.reco_z_from_muons.muon1_phi.at(iZ));
 //         reco_z_from_muons_.daughter0_charge .push_back(zfe.z_muon0.charge());
 //         reco_z_from_muons_.daughter1_charge .push_back(zfe.z_muon1.charge());

          // jpsi->mumu
          reco_jpsi_.jpsi_m                  .push_back(zfe.reco_jpsi.m.at(i));
          reco_jpsi_.jpsi_pt                 .push_back(zfe.reco_jpsi.pt.at(i));
          reco_jpsi_.jpsi_y                  .push_back(zfe.reco_jpsi.y.at(i));
          reco_jpsi_.jpsi_phi                .push_back(zfe.reco_jpsi.phi.at(i));
          reco_jpsi_.jpsi_eta                .push_back(zfe.reco_jpsi.eta.at(i));
          reco_jpsi_.jpsi_tau_xy             .push_back(zfe.reco_jpsi.tau_xy.at(i));
          reco_jpsi_.jpsi_tau_z              .push_back(zfe.reco_jpsi.tau_z.at(i));
          reco_jpsi_.jpsi_distance_xy        .push_back(zfe.reco_jpsi.distance_xy.at(i));
          reco_jpsi_.jpsi_distance_z         .push_back(zfe.reco_jpsi.distance_z.at(i));
          reco_jpsi_.jpsi_vtx_prob           .push_back(zfe.reco_jpsi.y.at(i));
          reco_jpsi_.jpsi_vtx_x              .push_back(zfe.reco_jpsi.vtx_x.at(i));
          reco_jpsi_.jpsi_vtx_y              .push_back(zfe.reco_jpsi.vtx_y.at(i));
          reco_jpsi_.jpsi_vtx_z              .push_back(zfe.reco_jpsi.vtx_z.at(i));
          reco_jpsi_.jpsi_eff                .push_back(zfe.reco_jpsi.jpsi_efficiency.at(i));
          reco_jpsi_.jpsi_acc_eff            .push_back(zfe.reco_jpsi.jpsi_acc_eff.at(i));
          reco_jpsi_.jpsi_scale_factor       .push_back(zfe.reco_jpsi.jpsi_scale_factor.at(i));
          reco_jpsi_.muon0_pt                .push_back(zfe.reco_jpsi.muon0.at(i).pt());
          reco_jpsi_.muon0_eta               .push_back(zfe.reco_jpsi.muon0.at(i).eta());
          reco_jpsi_.muon0_phi               .push_back(zfe.reco_jpsi.muon0.at(i).phi());
          reco_jpsi_.muon1_pt                .push_back(zfe.reco_jpsi.muon1.at(i).pt());
          reco_jpsi_.muon1_eta               .push_back(zfe.reco_jpsi.muon1.at(i).eta());
          reco_jpsi_.muon1_phi               .push_back(zfe.reco_jpsi.muon1.at(i).phi());
          reco_jpsi_.muon0_charge            .push_back(zfe.reco_jpsi.muon0.at(i).charge());
          reco_jpsi_.muon1_charge            .push_back(zfe.reco_jpsi.muon1.at(i).charge());
          reco_jpsi_.has_muons_in_eta_window .push_back(zfe.reco_jpsi.has_muons_in_eta_window.at(i));
          reco_jpsi_.has_high_pt_muons       .push_back(zfe.reco_jpsi.has_high_pt_muons.at(i));

          // event info
          event_.event_weight                                                .push_back(zfe.event_weight);
          event_.event_number                                                .push_back(zfe.id.event_num);
          event_.run_number                                                  .push_back(zfe.id.run_num);
          event_.n_verts                                                     .push_back(zfe.reco_vert.num);
          event_.found_good_muons_from_z                                     .push_back(zfe.found_good_muons_from_z);
          event_.found_dimuon_z_compatible_vertex                            .push_back(zfe.found_dimuon_z_compatible_vertex);
          event_.found_good_electrons_from_z                                 .push_back(zfe.found_good_electrons_from_z);
          event_.found_dielectron_z_compatible_vertex                        .push_back(zfe.found_dielectron_z_compatible_vertex);
          event_.found_dimuon_jpsi_with_soft_id_and_high_pt_muons            .push_back(zfe.found_dimuon_jpsi_with_soft_id_and_high_pt_muons.at(i));
          event_.found_dimuon_jpsi_with_good_muons_and_compatible_muon_vertex.push_back(zfe.found_dimuon_jpsi_with_good_muons_and_compatible_muon_vertex.at(i));
          event_.found_good_dimuon_jpsi_compatible_with_primary_vertex       .push_back(zfe.found_good_dimuon_jpsi_compatible_with_primary_vertex.at(i));

        }
      }
    }

    //---------------------------------------------------------------------------------------
    if (false) {
//      int n_jpsi_from_electrons = 0;
      for (unsigned int iZ = 0; iZ < zfe.reco_z_from_muons.muon0_pT.size() ; ++iZ ) {
        math::PtEtaPhiMLorentzVector Z_l1(zfe.reco_z_from_muons.muon0_pT.at(iZ), 
                                          zfe.reco_z_from_muons.muon0_eta.at(iZ),
                                          zfe.reco_z_from_muons.muon0_phi.at(iZ),
                                          0.1056583715);
        math::PtEtaPhiMLorentzVector Z_l2(zfe.reco_z_from_muons.muon1_pT.at(iZ), 
                                          zfe.reco_z_from_muons.muon1_eta.at(iZ),
                                          zfe.reco_z_from_muons.muon1_phi.at(iZ),
                                          0.1056583715);
        double fourMass = 0.;
        for (unsigned int i = 0; i < zfe.reco_jpsi_from_electrons.m.size() ; ++i ) {
          math::PtEtaPhiMLorentzVector Jpsi_l1(zfe.reco_jpsi_from_electrons.muon0_pT.at(i), 
                                               zfe.reco_jpsi_from_electrons.muon0_eta.at(i),
                                               zfe.reco_jpsi_from_electrons.muon0_phi.at(i),
                                               0.510998928/1000.);
          math::PtEtaPhiMLorentzVector Jpsi_l2(zfe.reco_jpsi_from_electrons.muon1_pT.at(i), 
                                               zfe.reco_jpsi_from_electrons.muon1_eta.at(i),
                                               zfe.reco_jpsi_from_electrons.muon1_phi.at(i),
                                               0.510998928/1000.);
          fourMass = (Z_l1+Z_l2+Jpsi_l1+Jpsi_l2).mass();
          if (fourMass<66. || fourMass>116.)
            continue;

          if (zfe.reco_jpsi_from_electrons.m.at(i)<2.6 || zfe.reco_jpsi_from_electrons.m.at(i)>3.6)
            continue;
 
          reco_jpsi_from_electrons_.jpsi_m                 .push_back(zfe.reco_jpsi_from_electrons.m.at(i));
          reco_jpsi_from_electrons_.jpsi_pt                .push_back(zfe.reco_jpsi_from_electrons.pt.at(i));
          reco_jpsi_from_electrons_.jpsi_y                 .push_back(zfe.reco_jpsi_from_electrons.y.at(i));
          reco_jpsi_from_electrons_.jpsi_phi               .push_back(zfe.reco_jpsi_from_electrons.phi.at(i));
          reco_jpsi_from_electrons_.jpsi_eta               .push_back(zfe.reco_jpsi_from_electrons.eta.at(i));
          reco_jpsi_from_electrons_.jpsi_tau_xy            .push_back(zfe.reco_jpsi_from_electrons.tau_xy.at(i));
          reco_jpsi_from_electrons_.jpsi_tau_z             .push_back(zfe.reco_jpsi_from_electrons.tau_z.at(i));
          reco_jpsi_from_electrons_.jpsi_distance_xy       .push_back(zfe.reco_jpsi_from_electrons.distance_xy.at(i));
          reco_jpsi_from_electrons_.jpsi_distance_z        .push_back(zfe.reco_jpsi_from_electrons.distance_z.at(i));
          reco_jpsi_from_electrons_.jpsi_vtx_prob          .push_back(zfe.reco_jpsi_from_electrons.y.at(i));
          reco_jpsi_from_electrons_.jpsi_vtx_x             .push_back(zfe.reco_jpsi_from_electrons.vtx_x.at(i));
          reco_jpsi_from_electrons_.jpsi_vtx_y             .push_back(zfe.reco_jpsi_from_electrons.vtx_y.at(i));
          reco_jpsi_from_electrons_.jpsi_vtx_z             .push_back(zfe.reco_jpsi_from_electrons.vtx_z.at(i));
          reco_jpsi_from_electrons_.muon0_pt               .push_back(zfe.reco_jpsi_from_electrons.muon0_pT.at(i));
          reco_jpsi_from_electrons_.muon0_eta              .push_back(zfe.reco_jpsi_from_electrons.muon0_eta.at(i));
          reco_jpsi_from_electrons_.muon0_phi              .push_back(zfe.reco_jpsi_from_electrons.muon0_phi.at(i));
          reco_jpsi_from_electrons_.muon1_pt               .push_back(zfe.reco_jpsi_from_electrons.muon1_pT.at(i));
          reco_jpsi_from_electrons_.muon1_eta              .push_back(zfe.reco_jpsi_from_electrons.muon1_eta.at(i));
          reco_jpsi_from_electrons_.muon1_phi              .push_back(zfe.reco_jpsi_from_electrons.muon1_phi.at(i));
          reco_jpsi_from_electrons_.has_muons_in_eta_window.push_back(zfe.reco_jpsi_from_electrons.has_muons_in_eta_window.at(i));
          reco_jpsi_from_electrons_.has_high_pt_muons      .push_back(zfe.reco_jpsi_from_electrons.has_high_pt_muons.at(i));

        }
      }
    }

    // Truth
//    if (!zfe.is_real_data) {
//      int n_jpsi = 0;
//      for (unsigned int i = 0; i < zfe.truth_jpsi.m.size() ; ++i ) {
//        n_jpsi++;
//        //TODO decide on how to do this for a tree (and histogram method too), for now just take first jpsi
//        //this is definitely a kludge, think of proper way to do this truth should only ever have 1 jpsi
//        if (n_jpsi > 1) {
//          continue;
//        }
//        truth_jpsi_.jpsi_m = zfe.truth_jpsi.m.at(i);
//        truth_jpsi_.jpsi_pt = zfe.truth_jpsi.pt.at(i);
//        truth_jpsi_.jpsi_y = zfe.truth_jpsi.y.at(i);
//        truth_jpsi_.jpsi_phi = zfe.truth_jpsi.phi.at(i);
//        truth_jpsi_.jpsi_eta = zfe.truth_jpsi.eta.at(i);
//        truth_jpsi_.jpsi_vtx_x = zfe.truth_jpsi.vtx_x.at(i);
//        truth_jpsi_.jpsi_vtx_y = zfe.truth_jpsi.vtx_y.at(i);
//        truth_jpsi_.jpsi_vtx_z = zfe.truth_jpsi.vtx_z.at(i);
//        truth_jpsi_.muon0_pt = zfe.jpsi_muon0.at(i)->pt();
//        truth_jpsi_.muon0_eta = zfe.jpsi_muon0.at(i)->eta();
//        truth_jpsi_.muon0_phi = zfe.jpsi_muon0.at(i)->phi();
//        truth_jpsi_.muon1_pt = zfe.jpsi_muon1.at(i)->pt();
//        truth_jpsi_.muon1_eta = zfe.jpsi_muon1.at(i)->eta();
//        truth_jpsi_.muon1_phi = zfe.jpsi_muon1.at(i)->phi();
//        truth_jpsi_.muon0_charge = zfe.jpsi_muon0.at(i)->charge();
//        truth_jpsi_.muon1_charge = zfe.jpsi_muon1.at(i)->charge();
//        truth_jpsi_.has_muons_in_eta_window = zfe.truth_jpsi.has_muons_in_eta_window.at(i);
//        truth_jpsi_.has_high_pt_muons = zfe.truth_jpsi.has_high_pt_muons.at(i);
//      }
//    }
    // Truth
    //if (IS_MC_ && !zf_event.is_real_data) {
    //  truth_.z_m = zf_event.truth_z.m;
    //  truth_.z_y = zf_event.truth_z.y;
    //  truth_.z_phistar_dressed = zf_event.truth_z.phistar;
    //  truth_.z_phistar_born = zf_event.truth_z.bornPhistar;
    //  truth_.z_phistar_naked = zf_event.truth_z.nakedPhistar;
    //  truth_.z_phistar_sc = zf_event.truth_z.scPhistar;
    //  truth_.z_pt = zf_event.truth_z.pt;
    //  truth_.z_eta = zf_event.truth_z.eta;
    //  truth_.n_verts = zf_event.truth_vert.num;
    //  truth_.n_true_pileup = zf_event.truth_vert.true_num;
    //  if (zf_event.e0_truth != nullptr) {
    //    truth_.e_pt[0] = zf_event.e0_truth->pt();
    //    truth_.e_eta[0] = zf_event.e0_truth->eta();
    //    truth_.e_phi[0] = zf_event.e0_truth->phi();
    //    truth_.e_rnine[0] = zf_event.e0_truth->r9();
    //    truth_.e_charge[0] = zf_event.e0_truth->charge();
    //  }
    //  if (zf_event.e1_truth != nullptr) {
    //    truth_.e_pt[1] = zf_event.e1_truth->pt();
    //    truth_.e_eta[1] = zf_event.e1_truth->eta();
    //    truth_.e_phi[1] = zf_event.e1_truth->phi();
    //    truth_.e_rnine[1] = zf_event.e1_truth->r9();
    //    truth_.e_charge[1] = zf_event.e1_truth->charge();
    //    if (zf_event.GetZDef(zdef_name_) != nullptr) {
    //      const cutlevel_vector* clv = zf_event.GetZDef(zdef_name_);
    //      truth_.t0tight = clv->back().second.t0p1_pass;
    //      truth_.t1tight = clv->back().second.t1p0_pass;
    //    }
    //  }
    //}

    // General Event info

    //------------
    event_.found_dimuon_jpsi_from_electrons_with_soft_id_and_high_pt_muons.push_back(zfe.found_dimuon_jpsi_from_electrons_with_soft_id_and_high_pt_muons);
    event_.found_dimuon_jpsi_from_electrons_with_good_muons_and_compatible_muon_vertex.push_back(zfe.found_dimuon_jpsi_from_electrons_with_good_muons_and_compatible_muon_vertex);
    event_.found_good_dimuon_jpsi_from_electrons_compatible_with_primary_vertex       .push_back(zfe.found_good_dimuon_jpsi_from_electrons_compatible_with_primary_vertex);
    //------------

    if (true) {
      tree_->Fill();
      //tree_->Write();
    }
  }

  TFile* ZFinderTree::GetCurrentFile() {
    return tree_->GetCurrentFile();
  }
}  // namespace zf
