#include "ZFinder/Event/interface/ZFinderTree.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

// Standard Library
#include <algorithm>
#include <vector>  // std::min


#include <typeinfo>


// ZFinder Code


namespace zf {
  using namespace edm;
  using namespace zf;

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
    tree_->Branch("reco_z_daughter0_d0",        &reco_z_.daughter0_d0);    
    tree_->Branch("reco_z_daughter0_dxy",       &reco_z_.daughter0_dxy);    
    tree_->Branch("reco_z_daughter0_dz",        &reco_z_.daughter0_dz);    
    tree_->Branch("reco_z_daughter1_d0",        &reco_z_.daughter1_d0);    
    tree_->Branch("reco_z_daughter1_dxy",       &reco_z_.daughter1_dxy);    
    tree_->Branch("reco_z_daughter1_dz",        &reco_z_.daughter1_dz);    
    tree_->Branch("reco_z_daughter0_d0err",     &reco_z_.daughter0_d0err);  
    tree_->Branch("reco_z_daughter0_dxyerr",    &reco_z_.daughter0_dxyerr);  
    tree_->Branch("reco_z_daughter0_dzerr",     &reco_z_.daughter0_dzerr);  
    tree_->Branch("reco_z_daughter1_d0err",     &reco_z_.daughter1_d0err);  
    tree_->Branch("reco_z_daughter1_dxyerr",    &reco_z_.daughter1_dxyerr);  
    tree_->Branch("reco_z_daughter1_dzerr",     &reco_z_.daughter1_dzerr);  
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
    tree_->Branch("reco_z_from_muons_daughter0_d0",        &reco_z_from_muons_.daughter0_d0);    
    tree_->Branch("reco_z_from_muons_daughter0_dxy",       &reco_z_from_muons_.daughter0_dxy);    
    tree_->Branch("reco_z_from_muons_daughter0_dz",        &reco_z_from_muons_.daughter0_dz);    
    tree_->Branch("reco_z_from_muons_daughter0_d0err",     &reco_z_from_muons_.daughter0_d0err);    
    tree_->Branch("reco_z_from_muons_daughter0_dxyerr",    &reco_z_from_muons_.daughter0_dxyerr);    
    tree_->Branch("reco_z_from_muons_daughter0_dzerr",     &reco_z_from_muons_.daughter0_dzerr);    
    tree_->Branch("reco_z_from_muons_daughter0_trkKink",   &reco_z_from_muons_.daughter0_trkKink);    
    tree_->Branch("reco_z_from_muons_daughter0_glbKink",   &reco_z_from_muons_.daughter0_glbKink);    
    tree_->Branch("reco_z_from_muons_daughter1_pt",        &reco_z_from_muons_.daughter1_pt);          
    tree_->Branch("reco_z_from_muons_daughter1_eta",       &reco_z_from_muons_.daughter1_eta);        
    tree_->Branch("reco_z_from_muons_daughter1_phi",       &reco_z_from_muons_.daughter1_phi);        
    tree_->Branch("reco_z_from_muons_daughter1_d0",        &reco_z_from_muons_.daughter1_d0);    
    tree_->Branch("reco_z_from_muons_daughter1_dxy",       &reco_z_from_muons_.daughter1_dxy);    
    tree_->Branch("reco_z_from_muons_daughter1_dz",        &reco_z_from_muons_.daughter1_dz);    
    tree_->Branch("reco_z_from_muons_daughter1_d0err",     &reco_z_from_muons_.daughter1_d0err);    
    tree_->Branch("reco_z_from_muons_daughter1_dxyerr",    &reco_z_from_muons_.daughter1_dxyerr);    
    tree_->Branch("reco_z_from_muons_daughter1_dzerr",     &reco_z_from_muons_.daughter1_dzerr);    
    tree_->Branch("reco_z_from_muons_daughter1_trkKink",   &reco_z_from_muons_.daughter1_trkKink);    
    tree_->Branch("reco_z_from_muons_daughter1_glbKink",   &reco_z_from_muons_.daughter1_glbKink);    
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
    tree_->Branch("reco_jpsi_muon0_d0",                &reco_jpsi_.muon0_d0);
    tree_->Branch("reco_jpsi_muon0_dxy",               &reco_jpsi_.muon0_dxy);
    tree_->Branch("reco_jpsi_muon0_dz",                &reco_jpsi_.muon0_dz);
    tree_->Branch("reco_jpsi_muon0_d0err",             &reco_jpsi_.muon0_d0err);
    tree_->Branch("reco_jpsi_muon0_dxyerr",            &reco_jpsi_.muon0_dxyerr);
    tree_->Branch("reco_jpsi_muon0_dzerr",             &reco_jpsi_.muon0_dzerr);
    tree_->Branch("reco_jpsi_muons0_trkKink",          &reco_jpsi_.muon0_trkKink);    
    tree_->Branch("reco_jpsi_muons0_glbKink",          &reco_jpsi_.muon0_glbKink);    
    tree_->Branch("reco_jpsi_muon1_pt",                &reco_jpsi_.muon1_pt);
    tree_->Branch("reco_jpsi_muon1_eta",               &reco_jpsi_.muon1_eta);
    tree_->Branch("reco_jpsi_muon1_phi",               &reco_jpsi_.muon1_phi);
    tree_->Branch("reco_jpsi_muon1_d0",                &reco_jpsi_.muon1_d0);
    tree_->Branch("reco_jpsi_muon1_dxy",               &reco_jpsi_.muon1_dxy);
    tree_->Branch("reco_jpsi_muon1_dz",                &reco_jpsi_.muon1_dz);
    tree_->Branch("reco_jpsi_muon1_d0err",             &reco_jpsi_.muon1_d0err);
    tree_->Branch("reco_jpsi_muon1_dxyerr",            &reco_jpsi_.muon1_dxyerr);
    tree_->Branch("reco_jpsi_muon1_dzerr",             &reco_jpsi_.muon1_dzerr);
    tree_->Branch("reco_jpsi_muons1_trkKink",          &reco_jpsi_.muon1_trkKink);    
    tree_->Branch("reco_jpsi_muons1_glbKink",          &reco_jpsi_.muon1_glbKink);    
    tree_->Branch("reco_jpsi_muon0_charge",            &reco_jpsi_.muon0_charge);
    tree_->Branch("reco_jpsi_muon1_charge",            &reco_jpsi_.muon1_charge);
    tree_->Branch("reco_jpsi_muon0_iso_sum_charged_hadron_pt",   &reco_jpsi_.muon0_iso_sum_charged_hadron_pt);
    tree_->Branch("reco_jpsi_muon0_iso_sum_charged_particle_pt", &reco_jpsi_.muon0_iso_sum_charged_particle_pt);
    tree_->Branch("reco_jpsi_muon0_iso_sum_neutral_hadron_et",   &reco_jpsi_.muon0_iso_sum_neutral_hadron_et);
    tree_->Branch("reco_jpsi_muon0_iso_sum_photon_et",           &reco_jpsi_.muon0_iso_sum_photon_et);
    tree_->Branch("reco_jpsi_muon0_iso_sum_pileup_pt",           &reco_jpsi_.muon0_iso_sum_pileup_pt);
    tree_->Branch("reco_jpsi_muon0_iso",                         &reco_jpsi_.muon0_iso);
    tree_->Branch("reco_jpsi_muon1_iso_sum_charged_hadron_pt",   &reco_jpsi_.muon1_iso_sum_charged_hadron_pt); 
    tree_->Branch("reco_jpsi_muon1_iso_sum_charged_particle_pt", &reco_jpsi_.muon1_iso_sum_charged_particle_pt);
    tree_->Branch("reco_jpsi_muon1_iso_sum_neutral_hadron_et",   &reco_jpsi_.muon1_iso_sum_neutral_hadron_et);
    tree_->Branch("reco_jpsi_muon1_iso_sum_photon_et",           &reco_jpsi_.muon1_iso_sum_photon_et);
    tree_->Branch("reco_jpsi_muon1_iso_sum_pileup_pt",           &reco_jpsi_.muon1_iso_sum_pileup_pt);
    tree_->Branch("reco_jpsi_muon1_iso",                         &reco_jpsi_.muon1_iso);
    tree_->Branch("reco_jpsi_muon0_charge",                          &reco_jpsi_.muon0_charge);
    tree_->Branch("reco_jpsi_muon1_charge",                          &reco_jpsi_.muon1_charge);
    tree_->Branch("reco_jpsi_has_muons_in_eta_window", &reco_jpsi_.has_muons_in_eta_window);
    tree_->Branch("reco_jpsi_has_high_pt_muons",       &reco_jpsi_.has_high_pt_muons);
    tree_->Branch("reco_jpsi_soft_id",                 &reco_jpsi_.reco_jpsi_soft_id); 

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
    tree_->Branch("reco_jpsi_from_electrons_soft_id",                 &reco_jpsi_from_electrons_.reco_jpsi_soft_id); 
    tree_->Branch("triggers_passed",                                  &reco_z_from_muons_.passed_triggers);

    tree_->Branch("four_lepton_vertex_muon0_pt",  &four_lepton_vertex_.muon0_pt);
    tree_->Branch("four_lepton_vertex_muon1_pt",  &four_lepton_vertex_.muon1_pt);
    tree_->Branch("four_lepton_vertex_muon2_pt",  &four_lepton_vertex_.muon2_pt);
    tree_->Branch("four_lepton_vertex_muon3_pt",  &four_lepton_vertex_.muon3_pt);
    tree_->Branch("four_lepton_vertex_muon0_eta", &four_lepton_vertex_.muon0_eta);
    tree_->Branch("four_lepton_vertex_muon1_eta", &four_lepton_vertex_.muon1_eta);
    tree_->Branch("four_lepton_vertex_muon2_eta", &four_lepton_vertex_.muon2_eta);
    tree_->Branch("four_lepton_vertex_muon3_eta", &four_lepton_vertex_.muon3_eta);
    tree_->Branch("four_lepton_vertex_muon0_phi", &four_lepton_vertex_.muon0_phi);
    tree_->Branch("four_lepton_vertex_muon1_phi", &four_lepton_vertex_.muon1_phi);
    tree_->Branch("four_lepton_vertex_muon2_phi", &four_lepton_vertex_.muon2_phi);
    tree_->Branch("four_lepton_vertex_muon3_phi", &four_lepton_vertex_.muon3_phi);
    tree_->Branch("four_lepton_vertex_vtx_chi2",  &four_lepton_vertex_.vtx_chi2);
    tree_->Branch("four_lepton_vertex_vtx_ndf",   &four_lepton_vertex_.vtx_ndf);
    tree_->Branch("four_lepton_vertex_vtx_prob",  &four_lepton_vertex_.vtx_prob);

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

      tree_->Branch("truth_pdg_id_all", &truth_z_muons_.truth_pdg_id);
    }

    //-----------------
    tree_->Branch("event_weight",                                                                &event_.event_weight);
    tree_->Branch("event_number",                                                                &event_.event_number);
    tree_->Branch("run_number",                                                                  &event_.run_number);  
    tree_->Branch("n_verts",                                                                     &event_.n_verts);     
    tree_->Branch("found_good_muons_from_z",                                                     &event_.found_good_muons_from_z);   
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
    four_lepton_vertex_.clear_values();

    // Reco
    
    if (is_Jpsimumu) {
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
//      for (unsigned int iZ = 0; iZ < zfe.reco_z_from_muons.muon0_pT.size() ; ++iZ ) {
//        math::PtEtaPhiMLorentzVector Z_l1(zfe.reco_z_from_muons.muon0_pT.at(iZ), 
//                                          zfe.reco_z_from_muons.muon0_eta.at(iZ),
//                                          zfe.reco_z_from_muons.muon0_phi.at(iZ),
//                                          0.1056583715);
//        math::PtEtaPhiMLorentzVector Z_l2(zfe.reco_z_from_muons.muon1_pT.at(iZ), 
//                                          zfe.reco_z_from_muons.muon1_eta.at(iZ),
//                                          zfe.reco_z_from_muons.muon1_phi.at(iZ),
//                                          0.1056583715);

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

          if (Z_l1.pt()==Jpsi_l1.pt() || Z_l1.pt()==Jpsi_l2.pt() || Z_l2.pt()==Jpsi_l1.pt() || Z_l2.pt()==Jpsi_l2.pt())
            continue;

          // Z->ee
          if (is_Zee) { 
            reco_z_.z_m              .push_back((Z_l1+Z_l2).mass());
            reco_z_.z_pt             .push_back((Z_l1+Z_l2).pt());
            reco_z_.z_y              .push_back(zfe.reco_z_from_muons.y);
            reco_z_.z_phi            .push_back(zfe.reco_z_from_muons.phi);
            reco_z_.z_phistar        .push_back(zfe.reco_z_from_muons.phistar);
            reco_z_.z_eta            .push_back((Z_l1+Z_l2).eta());
            reco_z_.z_vtx_prob       .push_back(zfe.reco_z.vtx_prob.at(iZ));
            reco_z_.z_vtx_x          .push_back(zfe.reco_z.vtx_x.at(iZ));
            reco_z_.z_vtx_y          .push_back(zfe.reco_z.vtx_y.at(iZ));
            reco_z_.z_vtx_z          .push_back(zfe.reco_z.vtx_z.at(iZ));
            reco_z_.daughter0_pt     .push_back(zfe.reco_z.muon0_pT.at(iZ));
            reco_z_.daughter0_eta    .push_back(zfe.reco_z.muon0_eta.at(iZ));
            reco_z_.daughter0_phi    .push_back(zfe.reco_z.muon0_phi.at(iZ));
            reco_z_.daughter0_d0     .push_back(zfe.reco_z.muon0_d0.at(iZ)); 
            reco_z_.daughter0_dxy    .push_back(zfe.reco_z.muon0_dxy.at(iZ));
            reco_z_.daughter0_dz     .push_back(zfe.reco_z.muon0_dz.at(iZ));
            reco_z_.daughter0_d0err  .push_back(zfe.reco_z.muon0_d0err.at(iZ)); 
            reco_z_.daughter0_dxyerr .push_back(zfe.reco_z.muon0_dxyerr.at(iZ));
            reco_z_.daughter0_dzerr  .push_back(zfe.reco_z.muon0_dzerr.at(iZ));
            reco_z_.daughter1_pt     .push_back(zfe.reco_z.muon1_pT.at(iZ));
            reco_z_.daughter1_eta    .push_back(zfe.reco_z.muon1_eta.at(iZ));
            reco_z_.daughter1_phi    .push_back(zfe.reco_z.muon1_phi.at(iZ));
            reco_z_.daughter1_d0     .push_back(zfe.reco_z.muon1_d0.at(iZ)); 
            reco_z_.daughter1_dxy    .push_back(zfe.reco_z.muon1_dxy.at(iZ));
            reco_z_.daughter1_dz     .push_back(zfe.reco_z.muon1_dz.at(iZ));
            reco_z_.daughter1_d0err  .push_back(zfe.reco_z.muon1_d0err.at(iZ)); 
            reco_z_.daughter1_dxyerr .push_back(zfe.reco_z.muon1_dxyerr.at(iZ));
            reco_z_.daughter1_dzerr  .push_back(zfe.reco_z.muon1_dzerr.at(iZ));
          }

          // z->mumu
          if (is_Zmumu) {
            reco_z_from_muons_.z_m        .push_back((Z_l1+Z_l2).mass());
            reco_z_from_muons_.z_pt       .push_back((Z_l1+Z_l2).pt());
            reco_z_from_muons_.z_eta      .push_back((Z_l1+Z_l2).eta());
            reco_z_from_muons_.z_y        .push_back(zfe.reco_z_from_muons.y);
            reco_z_from_muons_.z_phi      .push_back(zfe.reco_z_from_muons.phi);
            reco_z_from_muons_.z_phistar  .push_back(zfe.reco_z_from_muons.phistar);
            reco_z_from_muons_.z_vtx_prob .push_back(zfe.reco_z_from_muons.vtx_prob.at(iZ));
            reco_z_from_muons_.z_vtx_x    .push_back(zfe.reco_z_from_muons.vtx_x.at(iZ));
            reco_z_from_muons_.z_vtx_y    .push_back(zfe.reco_z_from_muons.vtx_y.at(iZ));
            reco_z_from_muons_.z_vtx_z    .push_back(zfe.reco_z_from_muons.vtx_z.at(iZ));
            reco_z_from_muons_.daughter0_pt     .push_back(zfe.reco_z_from_muons.muon0_pT.at(iZ));
            reco_z_from_muons_.daughter0_eta    .push_back(zfe.reco_z_from_muons.muon0_eta.at(iZ));
            reco_z_from_muons_.daughter0_phi    .push_back(zfe.reco_z_from_muons.muon0_phi.at(iZ));
            reco_z_from_muons_.daughter0_d0     .push_back(zfe.reco_z_from_muons.muon0_d0.at(iZ)); 
            reco_z_from_muons_.daughter0_dxy    .push_back(zfe.reco_z_from_muons.muon0_dxy.at(iZ));
            reco_z_from_muons_.daughter0_dz     .push_back(zfe.reco_z_from_muons.muon0_dz.at(iZ));
            reco_z_from_muons_.daughter0_d0err  .push_back(zfe.reco_z_from_muons.muon0_d0err.at(iZ)); 
            reco_z_from_muons_.daughter0_dxyerr .push_back(zfe.reco_z_from_muons.muon0_dxyerr.at(iZ));
            reco_z_from_muons_.daughter0_dzerr  .push_back(zfe.reco_z_from_muons.muon0_dzerr.at(iZ));
            reco_z_from_muons_.daughter0_trkKink.push_back(zfe.reco_z_from_muons.muon0_trkKink.at(iZ));
            reco_z_from_muons_.daughter0_glbKink.push_back(zfe.reco_z_from_muons.muon0_glbKink.at(iZ));
            reco_z_from_muons_.daughter1_pt     .push_back(zfe.reco_z_from_muons.muon1_pT.at(iZ));
            reco_z_from_muons_.daughter1_eta    .push_back(zfe.reco_z_from_muons.muon1_eta.at(iZ));
            reco_z_from_muons_.daughter1_phi    .push_back(zfe.reco_z_from_muons.muon1_phi.at(iZ));
            reco_z_from_muons_.daughter1_d0     .push_back(zfe.reco_z_from_muons.muon1_d0.at(iZ)); 
            reco_z_from_muons_.daughter1_dxy    .push_back(zfe.reco_z_from_muons.muon1_dxy.at(iZ));
            reco_z_from_muons_.daughter1_dz     .push_back(zfe.reco_z_from_muons.muon1_dz.at(iZ));
            reco_z_from_muons_.daughter1_d0err  .push_back(zfe.reco_z_from_muons.muon1_d0err.at(iZ)); 
            reco_z_from_muons_.daughter1_dxyerr .push_back(zfe.reco_z_from_muons.muon1_dxyerr.at(iZ));
            reco_z_from_muons_.daughter1_dzerr  .push_back(zfe.reco_z_from_muons.muon1_dzerr.at(iZ));

            reco_z_from_muons_.daughter1_trkKink.push_back(zfe.reco_z_from_muons.muon1_trkKink.at(iZ));
            reco_z_from_muons_.daughter1_glbKink.push_back(zfe.reco_z_from_muons.muon1_glbKink.at(iZ));
          }

          // 4-lepton vertex
          for (unsigned int ifound = 0; ifound < zfe.four_lepton_vertex.muon0_pt.size(); ++ifound) {
              if ((zfe.four_lepton_vertex.muon0_pt.at(ifound) == Z_l1.pt() ||
                   zfe.four_lepton_vertex.muon0_pt.at(ifound) == Z_l2.pt() ||
                   zfe.four_lepton_vertex.muon0_pt.at(ifound) == Jpsi_l1.pt() ||
                   zfe.four_lepton_vertex.muon0_pt.at(ifound) == Jpsi_l2.pt()) &&
                  (zfe.four_lepton_vertex.muon1_pt.at(ifound) == Z_l1.pt() ||
                   zfe.four_lepton_vertex.muon1_pt.at(ifound) == Z_l2.pt() ||
                   zfe.four_lepton_vertex.muon1_pt.at(ifound) == Jpsi_l1.pt() ||
                   zfe.four_lepton_vertex.muon1_pt.at(ifound) == Jpsi_l2.pt()) &&
                  (zfe.four_lepton_vertex.muon2_pt.at(ifound) == Z_l1.pt() ||
                   zfe.four_lepton_vertex.muon2_pt.at(ifound) == Z_l2.pt() ||
                   zfe.four_lepton_vertex.muon2_pt.at(ifound) == Jpsi_l1.pt() ||
                   zfe.four_lepton_vertex.muon2_pt.at(ifound) == Jpsi_l2.pt()) &&
                  (zfe.four_lepton_vertex.muon3_pt.at(ifound) == Z_l1.pt() ||
                   zfe.four_lepton_vertex.muon3_pt.at(ifound) == Z_l2.pt() ||
                   zfe.four_lepton_vertex.muon3_pt.at(ifound) == Jpsi_l1.pt() ||
                   zfe.four_lepton_vertex.muon3_pt.at(ifound) == Jpsi_l2.pt())) {                  
              four_lepton_vertex_.muon0_pt .push_back(zfe.four_lepton_vertex.muon0_pt.at(ifound)); 
              four_lepton_vertex_.muon1_pt .push_back(zfe.four_lepton_vertex.muon1_pt.at(ifound)); 
              four_lepton_vertex_.muon2_pt .push_back(zfe.four_lepton_vertex.muon2_pt.at(ifound)); 
              four_lepton_vertex_.muon3_pt .push_back(zfe.four_lepton_vertex.muon3_pt.at(ifound)); 
              four_lepton_vertex_.muon0_eta.push_back(zfe.four_lepton_vertex.muon0_pt.at(ifound));
              four_lepton_vertex_.muon1_eta.push_back(zfe.four_lepton_vertex.muon1_pt.at(ifound));
              four_lepton_vertex_.muon2_eta.push_back(zfe.four_lepton_vertex.muon2_pt.at(ifound));
              four_lepton_vertex_.muon3_eta.push_back(zfe.four_lepton_vertex.muon3_pt.at(ifound));
              four_lepton_vertex_.muon0_phi.push_back(zfe.four_lepton_vertex.muon0_pt.at(ifound));
              four_lepton_vertex_.muon1_phi.push_back(zfe.four_lepton_vertex.muon1_pt.at(ifound));
              four_lepton_vertex_.muon2_phi.push_back(zfe.four_lepton_vertex.muon2_pt.at(ifound));
              four_lepton_vertex_.muon3_phi.push_back(zfe.four_lepton_vertex.muon3_pt.at(ifound));
              four_lepton_vertex_.vtx_chi2 .push_back(zfe.four_lepton_vertex.vtx_chi2.at(ifound)); 
              four_lepton_vertex_.vtx_ndf  .push_back(zfe.four_lepton_vertex.vtx_ndf.at(ifound)); 
              four_lepton_vertex_.vtx_prob .push_back(zfe.four_lepton_vertex.vtx_prob.at(ifound)); 
            }
          }
          // end here

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
          reco_jpsi_.jpsi_vtx_prob           .push_back(zfe.reco_jpsi.vtx_prob.at(i));
          reco_jpsi_.jpsi_vtx_x              .push_back(zfe.reco_jpsi.vtx_x.at(i));
          reco_jpsi_.jpsi_vtx_y              .push_back(zfe.reco_jpsi.vtx_y.at(i));
          reco_jpsi_.jpsi_vtx_z              .push_back(zfe.reco_jpsi.vtx_z.at(i));
          reco_jpsi_.jpsi_eff                .push_back(zfe.reco_jpsi.jpsi_efficiency.at(i));
          reco_jpsi_.jpsi_acc_eff            .push_back(zfe.reco_jpsi.jpsi_acc_eff.at(i));
          reco_jpsi_.jpsi_scale_factor       .push_back(zfe.reco_jpsi.jpsi_scale_factor.at(i));
          reco_jpsi_.muon0_pt                .push_back(zfe.reco_jpsi.muon0.at(i).pt());
          reco_jpsi_.muon0_eta               .push_back(zfe.reco_jpsi.muon0.at(i).eta());
          reco_jpsi_.muon0_phi               .push_back(zfe.reco_jpsi.muon0.at(i).phi());
          reco_jpsi_.muon0_d0                .push_back(zfe.reco_jpsi.muon0_d0.at(i));
          reco_jpsi_.muon0_dxy               .push_back(zfe.reco_jpsi.muon0_dxy.at(i));
          reco_jpsi_.muon0_dz                .push_back(zfe.reco_jpsi.muon0_dz.at(i));
          reco_jpsi_.muon0_d0err             .push_back(zfe.reco_jpsi.muon0_d0err.at(i));
          reco_jpsi_.muon0_dxyerr            .push_back(zfe.reco_jpsi.muon0_dxyerr.at(i));
          reco_jpsi_.muon0_dzerr             .push_back(zfe.reco_jpsi.muon0_dzerr.at(i));
          reco_jpsi_.muon0_trkKink           .push_back(zfe.reco_jpsi.muon0.at(i).combinedQuality().trkKink);
          reco_jpsi_.muon0_glbKink           .push_back(zfe.reco_jpsi.muon0.at(i).combinedQuality().glbKink);
          reco_jpsi_.muon1_pt                .push_back(zfe.reco_jpsi.muon1.at(i).pt());
          reco_jpsi_.muon1_eta               .push_back(zfe.reco_jpsi.muon1.at(i).eta());
          reco_jpsi_.muon1_phi               .push_back(zfe.reco_jpsi.muon1.at(i).phi());
          reco_jpsi_.muon1_d0                .push_back(zfe.reco_jpsi.muon1_d0.at(i));
          reco_jpsi_.muon1_dxy               .push_back(zfe.reco_jpsi.muon1_dxy.at(i));
          reco_jpsi_.muon1_dz                .push_back(zfe.reco_jpsi.muon1_dz.at(i));
          reco_jpsi_.muon1_d0err             .push_back(zfe.reco_jpsi.muon1_d0err.at(i));
          reco_jpsi_.muon1_dxyerr            .push_back(zfe.reco_jpsi.muon1_dxyerr.at(i));
          reco_jpsi_.muon1_dzerr             .push_back(zfe.reco_jpsi.muon1_dzerr.at(i));
          reco_jpsi_.muon1_trkKink           .push_back(zfe.reco_jpsi.muon1.at(i).combinedQuality().trkKink);
          reco_jpsi_.muon1_glbKink           .push_back(zfe.reco_jpsi.muon1.at(i).combinedQuality().glbKink);
          reco_jpsi_.muon0_charge            .push_back(zfe.reco_jpsi.muon0.at(i).charge());
          reco_jpsi_.muon1_charge            .push_back(zfe.reco_jpsi.muon1.at(i).charge());
          reco_jpsi_.muon0_iso_sum_charged_hadron_pt     .push_back(zfe.reco_jpsi.iso_sum_charged_hadron_pt_mu0.at(i));
          reco_jpsi_.muon0_iso_sum_charged_particle_pt   .push_back(zfe.reco_jpsi.iso_sum_charged_particle_pt_mu0.at(i));
          reco_jpsi_.muon0_iso_sum_neutral_hadron_et     .push_back(zfe.reco_jpsi.iso_sum_neutral_hadron_et_mu0.at(i));   
          reco_jpsi_.muon0_iso_sum_photon_et             .push_back(zfe.reco_jpsi.iso_sum_photon_et_mu0.at(i));           
          reco_jpsi_.muon0_iso_sum_pileup_pt             .push_back(zfe.reco_jpsi.iso_sum_pileup_pt_mu0.at(i));         
          reco_jpsi_.muon0_iso                           .push_back(zfe.reco_jpsi.iso_mu0.at(i));                       
          reco_jpsi_.muon1_iso_sum_charged_hadron_pt     .push_back(zfe.reco_jpsi.iso_sum_charged_hadron_pt_mu1.at(i));   
          reco_jpsi_.muon1_iso_sum_charged_particle_pt   .push_back(zfe.reco_jpsi.iso_sum_charged_particle_pt_mu1.at(i)); 
          reco_jpsi_.muon1_iso_sum_neutral_hadron_et     .push_back(zfe.reco_jpsi.iso_sum_neutral_hadron_et_mu1.at(i));   
          reco_jpsi_.muon1_iso_sum_photon_et             .push_back(zfe.reco_jpsi.iso_sum_photon_et_mu1.at(i));           
          reco_jpsi_.muon1_iso_sum_pileup_pt             .push_back(zfe.reco_jpsi.iso_sum_pileup_pt_mu1.at(i));           
          reco_jpsi_.muon1_iso                           .push_back(zfe.reco_jpsi.iso_mu1.at(i));                         

          reco_jpsi_.has_muons_in_eta_window .push_back(zfe.reco_jpsi.has_muons_in_eta_window.at(i));
          reco_jpsi_.has_high_pt_muons       .push_back(zfe.reco_jpsi.has_high_pt_muons.at(i));
          reco_jpsi_.reco_jpsi_soft_id       .push_back(zfe.reco_jpsi.has_soft_id_muons.at(i));

          for (unsigned int iTrig = 0; iTrig < zfe.reco_z_from_muons.trigger_list.size(); ++iTrig)
            reco_z_from_muons_.passed_triggers.push_back(zfe.reco_z_from_muons.trigger_list.at(iTrig));

          // event info
          event_.event_weight                                                .push_back(zfe.event_weight);
          event_.event_number                                                .push_back(zfe.id.event_num);
          event_.run_number                                                  .push_back(zfe.id.run_num);
          event_.n_verts                                                     .push_back(zfe.reco_vert.num);
          event_.found_good_muons_from_z                                     .push_back(zfe.found_good_muons_from_z);
          event_.found_good_electrons_from_z                                 .push_back(zfe.found_good_electrons_from_z);

        }
      }
    }

    //---------------------------------------------------------------------------------------
    if (is_Jpsiee) {
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
//      for (unsigned int iZ = 0; iZ < zfe.reco_z_from_muons.muon0_pT.size() ; ++iZ ) {
//        math::PtEtaPhiMLorentzVector Z_l1(zfe.reco_z_from_muons.muon0_pT.at(iZ), 
//                                          zfe.reco_z_from_muons.muon0_eta.at(iZ),
//                                          zfe.reco_z_from_muons.muon0_phi.at(iZ),
//                                          0.1056583715);
//        math::PtEtaPhiMLorentzVector Z_l2(zfe.reco_z_from_muons.muon1_pT.at(iZ), 
//                                          zfe.reco_z_from_muons.muon1_eta.at(iZ),
//                                          zfe.reco_z_from_muons.muon1_phi.at(iZ),
//                                          0.1056583715);
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

          if (Z_l1.pt()==Jpsi_l1.pt() || Z_l1.pt()==Jpsi_l2.pt() || Z_l2.pt()==Jpsi_l1.pt() || Z_l2.pt()==Jpsi_l2.pt())
            continue;

          // 4-lepton vertex
          for (unsigned int ifound = 0; ifound < zfe.four_lepton_vertex.muon0_pt.size(); ++ifound) {
            if ((zfe.four_lepton_vertex.muon0_pt.at(ifound) == Z_l1.pt() ||
                 zfe.four_lepton_vertex.muon0_pt.at(ifound) == Z_l2.pt() ||
                 zfe.four_lepton_vertex.muon0_pt.at(ifound) == Jpsi_l1.pt() ||
                 zfe.four_lepton_vertex.muon0_pt.at(ifound) == Jpsi_l2.pt()) &&
                (zfe.four_lepton_vertex.muon1_pt.at(ifound) == Z_l1.pt() ||
                 zfe.four_lepton_vertex.muon1_pt.at(ifound) == Z_l2.pt() ||
                 zfe.four_lepton_vertex.muon1_pt.at(ifound) == Jpsi_l1.pt() ||
                 zfe.four_lepton_vertex.muon1_pt.at(ifound) == Jpsi_l2.pt()) &&
                (zfe.four_lepton_vertex.muon2_pt.at(ifound) == Z_l1.pt() ||
                 zfe.four_lepton_vertex.muon2_pt.at(ifound) == Z_l2.pt() ||
                 zfe.four_lepton_vertex.muon2_pt.at(ifound) == Jpsi_l1.pt() ||
                 zfe.four_lepton_vertex.muon2_pt.at(ifound) == Jpsi_l2.pt()) &&
                (zfe.four_lepton_vertex.muon3_pt.at(ifound) == Z_l1.pt() ||
                 zfe.four_lepton_vertex.muon3_pt.at(ifound) == Z_l2.pt() ||
                 zfe.four_lepton_vertex.muon3_pt.at(ifound) == Jpsi_l1.pt() ||
                 zfe.four_lepton_vertex.muon3_pt.at(ifound) == Jpsi_l2.pt())
            ) {                  
              four_lepton_vertex_.muon0_pt .push_back(zfe.four_lepton_vertex.muon0_pt.at(ifound)); 
              four_lepton_vertex_.muon1_pt .push_back(zfe.four_lepton_vertex.muon1_pt.at(ifound)); 
              four_lepton_vertex_.muon2_pt .push_back(zfe.four_lepton_vertex.muon2_pt.at(ifound)); 
              four_lepton_vertex_.muon3_pt .push_back(zfe.four_lepton_vertex.muon3_pt.at(ifound)); 
              four_lepton_vertex_.muon0_eta.push_back(zfe.four_lepton_vertex.muon0_pt.at(ifound));
              four_lepton_vertex_.muon1_eta.push_back(zfe.four_lepton_vertex.muon1_pt.at(ifound));
              four_lepton_vertex_.muon2_eta.push_back(zfe.four_lepton_vertex.muon2_pt.at(ifound));
              four_lepton_vertex_.muon3_eta.push_back(zfe.four_lepton_vertex.muon3_pt.at(ifound));
              four_lepton_vertex_.muon0_phi.push_back(zfe.four_lepton_vertex.muon0_pt.at(ifound));
              four_lepton_vertex_.muon1_phi.push_back(zfe.four_lepton_vertex.muon1_pt.at(ifound));
              four_lepton_vertex_.muon2_phi.push_back(zfe.four_lepton_vertex.muon2_pt.at(ifound));
              four_lepton_vertex_.muon3_phi.push_back(zfe.four_lepton_vertex.muon3_pt.at(ifound));
              four_lepton_vertex_.vtx_chi2 .push_back(zfe.four_lepton_vertex.vtx_chi2.at(ifound)); 
              four_lepton_vertex_.vtx_ndf  .push_back(zfe.four_lepton_vertex.vtx_ndf.at(ifound)); 
              four_lepton_vertex_.vtx_prob .push_back(zfe.four_lepton_vertex.vtx_prob.at(ifound)); 
            }
          }
          // end here

          // Z->ee
          if (is_Zee) {
            reco_z_.z_m              .push_back((Z_l1+Z_l2).mass());
            reco_z_.z_pt             .push_back((Z_l1+Z_l2).pt());
            reco_z_.z_eta            .push_back((Z_l1+Z_l2).eta());
            reco_z_.z_y              .push_back(zfe.reco_z.y);
            reco_z_.z_phi            .push_back(zfe.reco_z.phi);
            reco_z_.z_phistar        .push_back(zfe.reco_z.phistar);
            reco_z_.z_vtx_prob       .push_back(zfe.reco_z.vtx_prob.at(iZ));
            reco_z_.z_vtx_x          .push_back(zfe.reco_z.vtx_x.at(iZ));
            reco_z_.z_vtx_y          .push_back(zfe.reco_z.vtx_y.at(iZ));
            reco_z_.z_vtx_z          .push_back(zfe.reco_z.vtx_z.at(iZ));
            reco_z_.daughter0_pt     .push_back(zfe.reco_z.muon0_pT.at(iZ));
            reco_z_.daughter0_eta    .push_back(zfe.reco_z.muon0_eta.at(iZ));
            reco_z_.daughter0_phi    .push_back(zfe.reco_z.muon0_phi.at(iZ));
            reco_z_.daughter0_d0     .push_back(zfe.reco_z.muon0_d0.at(iZ)); 
            reco_z_.daughter0_dxy    .push_back(zfe.reco_z.muon0_dxy.at(iZ));
            reco_z_.daughter0_dz     .push_back(zfe.reco_z.muon0_dz.at(iZ));
            reco_z_.daughter0_d0err  .push_back(zfe.reco_z.muon0_d0err.at(iZ)); 
            reco_z_.daughter0_dxyerr .push_back(zfe.reco_z.muon0_dxyerr.at(iZ));
            reco_z_.daughter0_dzerr  .push_back(zfe.reco_z.muon0_dzerr.at(iZ));
            reco_z_.daughter1_pt     .push_back(zfe.reco_z.muon1_pT.at(iZ));
            reco_z_.daughter1_eta    .push_back(zfe.reco_z.muon1_eta.at(iZ));
            reco_z_.daughter1_phi    .push_back(zfe.reco_z.muon1_phi.at(iZ));
            reco_z_.daughter1_d0     .push_back(zfe.reco_z.muon1_d0.at(iZ)); 
            reco_z_.daughter1_dxy    .push_back(zfe.reco_z.muon1_dxy.at(iZ));
            reco_z_.daughter1_dz     .push_back(zfe.reco_z.muon1_dz.at(iZ));
            reco_z_.daughter1_d0err  .push_back(zfe.reco_z.muon1_d0err.at(iZ)); 
            reco_z_.daughter1_dxyerr .push_back(zfe.reco_z.muon1_dxyerr.at(iZ));
            reco_z_.daughter1_dzerr  .push_back(zfe.reco_z.muon1_dzerr.at(iZ));
          }

          // z->mumu
          if (is_Zmumu) {
            reco_z_from_muons_.z_m        .push_back((Z_l1+Z_l2).mass());
            reco_z_from_muons_.z_pt       .push_back((Z_l1+Z_l2).pt());
            reco_z_from_muons_.z_y        .push_back(zfe.reco_z_from_muons.y);
            reco_z_from_muons_.z_phi      .push_back(zfe.reco_z_from_muons.phi);
            reco_z_from_muons_.z_phistar  .push_back(zfe.reco_z_from_muons.phistar);
            reco_z_from_muons_.z_eta      .push_back((Z_l1+Z_l2).eta());
            reco_z_from_muons_.z_vtx_prob .push_back(zfe.reco_z_from_muons.vtx_prob.at(iZ));
            reco_z_from_muons_.z_vtx_x    .push_back(zfe.reco_z_from_muons.vtx_x.at(iZ));
            reco_z_from_muons_.z_vtx_y    .push_back(zfe.reco_z_from_muons.vtx_y.at(iZ));
            reco_z_from_muons_.z_vtx_z    .push_back(zfe.reco_z_from_muons.vtx_z.at(iZ));
            reco_z_from_muons_.daughter0_pt     .push_back(zfe.reco_z_from_muons.muon0_pT.at(iZ));
            reco_z_from_muons_.daughter0_eta    .push_back(zfe.reco_z_from_muons.muon0_eta.at(iZ));
            reco_z_from_muons_.daughter0_phi    .push_back(zfe.reco_z_from_muons.muon0_phi.at(iZ));
            reco_z_from_muons_.daughter0_d0     .push_back(zfe.reco_z_from_muons.muon0_d0.at(iZ)); 
            reco_z_from_muons_.daughter0_dxy    .push_back(zfe.reco_z_from_muons.muon0_dxy.at(iZ));
            reco_z_from_muons_.daughter0_dz     .push_back(zfe.reco_z_from_muons.muon0_dz.at(iZ));
            reco_z_from_muons_.daughter0_d0err  .push_back(zfe.reco_z_from_muons.muon0_d0err.at(iZ)); 
            reco_z_from_muons_.daughter0_dxyerr .push_back(zfe.reco_z_from_muons.muon0_dxyerr.at(iZ));
            reco_z_from_muons_.daughter0_dzerr  .push_back(zfe.reco_z_from_muons.muon0_dzerr.at(iZ));
            reco_z_from_muons_.daughter0_trkKink.push_back(zfe.reco_z_from_muons.muon0_trkKink.at(iZ));
            reco_z_from_muons_.daughter0_glbKink.push_back(zfe.reco_z_from_muons.muon0_glbKink.at(iZ));
            reco_z_from_muons_.daughter1_pt     .push_back(zfe.reco_z_from_muons.muon1_pT.at(iZ));
            reco_z_from_muons_.daughter1_eta    .push_back(zfe.reco_z_from_muons.muon1_eta.at(iZ));
            reco_z_from_muons_.daughter1_phi    .push_back(zfe.reco_z_from_muons.muon1_phi.at(iZ));
            reco_z_from_muons_.daughter1_d0     .push_back(zfe.reco_z_from_muons.muon1_d0.at(iZ)); 
            reco_z_from_muons_.daughter1_dxy    .push_back(zfe.reco_z_from_muons.muon1_dxy.at(iZ));
            reco_z_from_muons_.daughter1_dz     .push_back(zfe.reco_z_from_muons.muon1_dz.at(iZ));
            reco_z_from_muons_.daughter1_d0err  .push_back(zfe.reco_z_from_muons.muon1_d0err.at(iZ)); 
            reco_z_from_muons_.daughter1_dxyerr .push_back(zfe.reco_z_from_muons.muon1_dxyerr.at(iZ));
            reco_z_from_muons_.daughter1_dzerr  .push_back(zfe.reco_z_from_muons.muon1_dzerr.at(iZ));
            reco_z_from_muons_.daughter1_trkKink.push_back(zfe.reco_z_from_muons.muon1_trkKink.at(iZ));
            reco_z_from_muons_.daughter1_glbKink.push_back(zfe.reco_z_from_muons.muon1_glbKink.at(iZ));
            reco_z_from_muons_.daughter0_charge .push_back(zfe.z_muon0.charge());
            reco_z_from_muons_.daughter1_charge .push_back(zfe.z_muon1.charge());
          }

          // jpsi->ee
          reco_jpsi_from_electrons_.jpsi_m                  .push_back(zfe.reco_jpsi_from_electrons.m.at(i));
          reco_jpsi_from_electrons_.jpsi_pt                 .push_back(zfe.reco_jpsi_from_electrons.pt.at(i));
          reco_jpsi_from_electrons_.jpsi_y                  .push_back(zfe.reco_jpsi_from_electrons.y.at(i));
          reco_jpsi_from_electrons_.jpsi_phi                .push_back(zfe.reco_jpsi_from_electrons.phi.at(i));
          reco_jpsi_from_electrons_.jpsi_eta                .push_back(zfe.reco_jpsi_from_electrons.eta.at(i));
          reco_jpsi_from_electrons_.jpsi_tau_xy             .push_back(zfe.reco_jpsi_from_electrons.tau_xy.at(i));
          reco_jpsi_from_electrons_.jpsi_tau_z              .push_back(zfe.reco_jpsi_from_electrons.tau_z.at(i));
          reco_jpsi_from_electrons_.jpsi_distance_xy        .push_back(zfe.reco_jpsi_from_electrons.distance_xy.at(i));
          reco_jpsi_from_electrons_.jpsi_distance_z         .push_back(zfe.reco_jpsi_from_electrons.distance_z.at(i));
          reco_jpsi_from_electrons_.jpsi_vtx_prob           .push_back(zfe.reco_jpsi_from_electrons.vtx_prob.at(i));
          reco_jpsi_from_electrons_.jpsi_vtx_x              .push_back(zfe.reco_jpsi_from_electrons.vtx_x.at(i));
          reco_jpsi_from_electrons_.jpsi_vtx_y              .push_back(zfe.reco_jpsi_from_electrons.vtx_y.at(i));
          reco_jpsi_from_electrons_.jpsi_vtx_z              .push_back(zfe.reco_jpsi_from_electrons.vtx_z.at(i));
          reco_jpsi_from_electrons_.jpsi_eff                .push_back(zfe.reco_jpsi_from_electrons.jpsi_efficiency.at(i));
          reco_jpsi_from_electrons_.jpsi_acc_eff            .push_back(zfe.reco_jpsi_from_electrons.jpsi_acc_eff.at(i));
          reco_jpsi_from_electrons_.jpsi_scale_factor       .push_back(zfe.reco_jpsi_from_electrons.jpsi_scale_factor.at(i));
          reco_jpsi_from_electrons_.muon0_pt                .push_back(zfe.reco_jpsi_from_electrons.muon0_pT.at(i));
          reco_jpsi_from_electrons_.muon0_eta               .push_back(zfe.reco_jpsi_from_electrons.muon0_eta.at(i));
          reco_jpsi_from_electrons_.muon0_phi               .push_back(zfe.reco_jpsi_from_electrons.muon0_phi.at(i));
          reco_jpsi_from_electrons_.muon0_d0                .push_back(zfe.reco_jpsi_from_electrons.muon0_d0.at(i));
          reco_jpsi_from_electrons_.muon0_dxy               .push_back(zfe.reco_jpsi_from_electrons.muon0_dxy.at(i));
          reco_jpsi_from_electrons_.muon0_dz                .push_back(zfe.reco_jpsi_from_electrons.muon0_dz.at(i));
          reco_jpsi_from_electrons_.muon0_d0err             .push_back(zfe.reco_jpsi_from_electrons.muon0_d0err.at(i));
          reco_jpsi_from_electrons_.muon0_dxyerr            .push_back(zfe.reco_jpsi_from_electrons.muon0_dxyerr.at(i));
          reco_jpsi_from_electrons_.muon0_dzerr             .push_back(zfe.reco_jpsi_from_electrons.muon0_dzerr.at(i));
          reco_jpsi_from_electrons_.muon1_pt                .push_back(zfe.reco_jpsi_from_electrons.muon1_pT.at(i));
          reco_jpsi_from_electrons_.muon1_eta               .push_back(zfe.reco_jpsi_from_electrons.muon1_eta.at(i));
          reco_jpsi_from_electrons_.muon1_phi               .push_back(zfe.reco_jpsi_from_electrons.muon1_phi.at(i));
          reco_jpsi_from_electrons_.muon1_d0                .push_back(zfe.reco_jpsi_from_electrons.muon1_d0.at(i));
          reco_jpsi_from_electrons_.muon1_dxy               .push_back(zfe.reco_jpsi_from_electrons.muon1_dxy.at(i));
          reco_jpsi_from_electrons_.muon1_dz                .push_back(zfe.reco_jpsi_from_electrons.muon1_dz.at(i));
          reco_jpsi_from_electrons_.muon1_d0err             .push_back(zfe.reco_jpsi_from_electrons.muon1_d0err.at(i));
          reco_jpsi_from_electrons_.muon1_dxyerr            .push_back(zfe.reco_jpsi_from_electrons.muon1_dxyerr.at(i));
          reco_jpsi_from_electrons_.muon1_dzerr             .push_back(zfe.reco_jpsi_from_electrons.muon1_dzerr.at(i));
          reco_jpsi_from_electrons_.muon0_charge            .push_back(zfe.reco_jpsi_from_electrons.muon0_charge.at(i));
          reco_jpsi_from_electrons_.muon1_charge            .push_back(zfe.reco_jpsi_from_electrons.muon1_charge.at(i));
          reco_jpsi_from_electrons_.has_muons_in_eta_window .push_back(zfe.reco_jpsi_from_electrons.has_muons_in_eta_window.at(i));
          reco_jpsi_from_electrons_.has_high_pt_muons       .push_back(zfe.reco_jpsi_from_electrons.has_high_pt_muons.at(i));
          reco_jpsi_from_electrons_.reco_jpsi_soft_id       .push_back(zfe.reco_jpsi_from_electrons.has_soft_id_muons.at(i));

          for (unsigned int iTrig = 0; iTrig < zfe.reco_z_from_muons.trigger_list.size(); ++iTrig)
            reco_z_from_muons_.passed_triggers.push_back(zfe.reco_z_from_muons.trigger_list.at(iTrig));

          event_.event_weight                                                .push_back(zfe.event_weight);
          event_.event_number                                                .push_back(zfe.id.event_num);
          event_.run_number                                                  .push_back(zfe.id.run_num);
          event_.n_verts                                                     .push_back(zfe.reco_vert.num);
          event_.found_good_electrons_from_z                                 .push_back(zfe.found_good_electrons_from_z);

        }
      }
    }

    // Truth
    if (!zfe.is_real_data) {
//      int n_jpsi = 0;
//      for (unsigned int i = 0; i < zfe.truth_jpsi.m.size() ; ++i ) {
//        n_jpsi++;
//        //TODO decide on how to do this for a tree (and histogram method too), for now just take first jpsi
//        //this is definitely a kludge, think of proper way to do this truth should only ever have 1 jpsi
//        truth_jpsi_.jpsi_m       .push_back(zfe.truth_jpsi.m.at(i));
//        truth_jpsi_.jpsi_pt      .push_back(zfe.truth_jpsi.pt.at(i));
//        truth_jpsi_.jpsi_y       .push_back(zfe.truth_jpsi.y.at(i));
//        truth_jpsi_.jpsi_phi     .push_back(zfe.truth_jpsi.phi.at(i));
//        truth_jpsi_.jpsi_eta     .push_back(zfe.truth_jpsi.eta.at(i));
//        truth_jpsi_.jpsi_vtx_x   .push_back(zfe.truth_jpsi.vtx_x.at(i));
//        truth_jpsi_.jpsi_vtx_y   .push_back(zfe.truth_jpsi.vtx_y.at(i));
//        truth_jpsi_.jpsi_vtx_z   .push_back(zfe.truth_jpsi.vtx_z.at(i));
//        truth_jpsi_.muon0_pt     .push_back(zfe.truth_jpsi.muon0_pT.at(i));  
//        truth_jpsi_.muon0_eta    .push_back(zfe.truth_jpsi.muon0_eta.at(i)); 
//        truth_jpsi_.muon0_phi    .push_back(zfe.truth_jpsi.muon0_phi.at(i)); 
//        truth_jpsi_.muon1_pt     .push_back(zfe.truth_jpsi.muon1_pT.at(i));  
//        truth_jpsi_.muon1_eta    .push_back(zfe.truth_jpsi.muon1_eta.at(i)); 
//        truth_jpsi_.muon1_phi    .push_back(zfe.truth_jpsi.muon1_phi.at(i)); 
//        truth_jpsi_.muon0_charge .push_back(zfe.truth_jpsi.muon0_pT.at(i));  
//        truth_jpsi_.muon1_charge .push_back(zfe.truth_jpsi.muon0_pT.at(i));  
//        truth_jpsi_.has_muons_in_eta_window.push_back(zfe.truth_jpsi.has_muons_in_eta_window.at(i));
//        truth_jpsi_.has_high_pt_muons      .push_back(zfe.truth_jpsi.has_high_pt_muons.at(i));
//      }
      for (unsigned int i = 0; i < zfe.truth_z_muons.muon0_pT.size() ; ++i ) {
        truth_z_muons_.z_m  .push_back(zfe.truth_z_muons.m);
        truth_z_muons_.z_y  .push_back(zfe.truth_z_muons.y);
        truth_z_muons_.z_pt .push_back(zfe.truth_z_muons.pt);
        truth_z_muons_.z_eta.push_back(zfe.truth_z_muons.eta);
        truth_z_muons_.z_phi.push_back(zfe.truth_z_muons.phi);
        truth_z_muons_.daughter0_pt.push_back(zfe.truth_z_muons.muon0_pT.at(i));
        truth_z_muons_.daughter1_pt.push_back(zfe.truth_z_muons.muon1_pT.at(i));
        truth_z_muons_.daughter0_eta.push_back(zfe.truth_z_muons.muon0_eta.at(i));
        truth_z_muons_.daughter1_eta.push_back(zfe.truth_z_muons.muon1_eta.at(i));
        truth_z_muons_.daughter0_phi.push_back(zfe.truth_z_muons.muon0_phi.at(i));
        truth_z_muons_.daughter1_phi.push_back(zfe.truth_z_muons.muon1_phi.at(i));
        truth_z_muons_.z_vtx_x.push_back(zfe.truth_z_muons.vtx_x.at(i));
        truth_z_muons_.z_vtx_y.push_back(zfe.truth_z_muons.vtx_y.at(i));
        truth_z_muons_.z_vtx_z.push_back(zfe.truth_z_muons.vtx_z.at(i));
      }
      for (unsigned int i = 0; i < zfe.truth_z_electrons.muon0_pT.size() ; ++i ) {
        truth_z_electrons_.z_m  .push_back(zfe.truth_z_electrons.m);
        truth_z_electrons_.z_y  .push_back(zfe.truth_z_electrons.y);
        truth_z_electrons_.z_pt .push_back(zfe.truth_z_electrons.pt);
        truth_z_electrons_.z_eta.push_back(zfe.truth_z_electrons.eta);
        truth_z_electrons_.z_phi.push_back(zfe.truth_z_electrons.phi);
        truth_z_electrons_.daughter0_pt.push_back (zfe.truth_z_electrons.muon0_pT.at(i));
        truth_z_electrons_.daughter1_pt.push_back (zfe.truth_z_electrons.muon1_pT.at(i));
        truth_z_electrons_.daughter0_eta.push_back(zfe.truth_z_electrons.muon0_eta.at(i));
        truth_z_electrons_.daughter1_eta.push_back(zfe.truth_z_electrons.muon1_eta.at(i));
        truth_z_electrons_.daughter0_phi.push_back(zfe.truth_z_electrons.muon0_phi.at(i));
        truth_z_electrons_.daughter1_phi.push_back(zfe.truth_z_electrons.muon1_phi.at(i));
        truth_z_electrons_.z_vtx_x.push_back      (zfe.truth_z_electrons.vtx_x.at(i));
        truth_z_electrons_.z_vtx_y.push_back      (zfe.truth_z_electrons.vtx_y.at(i));
        truth_z_electrons_.z_vtx_z.push_back      (zfe.truth_z_electrons.vtx_z.at(i));
      }
    }

    for (unsigned int i = 0; i < zfe.PDG_id.size() ; ++i ) {
      truth_z_muons_.truth_pdg_id.push_back(zfe.PDG_id.at(i));
    }
    // General Event info

    //------------
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
