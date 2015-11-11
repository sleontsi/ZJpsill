#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
void make_images2 (string file_name_mc, string file_name_data_zmumu, string file_name_data_zee)
{
  
  //TODO give this better documentation/ a more descriptive name
  //TFile *file_mc = new TFile("/home/user1/turkewitz/Work/CMSSW_5_3_13_ZJPsi/src/test9d10.root");
  TFile *file_mc = new TFile( file_name_mc.c_str() );
  TFile *file_data_zmumu = new TFile( file_name_data_zmumu.c_str() );
  TFile *file_data_zee = new TFile( file_name_data_zee.c_str() );
  std::string image_name1 = "";
  std::string image_name2 = "";
  std::string image_name3 = "";
  std::string image_name4 = "";
  std::string image_name5 = "";
  std::string image_name6 = "";
  std::string image_name7 = "";
  std::string image_name8 = "";
  std::string image_name9 = "";
  std::string image_name10 = "";
  std::string image_name11 = "";
  std::string image_name12 = "";
  std::string image_name13 = "";
  std::string image_name14 = "";
  std::string image_name15 = "";
  std::string image_name16 = "";
  path = "/home/user1/turkewitz/public_html/ZPhysics/tmp/Test50/";
  image_name1.append(path);
  image_name2.append(path);
  image_name3.append(path);
  image_name4.append(path);
  image_name5.append(path);
  image_name6.append(path);
  image_name7.append(path);
  image_name8.append(path);
  image_name9.append(path);
  image_name10.append(path);
  image_name11.append(path);
  image_name12.append(path);
  image_name13.append(path);
  image_name14.append(path);
  image_name15.append(path);
  image_name16.append(path);

  //TODO put this in rootrc
  gStyle->SetLineWidth(2.);
  gStyle->SetHistLineWidth(2.5);
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  //Jpsi
  TH1D *h_truth_vtx_minus_reco_vtx = (TH1D*) file_mc->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_truth_vtx_z_minus_jpsi_reco_vtx_z_");

  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->cd();
  c1->SetLogy();
  h_truth_vtx_minus_reco_vtx->GetXaxis()->SetRangeUser(-1,1);
  h_truth_vtx_minus_reco_vtx->SetTitle("J/#psi MC");
  h_truth_vtx_minus_reco_vtx->Draw();
  image_name1 = image_name1.append("jpsi_truth_vtx_z_minus_jpsi_reco_vtx_z");
  image_name1.append(".png");
  c1->Print(image_name1.c_str() , "png");
  //c1->Close();
  //delete c1;

  //TODO testing
  TFile *file_mc_zjpsi_mpi = new TFile("zjpsi_mpi_hadded.root"  );
  TH1D *h_jpsi_pt_zjpsi_mpi = (TH1D*) file_mc_zjpsi_mpi->Get("ZFinder/Z_To_Muons_And_Jpsi/jpsi p_{T}");
  TH1D *h_jpsi_pt = (TH1D*) file_mc->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi p_{T}");
  TH1D *h_jpsi_pt_zmumu_jpsi = (TH1D*) file_data_zmumu->Get("ZFinder/Z_To_Muons_And_Jpsi/jpsi p_{T}");
  TH1D *h_jpsi_pt_zee_jpsi = (TH1D*) file_data_zee->Get("ZFinder/Z_To_Electrons_And_Jpsi/jpsi p_{T}");
  h_jpsi_pt_zmumu_jpsi->Add(h_jpsi_pt_zee_jpsi);
  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->cd();
  c2->SetLogy();
  h_jpsi_pt_zjpsi_mpi->Rebin(5);
  h_jpsi_pt->Rebin(5);
  h_jpsi_pt_zmumu_jpsi->Rebin(5);
  h_jpsi_pt->GetXaxis()->SetRangeUser(8.5,100);
  h_jpsi_pt->GetXaxis()->SetTitle("J/#psi p_{T}");
  h_jpsi_pt->GetYaxis()->SetTitle("Norm. Events / 5 GeV");
  h_jpsi_pt->SetTitle("J/#psi p_{T}");
  double norm_mc_zjpsi = h_jpsi_pt_zjpsi_mpi->GetEntries();
  double norm_mc = h_jpsi_pt->GetEntries();
  double norm_data = h_jpsi_pt_zmumu_jpsi->GetEntries();
  //TODO testing
  h_jpsi_pt->Scale(1/norm_mc*norm_data);
  //h_jpsi_pt->GetYaxis()->SetRangeUser(1e-6,0.9);
  h_jpsi_pt->Draw();
  //h_jpsi_pt_zmumu_jpsi->Scale(1/norm_data);
  h_jpsi_pt_zmumu_jpsi->SetLineColor(kRed-2);
  h_jpsi_pt_zmumu_jpsi->Draw("same");


  h_jpsi_pt_zjpsi_mpi->Scale(1/norm_mc_zjpsi*norm_data);
  h_jpsi_pt_zjpsi_mpi->SetLineColor(kGreen-2);
  h_jpsi_pt_zjpsi_mpi->Draw("same");

  Double_t xl1=.55, yl1=0.70, xl2=xl1+.30, yl2=yl1+.20;
  TLegend *leg2 = new TLegend(xl1,yl1,xl2,yl2);
  leg2->AddEntry(h_jpsi_pt_zmumu_jpsi,"Z->ll + J/#psi Data","lep");
  leg2->AddEntry(h_jpsi_pt,"J/#psi MC","lep");
  leg2->AddEntry(h_jpsi_pt_zjpsi_mpi,"J/#psi MC (Z+J/#psi MPI)","lep");
  leg2->Draw();

  image_name2 = image_name2.append("jpsi_pt");
  image_name2.append(".png");
  c2->Print(image_name2.c_str() , "png");

  //TODO fix this
  TH1D *h_jpsi_rap_zjpsi_mpi = (TH1D*) file_mc_zjpsi_mpi->Get("ZFinder/Z_To_Muons_And_Jpsi/jpsi Rapidity");
  TH1D *h_jpsi_rap = (TH1D*) file_mc->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi Rapidity");
  TH1D *h_jpsi_rap_zmumu = (TH1D*) file_data_zmumu->Get("ZFinder/Z_To_Muons_And_Jpsi/jpsi Rapidity");
  TH1D *h_jpsi_rap_zee = (TH1D*) file_data_zee->Get("ZFinder/Z_To_Electrons_And_Jpsi/jpsi Rapidity");
  h_jpsi_rap_zmumu->Add(h_jpsi_rap_zee);
  TCanvas *c3 = new TCanvas("c3", "c3");
  c3->cd();
  c3->SetLogy();
  h_jpsi_rap->SetTitle("J/#psi Rapidity");
  h_jpsi_rap->GetXaxis()->SetTitle("J/#psi Rap");
  h_jpsi_rap->GetXaxis()->SetRangeUser(-2.1,2.1);
  h_jpsi_rap->GetYaxis()->SetTitle("Norm. Events / 0.1");

  h_jpsi_rap_zjpsi_mpi->Rebin(5);
  h_jpsi_rap->Rebin(5);
  h_jpsi_rap_zmumu->Rebin(5);

  double norm_mc_zjpsi = h_jpsi_rap_zjpsi_mpi->GetEntries();
  double norm_mc = h_jpsi_rap->GetEntries();
  double norm_data = h_jpsi_rap_zmumu->GetEntries();
  h_jpsi_rap->Scale(norm_data/norm_mc);
  //h_jpsi_rap->GetYaxis()->SetRangeUser(1e-4,0.9);
  h_jpsi_rap->Draw();
  h_jpsi_rap_zjpsi_mpi->SetLineColor(kGreen-2);
  h_jpsi_rap_zjpsi_mpi->Scale(norm_data/norm_mc_zjpsi);
  h_jpsi_rap_zjpsi_mpi->Draw("same");
  //h_jpsi_rap_zmumu->Scale(1/);
  h_jpsi_rap_zmumu->SetLineColor(kRed-2);
  h_jpsi_rap_zmumu->Draw("same");

  Double_t xl1=.55, yl1=0.70, xl2=xl1+.30, yl2=yl1+.20;
  TLegend *leg3 = new TLegend(xl1,yl1,xl2,yl2);
  leg3->AddEntry(h_jpsi_rap_zmumu,"Z->ll + J/#psi Data","lep");
  leg3->AddEntry(h_jpsi_rap,"J/#psi MC","lep");
  leg3->AddEntry(h_jpsi_pt_zjpsi_mpi,"J/#psi MC (Z+J/#psi MPI)","lep");
  leg3->Draw();
  image_name3.append("jpsi_rap");
  image_name3.append(".png");
  c3->Print(image_name3.c_str() , "png");

  TH1D *h_jpsi_mass = (TH1D*) file_mc->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_mass");
  TCanvas *c4 = new TCanvas("c4", "c4");
  c4->cd();
  c4->SetLogy();
  //h_jpsi_pt->GetXaxis()->SetRangeUser(-1,1);
  h_jpsi_mass->SetTitle("J/#psi MC");
  h_jpsi_mass->Draw();
  image_name4.append("jpsi_mass");
  image_name4.append(".png");
  c4->Print(image_name4.c_str() , "png");

  TH1D *h_jpsi_cos_mu_plus = (TH1D*) file_mc->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_cos_mu_plus");
  TCanvas *c5 = new TCanvas("c5", "c5");
  c5->cd();
  //c5->SetLogy();
  //h_jpsi_pt->GetXaxis()->SetRangeUser(-1,1);
  h_jpsi_cos_mu_plus->SetTitle("J/#psi MC");
  h_jpsi_cos_mu_plus->Draw();
  image_name5.append("jpsi_cos_mu_plus");
  image_name5.append(".png");
  c5->Print(image_name5.c_str() , "png");

  TH1D *h_jpsi_tau_xy = (TH1D*) file_mc->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_tau_xy_very_fine_all");
  TCanvas *c6 = new TCanvas("c6", "c6");
  c6->cd();
  c6->SetLogy();
  //h_jpsi_pt->GetXaxis()->SetRangeUser(-1,1);
  h_jpsi_tau_xy->SetTitle("J/#psi MC");
  h_jpsi_tau_xy->Draw();
  image_name6.append("jpsi_tau_xy_very_fine_al");
  image_name6.append(".png");
  c6->Print(image_name6.c_str() , "png");

  TH2D *h_dimuon_mass_vs_tau_xy = (TH2D*) file_mc->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/dimuon_mass_vs_dimuon_tau_xy_fine");
  TCanvas *c7 = new TCanvas("c7", "c7");
  c7->cd();
  //h_jpsi_pt->GetXaxis()->SetRangeUser(-1,1);
  c7->SetLogz();
  h_dimuon_mass_vs_tau_xy->SetTitle("J/#psi MC");
  h_dimuon_mass_vs_tau_xy->Draw("colz");
  image_name7.append("dimuon_mass_vs_dimuon_tau_xy");
  image_name7.append(".png");
  c7->Print(image_name7.c_str() , "png");

  //TH1D *h_mu0_pt = (TH1D*) file_mc->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/p_{T,mu_{0}}");
  //TCanvas *c8 = new TCanvas("c8", "c8");
  //c8->cd();
  //c8->SetLogy();
  //h_mu0_pt->SetTitle("J/#psi MC");
  //h_mu0_pt->Draw("");
  //image_name8.append("mu0_pT");
  //image_name8.append(".png");
  //c8->Print(image_name8.c_str() , "png");

  TH1D *h_mu0_pt = (TH1D*) file_mc->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/p_{T,mu_{0}}");
  TH1D *h_mu0_pt_zmumu = (TH1D*) file_data_zmumu->Get("ZFinder/Z_To_Muons_And_Jpsi/p_{T,mu_{0}}");
  TH1D *h_mu0_pt_zee = (TH1D*) file_data_zee->Get("ZFinder/Z_To_Electrons_And_Jpsi/p_{T,mu_{0}}");
  h_mu0_pt_zmumu->Add(h_mu0_pt_zee);
  TCanvas *c8 = new TCanvas("c8", "c8");
  c8->cd();
  c8->SetLogy();
  h_mu0_pt->SetTitle("J/#psi #mu_{0} p_{T}");
  h_mu0_pt->GetXaxis()->SetTitle("#mu_{0} p_{T}");
  h_mu0_pt->GetXaxis()->SetRangeUser(3.5,100);
  h_mu0_pt->GetYaxis()->SetTitle("Norm. Events / 5 GeV");
  h_mu0_pt->Rebin(10);
  h_mu0_pt_zmumu->Rebin(10);
  h_mu0_pt->Draw();

  double norm_mc = h_jpsi_rap->GetEntries();
  h_mu0_pt->Scale(1/norm_mc);
  h_mu0_pt->GetYaxis()->SetRangeUser(1e-5,0.9);
  h_mu0_pt->Draw();
  double norm_data = h_mu0_pt_zmumu->GetEntries();
  h_mu0_pt_zmumu->Scale(1/norm_data);
  h_mu0_pt_zmumu->SetLineColor(kRed-2);
  h_mu0_pt_zmumu->Draw("same");

  Double_t xl1=.55, yl1=0.70, xl2=xl1+.30, yl2=yl1+.20;
  TLegend *leg8 = new TLegend(xl1,yl1,xl2,yl2);
  leg8->AddEntry(h_jpsi_rap_zmumu,"Z->ll + J/#psi Data","lep");
  leg8->AddEntry(h_jpsi_rap,"J/#psi MC","lep");
  leg8->Draw();
  image_name8.append("mu0_pT");
  image_name8.append(".png");
  c8->Print(image_name8.c_str() , "png");

  //TH1D *h_mu1_pt = (TH1D*) file_mc->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/p_{T,mu_{1}}");
  //TCanvas *c9 = new TCanvas("c9", "c9");
  //c9->cd();
  ////h_jpsi_pt->GetXaxis()->SetRangeUser(-1,1);
  //c9->SetLogy();
  //h_mu1_pt->SetTitle("J/#psi MC");
  //h_mu1_pt->Draw("");
  //image_name9.append("mu1_pT");
  //image_name9.append(".png");
  //c9->Print(image_name9.c_str() , "png");

  TH1D *h_mu1_pt = (TH1D*) file_mc->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/p_{T,mu_{1}}");
  TH1D *h_mu1_pt_zmumu = (TH1D*) file_data_zmumu->Get("ZFinder/Z_To_Muons_And_Jpsi/p_{T,mu_{1}}");
  TH1D *h_mu1_pt_zee = (TH1D*) file_data_zee->Get("ZFinder/Z_To_Electrons_And_Jpsi/p_{T,mu_{1}}");
  h_mu1_pt_zmumu->Add(h_mu1_pt_zee);
  TCanvas *c9 = new TCanvas("c9", "c9");
  c9->cd();
  c9->SetLogy();
  h_mu1_pt->SetTitle("J/#psi #mu_{1} p_{T}");
  h_mu1_pt->GetXaxis()->SetTitle("#mu_{1} p_{T}");
  h_mu1_pt->GetXaxis()->SetRangeUser(3.5,100);
  h_mu1_pt->GetYaxis()->SetTitle("Norm. Events / GeV");
  h_mu1_pt->Rebin(2);
  h_mu1_pt_zmumu->Rebin(2);
  h_mu1_pt->Draw();

  double norm_mc = h_jpsi_rap->GetEntries();
  h_mu1_pt->Scale(1/norm_mc);
  h_mu1_pt->GetYaxis()->SetRangeUser(1e-5,0.9);
  h_mu1_pt->Draw();
  double norm_data = h_mu1_pt_zmumu->GetEntries();
  h_mu1_pt_zmumu->Scale(1/norm_data);
  h_mu1_pt_zmumu->SetLineColor(kRed-2);
  h_mu1_pt_zmumu->Draw("same");

  Double_t xl1=.55, yl1=0.70, xl2=xl1+.30, yl2=yl1+.20;
  TLegend *leg9 = new TLegend(xl1,yl1,xl2,yl2);
  leg9->AddEntry(h_jpsi_rap_zmumu,"Z->ll + J/#psi Data","lep");
  leg9->AddEntry(h_jpsi_rap,"J/#psi MC","lep");
  leg9->Draw();
  image_name9.append("mu1_pT");
  image_name9.append(".png");
  c9->Print(image_name9.c_str() , "png");

  TH1D *h_mu0_eta = (TH1D*) file_mc->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/#eta_{mu_{0}}");
  TCanvas *c10 = new TCanvas("c10", "c10");
  c10->cd();
  //h_jpsi_pt->GetXaxis()->SetRangeUser(-1,1);
  c10->SetLogy();
  h_mu0_eta->SetTitle("J/#psi MC");
  h_mu0_eta->Draw("");
  image_name10.append("mu0_eta");
  image_name10.append(".png");
  c10->Print(image_name10.c_str() , "png");

  TH1D *h_mu1_eta = (TH1D*) file_mc->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/#eta_{mu_{1}}");
  TCanvas *c11 = new TCanvas("c11", "c11");
  c11->cd();
  //h_jpsi_pt->GetXaxis()->SetRangeUser(-1,1);
  c11->SetLogy();
  h_mu1_eta->SetTitle("J/#psi MC");
  h_mu1_eta->Draw("");
  image_name11.append("mu1_eta");
  image_name11.append(".png");
  c11->Print(image_name11.c_str() , "png");
}
