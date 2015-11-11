#include <string.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TTree.h>
#include <TGraphAsymmErrors.h>
#include <TKey.h>
#include <TEfficiency.h>
#include <Riostream.h>
void calculate_jpsi_efficiencies (string file_name )
{
  
  //TODO maybe make the file name an input argument?
  //TFile *theFile0 = new TFile("/home/user1/turkewitz/Work/CMSSW_5_3_13_ZJPsi/src/test9d10.root");
  TFile *theFile0 = new TFile( file_name.c_str() );

  //TODO put this in rootrc
  gStyle->SetLineWidth(2.);
  gStyle->SetHistLineWidth(2.5);
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  //TH1D *h_n_truth_matched_muons_all = (TH1D*) theFile0->Get("ZFinder/All/N_truth_matched_jpsi_muons");
  //TH1D *h_n_truth_matched_muons_dimuon = (TH1D*) theFile0->Get("ZFinder/Dimuon_Jpsi/N_truth_matched_jpsi_muons");
  //TH1D *h_n_truth_matched_muons_dimuon_soft = (TH1D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Soft/N_truth_matched_jpsi_muons");
  //TH1D *h_n_truth_matched_muons_dimuon_vtx_comp = (TH1D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/N_truth_matched_jpsi_muons");
  //TH1D *h_n_truth_matched_muons_dimuon_primary_vert = (TH1D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Primary_Vertex/N_truth_matched_jpsi_muons");
  //TH1D *h_n_truth_matched_muons_jpsi = (TH1D*) theFile0->Get("ZFinder/Jpsi/N_truth_matched_jpsi_muons");

  //TH1D *h_jpsi_mass_fine_jpsi_mc = (TH1D*) theFile0->Get("ZFinder/MC_Jpsi/jpsi Mass: Fine");

  //double mc_events = h_jpsi_mass_fine_jpsi_mc->Integral();

  //double two_truth_matched_muons_dimuon = h_n_truth_matched_muons_dimuon->GetBinContent(3) ;
  //double two_truth_matched_muons_dimuon_soft = h_n_truth_matched_muons_dimuon_soft->GetBinContent(3);
  //double two_truth_matched_muons_dimuon_vtx_comp = h_n_truth_matched_muons_dimuon_vtx_comp->GetBinContent(3);
  //double two_truth_matched_muons_dimuon_primary_vert  = h_n_truth_matched_muons_dimuon_primary_vert->GetBinContent(3);
  //double two_truth_matched_muons_jpsi = h_n_truth_matched_muons_jpsi->GetBinContent(3);

  //std::cout << "mc_events " << mc_events << std::endl;
  //std::cout << "two_truth_matched_muons_dimuon " << two_truth_matched_muons_dimuon << std::endl;
  //std::cout << "two_truth_matched_muons_dimuon_soft " << two_truth_matched_muons_dimuon_soft  << std::endl;
  //std::cout << "two_truth_matched_muons_dimuon_vtx_comp " << two_truth_matched_muons_dimuon_vtx_comp  << std::endl;
  //std::cout << "two_truth_matched_muons_dimuon_primary_vert " << two_truth_matched_muons_dimuon_primary_vert  << std::endl;
  //std::cout << "two_truth_matched_muons_jpsi " << two_truth_matched_muons_jpsi  << std::endl;


  //TH2D *jpsi_pt_vs_rap_mc = (TH2D*) theFile0->Get("ZFinder/MC_All/jpsi_pt_vs_rap");
  //TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_pt_vs_rap");
  //TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get("ZFinder/Jpsi/jpsi_pt_vs_rap");
  //TH2D *jpsi_pt_vs_rap_mc = (TH2D*) theFile0->Get("ZFinder/MC_All/jpsi_pt_vs_rap_polarization_long");
  //TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get("ZFinder/Jpsi/jpsi_pt_vs_rap_polarization_long");
  //TH2D *jpsi_pt_vs_rap_mc = (TH2D*) theFile0->Get("ZFinder/MC_All/jpsi_pt_vs_rap_polarization_TPlusZero");
  //TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get("ZFinder/Jpsi/jpsi_pt_vs_rap_polarization_TPlusZero");


  TH2D *jpsi_pt_vs_rap_mc = (TH2D*) theFile0->Get("ZFinder/MC_All/jpsi_pt_vs_rap");
  TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_pt_vs_rap");
  //TH2D *jpsi_pt_vs_rap_mc = (TH2D*) theFile0->Get("ZFinder/MC_All/jpsi_pt_vs_rap_polarization_long");
  //TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_pt_vs_rap_polarization_long");
  //TH2D *jpsi_pt_vs_rap_mc = (TH2D*) theFile0->Get("ZFinder/MC_All/jpsi_pt_vs_rap_polarization_TPlusZero");
  //TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_pt_vs_rap_polarization_TPlusZero");

  //TODO uncomment following two lines as we only use 1 bin of pT
  //jpsi_pt_vs_rap_mc->Rebin2D(21,1);
  //jpsi_pt_vs_rap_jpsi->Rebin2D(21,1);
  std::cout << jpsi_pt_vs_rap_mc->Integral() << std::endl;
  std::cout << jpsi_pt_vs_rap_jpsi->Integral() << std::endl;
  jpsi_pt_vs_rap_mc->Sumw2();
  jpsi_pt_vs_rap_jpsi->Sumw2();
  TH2D *acc_eff_map = jpsi_pt_vs_rap_mc->Clone();
  acc_eff_map->Divide(jpsi_pt_vs_rap_jpsi, jpsi_pt_vs_rap_mc, 1.0, 1.0, "B");
  //acc_eff_map->Divide(jpsi_pt_vs_rap_jpsi, jpsi_pt_vs_rap_mc, 1.0, 1.0);
  //acc_eff_map->Divide(jpsi_pt_vs_rap_jpsi, jpsi_pt_vs_rap_mc, "cl=0.683 b(1,1) mode");



  acc_eff_map->Draw("colz");
  

  //TEfficiency * pEff = 0;
  //if (TEfficiency::CheckConsistency(*jpsi_pt_vs_rap_jpsi, *jpsi_pt_vs_rap_mc, "w")) { 
  //  pEff = new TEfficiency(*jpsi_pt_vs_rap_jpsi, *jpsi_pt_vs_rap_mc);
  //}
  //pEff->Draw("colz");

  int nxbins = acc_eff_map->GetNbinsX();
  int nybins = acc_eff_map->GetNbinsY();
  for (int x=1;x<=nxbins;++x)
  {
    for (int y=1;y<=nybins;++y)
    {
      double acc_eff = acc_eff_map->GetBinContent(x,y);
      double acc_error_low = acc_eff_map->GetBinError(x,y);
      double acc_error_high = acc_eff_map->GetBinError(x,y);
      //double acc_error_low = acc_eff_map->GetBinErrorLow(x,y);
      //double acc_error_low = acc_eff_map->GetBinErrorLow(x,y);
      //double acc_error_high = acc_eff_map->GetBinErrorHigh(x,y);
      TAxis *xaxis = acc_eff_map->GetXaxis();  
      double xlow = xaxis->GetBinLowEdge(x);
      double xwidth = xaxis->GetBinWidth(x);
      double xhigh = xlow + xwidth;
      TAxis *yaxis = acc_eff_map->GetYaxis();  
      double ylow = yaxis->GetBinLowEdge(y);
      double ywidth = yaxis->GetBinWidth(y);
      double yhigh = ylow + ywidth;
      //std::cout << "xlow " << xlow  << " xhigh " << xhigh << " ylow " << ylow << " yhigh " << yhigh << " acc_eff " << acc_eff <<
      //  " acc_error_low " << acc_error_low << " acc_error_high " << acc_error_high << std::endl;
      std::cout << "{" << xlow << ",  " << xhigh << ",  " << ylow << ",  " << yhigh << ",  " <<
        acc_eff << ",  " << acc_error_low << ",  " << acc_error_high << "   }," << std::endl;
    }
  }

  //TODO verify just dividing seems to give similar/same results instead of using TEfficiency
  //pEff->SetStatisticOption(1);
  //int nxbins = jpsi_pt_vs_rap_mc->GetNbinsX();
  //int nybins = jpsi_pt_vs_rap_mc->GetNbinsY();
  //for (int x=1;x<=nxbins;++x)
  //{
  //  for (int y=1;y<=nybins;++y)
  //  {
  //    int bin = pEff->GetGlobalBin(x,y);
  //    double acc_eff = pEff->GetEfficiency(bin);
  //    double acc_error_low = pEff->GetEfficiencyErrorLow(bin);
  //    double acc_error_high = pEff->GetEfficiencyErrorUp(bin);
  //    TAxis *xaxis = jpsi_pt_vs_rap_mc->GetXaxis();  
  //    double xlow = xaxis->GetBinLowEdge(x);
  //    double xwidth = xaxis->GetBinWidth(x);
  //    double xhigh = xlow + xwidth;
  //    TAxis *yaxis = jpsi_pt_vs_rap_mc->GetYaxis();  
  //    double ylow = yaxis->GetBinLowEdge(y);
  //    double ywidth = yaxis->GetBinWidth(y);
  //    double yhigh = ylow + ywidth;
  //    //std::cout << "xlow " << xlow  << " xhigh " << xhigh << " ylow " << ylow << " yhigh " << yhigh << " acc_eff " << acc_eff <<
  //    //  " acc_error_low " << acc_error_low << " acc_error_high " << acc_error_high << std::endl;
  //    //eta low eta high pt_low pt_high acc_eff acc_eff_err_low acc_eff_err_high
  //    std::cout << "{" << xlow << ",  " << xhigh << ",  " << ylow << ",  " << yhigh << ",  " <<
  //      acc_eff << ",  " << acc_error_low << ",  " << acc_error_high << "   }," << std::endl;
  //  }
  //}

  TFile output("acc_eff_map.root","new");
  acc_eff_map->Write();
}
