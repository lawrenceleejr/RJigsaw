#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>
#include <TF1.h>
#include <iostream>
#include <string>
#include <vector>
#include <TLegend.h>

#include "../RJigsaw/TRJigsaw.h"

using namespace std; 

void setstyle();
TColor *icolor[9][2];
int color_list[10];

///////////stuff for getting gamma distribution from PDF's
int type; //0-q qbar 1-gg

double rootS, M, C;

double eta_1 = 0.27871;
double eta_2 = 3.3627;
double eps_u = 4.4343;
double g_u = 38.599;

double del_S = -0.11912;
double eta_S = 9.4189;
double eps_S = -2.6287;
double g_S = 18.065;

double A_g = 3.4055;
double del_g = -0.12178;
double eta_g = 2.9278;
double eps_g = -2.3210;
double g_g = 1.9233;
double A_g1 = -1.6189;
double del_g1 = -0.23999;
double eta_g1 = 24.792;

double PDF_v(double x);
double PDF_S(double x);
double PDF_g(double x);
double PDF_TOT(double x);
double Calc_dsigma_dgamma(double gamma);

///////// particle 4-vectors
TLorentzVector P[2]; //parent particles 
TLorentzVector C_1[2], C_2[2]; //first set of child particles, indexed by hemisphere
TLorentzVector C_1_1[2], C_1_2[2]; //second set of child particles, indexed by hemisphere
TLorentzVector C_1_1_1[2], C_1_1_2[2]; //third set of child particles, indexed by hemisphere
TLorentzVector C_1_1_1_1[2], C_1_1_1_2[2]; //fourth set of child particles, indexed by hemisphere
TLorentzVector vMiss[2];

void BoostToLabFrame(double gamma_eff);
void BoostHem(int iHem, TVector3 vBETA);
void PtBoost(double PT);

void GetLepTopHem(int iHem, double Mtop, double MW, double Mb, double Mlep, double Mnu);

TCanvas* Plot_Me_2D(char *titlecan, TH2D* histo, char *titleX, char *titleY);

void RJ_ttbar(){
  setstyle();
	
  //give transverse momenta to CM system in lab frame?
  double PT = 0.1; //In units of sqrt{shat}
	
  //gamma factor associated with 'off-threshold-ness' of tops
  double gamma = 1.2;
	
  //Now, we also have the option to take gamma, event-by-event, 
  //from a more realistic distribution
  //to do this, set 'b_gamma' to true and the rest
  //of the variables below appropriately
  bool b_gamma = true;
  rootS = 13.;
  type = 1; //0 quark-antiquark  1 gluon-gluon
  M = 175./1000.; //TeV units
	
  //Number of toy events to throw
  int N = 100;

  Root::TRJigsaw* RJTool = new Root::TRJigsaw();
	
  //setup histograms
  RJTool->resetHists();

  RJTool->bookHist(1,"dphiTT_-1_vs_dphiTT_0",50, 0.0, TMath::Pi(), 50, 0.0,TMath::Pi());
  RJTool->bookHist(1,"dphiTT_-1_vs_dphiTT_1",50, 0.0, TMath::Pi(), 50, 0.0,TMath::Pi());
  RJTool->bookHist(1,"dphiM_-1_vs_dphiM_0",50, 0.0, TMath::Pi(), 50, 0.0,TMath::Pi());
  RJTool->bookHist(1,"dphiM_-1_vs_dphiM_1",50, 0.0, TMath::Pi(), 50, 0.0,TMath::Pi());
  RJTool->bookHist(1,"dphiTT_-1_vs_dphiM_-1",50, 0.0, TMath::Pi(), 50, 0.0,TMath::Pi());
  RJTool->bookHist(1,"dphiTT_0_vs_dphiM_0",50, 0.0, TMath::Pi(), 50, 0.0,TMath::Pi());
  RJTool->bookHist(1,"dphiTT_1_vs_dphiM_1",50, 0.0, TMath::Pi(), 50, 0.0,TMath::Pi());

  RJTool->bookHist(1,"MTT_0_vs_dphiM_0",50, 0.0, TMath::Pi(), 50, 0.0,TMath::Pi());
  RJTool->bookHist(1,"MTT_1_vs_dphiM_1",50, 0.0, TMath::Pi(), 50, 0.0,TMath::Pi());

  RJTool->bookHist(1,"costhetaT1_-1_vs_costhetaT2_-1",50, -1., 1., 50, -1.,1.);
  RJTool->bookHist(1,"costhetaT1_0_vs_costhetaT2_0",50, -1., 1., 50, -1.,1.);
  RJTool->bookHist(1,"costhetaT1_1_vs_costhetaT2_1",50, -1., 1., 50, -1.,1.);


  TH2D *hist_dcosthetaT1_costhetaT1_top = new TH2D("hist_dcosthetaT1_costhetaT1_top","hist_dcosthetaT1_costhetaT1_top",50, -1., 1., 50, -1.0,1.);
  TH2D *hist_dcosthetaT1_costhetaT1_W = new TH2D("hist_dcosthetaT1_costhetaT1_W","hist_dcosthetaT1_costhetaT1_W",50, -1., 1., 50, -1.0,1.);
  TH2D *hist_dcosthetaT1_dcosthetaT2_top = new TH2D("hist_dcosthetaT1_dcosthetaT1_top","hist_dcosthetaT1_dcosthetaT1_top",50, -1., 1., 50, -1.0,1.);
  TH2D *hist_dcosthetaT1_dcosthetaT2_W = new TH2D("hist_dcosthetaT1_dcosthetaT1_W","hist_dcosthetaT1_dcosthetaT1_W",50, -1., 1., 50, -1.0,1.);
  TH2D *hist_gammaTrue_gamma_top = new TH2D("hist_gammaTrue_gamma_top","hist_gammaTrue_gamma_top",50,1.,6,50,1.,6.);
  TH2D *hist_gammaTrue_gamma_W = new TH2D("hist_gammaTrue_gamma_W","hist_gammaTrue_gamma_W",50,1.,6,50,1.,6.);
  TH2D *hist_dcosthetaT1_gammaT_top = new TH2D("hist_dcosthetaT1_gammaT_top","hist_dcosthetaT1_gammaT_top",50, -1., 1., 50, 1.,6.);
  TH2D *hist_dcosthetaT1_gammaT_W = new TH2D("hist_dcosthetaT1_gammaT_W","hist_dcosthetaT1_gammaT_W",50, -1., 1., 50, 1.,6.);
  TH2D *hist_Eb1_gammaT_top = new TH2D("hist_Eb1_gammaT_top","hist_Eb1_gammaT_top",50, 0., 3., 50, 1.,6.);
  TH2D *hist_Eb1_gammaT_W = new TH2D("hist_Eb1_gammaT_W","hist_Eb1_gammaT_W",50, 0., 3., 50, 1.,6.);
  TH2D *hist_Eb1_costhetaT_top = new TH2D("hist_Eb1_costhetaT_top","hist_Eb1_costhetaT_top",50, 0., 3., 50, -1.,1.);
  TH2D *hist_Eb1_costhetaT_W = new TH2D("hist_Eb1_costhetaT_W","hist_Eb1_costhetaT_W",50, 0., 3., 50, -1.,1.);
  TH2D *hist_Eb1_dcosthetaT_top = new TH2D("hist_Eb1_dcosthetaT_top","hist_Eb1_dcosthetaT_top",50, 0., 3., 50, -1.,1.);
  TH2D *hist_Eb1_dcosthetaT_W = new TH2D("hist_Eb1_dcosthetaT_W","hist_Eb1_dcosthetaT_W",50, 0., 3., 50, -1.,1.);
  TH2D *hist_Eb1_Eb2_top = new TH2D("hist_Eb1_Eb2_top","hist_Eb1_Eb2_top",50, 0., 3., 50, 0.,3.);
  TH2D *hist_Eb1_Eb2_W = new TH2D("hist_Eb1_Eb2_W","hist_Eb1_Eb2_W",50, 0., 3., 50, 0.,3.);
  TH2D *hist_MTT_Eb1_top = new TH2D("hist_MTT_Eb1_top","hist_MTT_Eb1_top",50, 0., 3., 50, 0.,3.);
  TH2D *hist_MTT_Eb1_W = new TH2D("hist_MTT_Eb1_W","hist_MTT_Eb1_W",50, 0., 3., 50, 0.,3.);
  TH2D *hist_MT1_MT2_W = new TH2D("hist_MT1_MT2_W","hist_MT1_MT2_W",50, 0., 3., 50, 0.,3.);
  TH2D *hist_costhetaW1_costhetaW2 = new TH2D("hist_costhetaW1_costhetaW2","hist_costhetaW1_costhetaW2",50, -1., 1., 50, -1.0,1.);
  TH2D *hist_costhetaW1_costhetaW2_top = new TH2D("hist_costhetaW1_costhetaW2_top","hist_costhetaW1_costhetaW2_top",50, -1., 1., 50, -1.0,1.);
  TH2D *hist_costhetaW1_costhetaW2_W = new TH2D("hist_costhetaW1_costhetaW2_W","hist_costhetaW1_costhetaW2_W",50, -1., 1., 50, -1.0,1.);
  TH2D *hist_dcosthetaW1_costhetaW1_top = new TH2D("hist_dcosthetaW1_costhetaW1_top","hist_dcosthetaW1_costhetaW1_top",50, -1., 1., 50, -1.0,1.);
  TH2D *hist_dcosthetaW1_costhetaW1_W = new TH2D("hist_dcosthetaW1_costhetaW1_W","hist_dcosthetaW1_costhetaW1_W",50, -1., 1., 50, -1.0,1.);
  TH2D *hist_dcosthetaW1_dcosthetaW2_top = new TH2D("hist_dcosthetaW1_dcosthetaW1_top","hist_dcosthetaW1_dcosthetaW1_top",50, -1., 1., 50, -1.0,1.);
  TH2D *hist_dcosthetaW1_dcosthetaW2_W = new TH2D("hist_dcosthetaW1_dcosthetaW1_W","hist_dcosthetaW1_dcosthetaW1_W",50, -1., 1., 50, -1.0,1.);
  TH2D *hist_El1_gammaT_top = new TH2D("hist_El1_gammaT_top","hist_El1_gammaT_top",50, 0., 3., 50, 1.,6.);
  TH2D *hist_El1_gammaT_W = new TH2D("hist_El1_gammaT_W","hist_El1_gammaT_W",50, 0., 3., 50, 1.,6.);
  TH2D *hist_El1_costhetaW_top = new TH2D("hist_El1_costhetaW_top","hist_El1_costhetaW_top",50, 0., 3., 50, -1.,1.);
  TH2D *hist_El1_costhetaW_W = new TH2D("hist_El1_costhetaW_W","hist_El1_costhetaW_W",50, 0., 3., 50, -1.,1.);
  TH2D *hist_El1_dcosthetaW_top = new TH2D("hist_El1_dcosthetaW_top","hist_El1_dcosthetaW_top",50, 0., 3., 50, -1.,1.);
  TH2D *hist_El1_dcosthetaW_W = new TH2D("hist_El1_dcosthetaT_W","hist_El1_dcosthetaW_W",50, 0., 3., 50, -1.,1.);
  TH2D *hist_El1_El2_top = new TH2D("hist_El1_El2_top","hist_El1_El2_top",50, 0., 3., 50, 0.,3.);
  TH2D *hist_El1_El2_W = new TH2D("hist_El1_El2_W","hist_El1_El2_W",50, 0., 3., 50, 0.,3.);
  TH2D *hist_MTT_El1_top = new TH2D("hist_MTT_El1_top","hist_MTT_El1_top",50, 0., 3., 50, 0.,3.);
  TH2D *hist_MTT_El1_W = new TH2D("hist_MTT_El1_W","hist_MTT_El1_W",50, 0., 3., 50, 0.,3.);
  TH2D *hist_Eb1_El1_top = new TH2D("hist_Eb1_El1_top","hist_Eb1_El1_top",50, 0., 3., 50, 0.,3.);
  TH2D *hist_Eb1_El1_W = new TH2D("hist_Eb1_El1_W","hist_Eb1_El1_W",50, 0., 3., 50, 0.,3.);
  TH2D *hist_Eb1_El2_top = new TH2D("hist_Eb1_El2_top","hist_Eb1_El2_top",50, 0., 3., 50, 0.,3.);
  TH2D *hist_Eb1_El2_W = new TH2D("hist_Eb1_El2_W","hist_Eb1_El2_W",50, 0., 3., 50, 0.,3.);
  TH2D *hist_El1_Mt1_top = new TH2D("hist_El1_Mt1_top","hist_El1_Mt1_top",50, 0., 3., 50, 0.,3.);
  TH2D *hist_El1_Mt1_W = new TH2D("hist_El1_Mt1_W","hist_El1_Mt1_W",50, 0., 3., 50, 0.,3.);
  TH2D *hist_El1_costhetaT1_top = new TH2D("hist_El1_costhetaT1_top","hist_El1_costhetaT1_top",50, 0., 3., 50, -1.,1.);
  TH2D *hist_El1_costhetaT1_W = new TH2D("hist_El1_costhetaT1_W","hist_El1_costhetaT1_W",50, 0., 3., 50, -1.,1.);

  TH2D *hist_dphi_W_T1_dphi_W_T2 = new TH2D("dphi_W_T1_dphi_W_T2","dphi_W_T1_dphi_W_T2",50, 0., 2.*TMath::Pi(), 50, 0.,2.*TMath::Pi());
  TH2D *hist_dphi_W_T1_dphi_W_T2_top = new TH2D("dphi_W_T1_dphi_W_T2_top","dphi_W_T1_dphi_W_T2_top",50, 0., 2.*TMath::Pi(), 50, 0.,2.*TMath::Pi());
  TH2D *hist_dphi_W_T1_dphi_W_T2_W = new TH2D("dphi_W_T1_dphi_W_T2_W","dphi_W_T1_dphi_W_T2_W",50, 0., 2.*TMath::Pi(), 50, 0.,2.*TMath::Pi());
  TH2D *hist_ddphi_W_T1_ddphi_W_T2_top = new TH2D("ddphi_W_T1_ddphi_W_T2_top","ddphi_W_T1_ddphi_W_T2_top",50, -1., 1., 50, -1.,1.);
  TH2D *hist_ddphi_W_T1_ddphi_W_T2_W = new TH2D("ddphi_W_T1_ddphi_W_T2_W","ddphi_W_T1_ddphi_W_T2_W",50, -1., 1., 50, -1.,1.);
  TH2D *hist_dphi_W_T1_costhetaT1_top = new TH2D("dphi_W_T1_costhetaT1_top","dphi_W_T1_costhetaT1_top",50, 0., 2.*TMath::Pi(), 50, -1.,1.);
  TH2D *hist_dphi_W_T1_costhetaT1_W = new TH2D("dphi_W_T1_costhetaT1_W","dphi_W_T1_costhetaT1_W",50, 0., 2.*TMath::Pi(), 50, -1.,1.);
  TH2D *hist_dphi_T1_T2_dphi_W_T1 = new TH2D("dphi_T1_T2_dphi_W_T1","dphi_T1_T2_dphi_W_T1",50,0.,2.*TMath::Pi(),50,0.0,2.*TMath::Pi());
  TH2D *hist_dphi_T1_T2_dphi_W_T1_top = new TH2D("dphi_T1_T2_dphi_W_T1_top","dphi_T1_T2_dphi_W_T1_top",50,0.,2.*TMath::Pi(),50,0.0,2.*TMath::Pi());
  TH2D *hist_dphi_T1_T2_dphi_W_T1_W = new TH2D("dphi_T1_T2_dphi_W_T1_W","dphi_T1_T2_dphi_W_T1_W",50,0.,2.*TMath::Pi(),50,0.0,2.*TMath::Pi());
    TH2D *hist_dphi_T1_T2_costhetaT1_top = new TH2D("dphi_T1_T2_costhetaT1_top","dphi_T1_T2_costhetaT1_top",50,0.,2.*TMath::Pi(),50,-1.,1.);
    TH2D *hist_dphi_T1_T2_costhetaT1_W = new TH2D("dphi_T1_T2_costhetaT1_W","dphi_T1_T2_costhetaT1_W",50,0.,2.*TMath::Pi(),50,-1.,1.);
    TH2D *hist_dphi_T1_T2_costhetaW1_top = new TH2D("dphi_T1_T2_costhetaW1_top","dphi_T1_T2_costhetaW1_top",50,0.,2.*TMath::Pi(),50,-1.,1.);
    TH2D *hist_dphi_T1_T2_costhetaW1_W = new TH2D("dphi_T1_T2_costhetaW1_W","dphi_T1_T2_costhetaW1_W",50,0.,2.*TMath::Pi(),50,-1.,1.);
    TH2D *hist_dphi_W_T1_costhetaT2_top = new TH2D("dphi_W_T1_costhetaT2_top","dphi_W_T1_costhetaT2_top",50, 0., 2.*TMath::Pi(), 50, -1.,1.);
    TH2D *hist_dphi_W_T1_costhetaT2_W = new TH2D("dphi_W_T1_costhetaT2_W","dphi_W_T1_costhetaT2_W",50, 0., 2.*TMath::Pi(), 50, -1.,1.);
    TH2D *hist_dphi_W_T1_costhetaW1_top = new TH2D("dphi_W_T1_costhetaW1_top","dphi_W_T1_costhetaW1_top",50, 0., 2.*TMath::Pi(), 50, -1.,1.);
    TH2D *hist_dphi_W_T1_costhetaW1_W = new TH2D("dphi_W_T1_costhetaW1_W","dphi_W_T1_costhetaW1_W",50, 0., 2.*TMath::Pi(), 50, -1.,1.);
    TH2D *hist_dphi_W_T1_costhetaW2_top = new TH2D("dphi_W_T1_costhetaW1_top","dphi_W_T2_costhetaW1_top",50, 0., 2.*TMath::Pi(), 50, -1.,1.);
    TH2D *hist_dphi_W_T1_costhetaW2_W = new TH2D("dphi_W_T1_costhetaW1_W","dphi_W_T2_costhetaW1_W",50, 0., 2.*TMath::Pi(), 50, -1.,1.);
    
    TH2D *hist_Eb1_dphi_T1_T2_top = new TH2D("hist_Eb1_dphi_T1_T2_top","hist_Eb1_dphi_T1_T2_top",50, 0., 3., 50, 0., 2.*TMath::Pi());
    TH2D *hist_Eb1_dphi_T1_T2_W = new TH2D("hist_Eb1_dphi_T1_T2_W","hist_Eb1_dphi_T1_T2_W",50, 0., 3., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1_dphi_T1_T2_top = new TH2D("hist_El1_dphi_T1_T2_top","hist_El1_dphi_T1_T2_top",50, 0., 3., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1_dphi_T1_T2_W = new TH2D("hist_El1_dphi_T1_T2_W","hist_El1_dphi_T1_T2_W",50, 0., 3., 50, 0., 2.*TMath::Pi());
    TH2D *hist_MTT_dphi_T1_T2_top = new TH2D("hist_El1_dphi_T1_T2_top","hist_El1_dphi_T1_T2_top",50, 0., 3., 50, 0., 2.*TMath::Pi());
    TH2D *hist_MTT_dphi_T1_T2_W = new TH2D("hist_El1_dphi_T1_T2_W","hist_El1_dphi_T1_T2_W",50, 0., 3., 50, 0., 2.*TMath::Pi());
    TH2D *hist_Eb1_dphi_W_T1_top = new TH2D("hist_Eb1_dphi_W_T1_top","hist_Eb1_dphi_W_T1_top",50, 0., 3., 50, 0., 2.*TMath::Pi());
    TH2D *hist_Eb1_dphi_W_T1_W = new TH2D("hist_Eb1_dphi_W_T1_W","hist_Eb1_dphi_W_T1_W",50, 0., 3., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1_dphi_W_T1_top = new TH2D("hist_El1_dphi_W_T1_top","hist_El1_dphi_W_T1_top",50, 0., 3., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1_dphi_W_T1_W = new TH2D("hist_El1_dphi_W_T1_W","hist_El1_dphi_W_T1_W",50, 0., 3., 50, 0., 2.*TMath::Pi());
    TH2D *hist_MTT_dphi_W_T1_top = new TH2D("hist_El1_dphi_W_T1_top","hist_El1_dphi_W_T1_top",50, 0., 3., 50, 0., 2.*TMath::Pi());
    TH2D *hist_MTT_dphi_W_T1_W = new TH2D("hist_El1_dphi_W_T1_W","hist_El1_dphi_W_T1_W",50, 0., 3., 50, 0., 2.*TMath::Pi());
    TH2D *hist_Eb1_dphi_W_T2_top = new TH2D("hist_Eb1_dphi_W_T2_top","hist_Eb1_dphi_W_T2_top",50, 0., 3., 50, 0., 2.*TMath::Pi());
    TH2D *hist_Eb1_dphi_W_T2_W = new TH2D("hist_Eb1_dphi_W_T2_W","hist_Eb1_dphi_W_T2_W",50, 0., 3., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1_dphi_W_T2_top = new TH2D("hist_El1_dphi_W_T2_top","hist_El1_dphi_W_T2_top",50, 0., 3., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1_dphi_W_T2_W = new TH2D("hist_El1_dphi_W_T2_W","hist_El1_dphi_W_T2_W",50, 0., 3., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1Eb1_dphi_T1_T2_top = new TH2D("hist_El1Eb1_dphi_T1_T2_top","histEl1Eb1_dphi_T1_T2_top",50, 0., 1., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1Eb1_dphi_T1_T2_W = new TH2D("hist_El1Eb1_dphi_T1_T2_W","histEl1Eb1_dphi_T1_T2_W",50, 0., 1., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1Eb2_dphi_T1_T2_top = new TH2D("hist_El1Eb2_dphi_T1_T2_top","histEl1Eb1_dphi_T1_T2_top",50, 0., 1., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1Eb2_dphi_T1_T2_W = new TH2D("hist_El1Eb2_dphi_T1_T2_W","histEl1Eb1_dphi_T1_T2_W",50, 0., 1., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1Eb1_costhetaT1_top = new TH2D("hist_El1Eb1_costhetaT1_top","histEl1Eb1_costhetaT1_top",50, 0., 1., 50, -1.,1.);
    TH2D *hist_El1Eb1_costhetaT1_W = new TH2D("hist_El1Eb1_costhetaT1_W","histEl1Eb1_costhetaT1_W",50, 0., 1., 50, -1.,1.);
    TH2D *hist_El1Eb2_costhetaT1_top = new TH2D("hist_El1Eb2_costhetaT1_top","histEl1Eb1_costhetaT1_top",50, 0., 1., 50, -1.,1.);
    TH2D *hist_El1Eb2_costhetaT1_W = new TH2D("hist_El1Eb2_costhetaT1_W","histEl1Eb1_costhetaT1_W",50, 0., 1., 50, -1.,1.);
    TH2D *hist_El1Eb1_costhetaW1_top = new TH2D("hist_El1Eb1_costhetaW1_top","histEl1Eb1_costhetaW1_top",50, 0., 1., 50, -1.,1.);
    TH2D *hist_El1Eb1_costhetaW1_W = new TH2D("hist_El1Eb1_costhetaW1_W","histEl1Eb1_costhetaW1_W",50, 0., 1., 50, -1.,1.);
    TH2D *hist_El1Eb2_costhetaW1_top = new TH2D("hist_El1Eb2_costhetaW1_top","histEl1Eb1_costhetaW1_top",50, 0., 1., 50, -1.,1.);
    TH2D *hist_El1Eb2_costhetaW1_W = new TH2D("hist_El1Eb2_costhetaW1_W","histEl1Eb1_costhetaW1_W",50, 0., 1., 50, -1.,1.);
    TH2D *hist_El1Eb1_dphi_W_T1_top = new TH2D("hist_El1Eb1_dphi_W_T1_top","histEl1Eb1_dphi_W_T1_top",50, 0., 1., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1Eb1_dphi_W_T1_W = new TH2D("hist_El1Eb1_dphi_W_T1_W","histEl1Eb1_dphi_W_T1_W",50, 0., 1., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1Eb2_dphi_W_T1_top = new TH2D("hist_El1Eb2_dphi_W_T1_top","histEl1Eb1_dphi_W_T1_top",50, 0., 1., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1Eb2_dphi_W_T1_W = new TH2D("hist_El1Eb2_dphi_W_T1_W","histEl1Eb1_dphi_W_T1_W",50, 0., 1., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1Eb1_costhetaT2_top = new TH2D("hist_El1Eb1_costhetaT2_top","histEl1Eb1_costhetaT2_top",50, 0., 1., 50, -1.,1.);
    TH2D *hist_El1Eb1_costhetaT2_W = new TH2D("hist_El1Eb1_costhetaT2_W","histEl1Eb1_costhetaT2_W",50, 0., 1., 50, -1.,1.);
    TH2D *hist_El1Eb2_costhetaT2_top = new TH2D("hist_El1Eb2_costhetaT2_top","histEl1Eb1_costhetaT2_top",50, 0., 1., 50, -1.,1.);
    TH2D *hist_El1Eb2_costhetaT2_W = new TH2D("hist_El1Eb2_costhetaT2_W","histEl1Eb1_costhetaT2_W",50, 0., 1., 50, -1.,1.);
    TH2D *hist_El1Eb1_costhetaW2_top = new TH2D("hist_El1Eb1_costhetaW2_top","histEl1Eb1_costhetaW2_top",50, 0., 1., 50, -1.,1.);
    TH2D *hist_El1Eb1_costhetaW2_W = new TH2D("hist_El1Eb1_costhetaW2_W","histEl1Eb1_costhetaW2_W",50, 0., 1., 50, -1.,1.);
    TH2D *hist_El1Eb2_costhetaW2_top = new TH2D("hist_El1Eb2_costhetaW2_top","histEl1Eb1_costhetaW2_top",50, 0., 1., 50, -1.,1.);
    TH2D *hist_El1Eb2_costhetaW2_W = new TH2D("hist_El1Eb2_costhetaW2_W","histEl1Eb1_costhetaW2_W",50, 0., 1., 50, -1.,1.);
    TH2D *hist_El1Eb1_dphi_W_T2_top = new TH2D("hist_El1Eb1_dphi_W_T2_top","histEl1Eb1_dphi_W_T2_top",50, 0., 1., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1Eb1_dphi_W_T2_W = new TH2D("hist_El1Eb1_dphi_W_T2_W","histEl1Eb1_dphi_W_T2_W",50, 0., 1., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1Eb2_dphi_W_T2_top = new TH2D("hist_El1Eb2_dphi_W_T2_top","histEl1Eb1_dphi_W_T2_top",50, 0., 1., 50, 0., 2.*TMath::Pi());
    TH2D *hist_El1Eb2_dphi_W_T2_W = new TH2D("hist_El1Eb2_dphi_W_T2_W","histEl1Eb1_dphi_W_T2_W",50, 0., 1., 50, 0., 2.*TMath::Pi());

    TH2D *hist_Eb1_costhetaT2_top = new TH2D("hist_Eb1_costhetaT2_top","hist_Eb1_costhetaT2_top",50, 0., 3., 50, -1.,1.);
    TH2D *hist_Eb1_costhetaT2_W = new TH2D("hist_Eb1_costhetaT2_W","hist_Eb1_costhetaT2_W",50, 0., 3., 50, -1.,1.);
    TH2D *hist_Eb1_costhetaW1_top = new TH2D("hist_Eb1_costhetaW1_top","hist_Eb1_costhetaW1_top",50, 0., 3., 50, -1.,1.);
    TH2D *hist_Eb1_costhetaW1_W = new TH2D("hist_Eb1_costhetaW1_W","hist_Eb1_costhetaW1_W",50, 0., 3., 50, -1.,1.);
    TH2D *hist_Eb1_costhetaW2_top = new TH2D("hist_Eb1_costhetaW2_top","hist_Eb1_costhetaW2_top",50, 0., 3., 50, -1.,1.);
    TH2D *hist_Eb1_costhetaW2_W = new TH2D("hist_Eb1_costhetaW2_W","hist_Eb1_costhetaW2_W",50, 0., 3., 50, -1.,1.);
    TH2D *hist_El1_costhetaT2_top = new TH2D("hist_El1_costhetaT2_top","hist_El1_costhetaT2_top",50, 0., 3., 50, -1.,1.);
    TH2D *hist_El1_costhetaT2_W = new TH2D("hist_El1_costhetaT2_W","hist_El1_costhetaT2_W",50, 0., 3., 50, -1.,1.);
    TH2D *hist_El1_costhetaW2_top = new TH2D("hist_El1_costhetaW2_top","hist_El1_costhetaW2_top",50, 0., 3., 50, -1.,1.);
    TH2D *hist_El1_costhetaW2_W = new TH2D("hist_El1_costhetaW2_W","hist_El1_costhetaW2_W",50, 0., 3., 50, -1.,1.);

    TH2D *hist_costhetaT1_costhetaW1_top = new TH2D("costhetaT1_costhetaW1_top","costhetaT1_costhetaW1_top",50, -1., 1., 50, -1.,1.);
    TH2D *hist_costhetaT1_costhetaW2_top = new TH2D("costhetaT1_costhetaW2_top","costhetaT1_costhetaW2_top",50, -1., 1., 50, -1.,1.);
    TH2D *hist_costhetaT1_costhetaW1_W = new TH2D("costhetaT1_costhetaW1_W","costhetaT1_costhetaW1_W",50, -1., 1., 50, -1.,1.);
    TH2D *hist_costhetaT1_costhetaW2_W = new TH2D("costhetaT1_costhetaW2_W","costhetaT1_costhetaW2_W",50, -1., 1., 50, -1.,1.);
    TH2D *hist_costhetaT1_costhetaW1 = new TH2D("costhetaT1_costhetaW1","costhetaT1_costhetaW1",50, -1., 1., 50, -1.,1.);
    TH2D *hist_costhetaT1_costhetaW2 = new TH2D("costhetaT1_costhetaW2","costhetaT1_costhetaW2",50, -1., 1., 50, -1.,1.);
    TH2D *hist_dphi_T1_W1_costhetaW1 = new TH2D("dphi_T1_W1_costhetaW1","dphi_T1_W1_costhetaW1",50, 0.0, TMath::Pi(), 50, -1.,1.);
    TH2D *hist_dphi_T1_W1_costhetaW2 = new TH2D("dphi_T1_W1_costhetaW2","dphi_T1_W1_costhetaW2",50,0.0, TMath::Pi(), 50, -1.,1.);
    TH2D *hist_dphi_T1_W1_costhetaT1 = new TH2D("dphi_T1_W1_costhetaT1","dphi_T1_W1_costhetaT1",50, 0.0, TMath::Pi(), 50, -1.,1.);
    TH2D *hist_dphi_T1_W1_costhetaT2 = new TH2D("dphi_T1_W1_costhetaT2","dphi_T1_W1_costhetaT2",50, 0.0, TMath::Pi(), 50, -1.,1.);

  //
  // Generate fake events taking flat ME's for all decay angles
  //
	
  //here, we set up the stuff for dynamic gamma
  TH1D *h_gamma = (TH1D*) new TH1D("newgamma","newgamma",500,1.0,6.0);
  if(b_gamma){
    cout << "generating gamma distribution" << endl;
    for(int ibin = 1; ibin <= 500; ibin++){
      double g = h_gamma->GetBinCenter(ibin);
      double entry = Calc_dsigma_dgamma(g);
    
      if(entry > 0.)
	h_gamma->SetBinContent(ibin, entry);
    }
    cout << "done" << endl;
  }
		
  for(int i = 0; i < N; i++){
    if(b_gamma){
      gamma = h_gamma->GetRandom();
    }
   	
    //////////////////////////////////////////////////////////////
    // BEGIN CODE TO generate toy GEN LEVEL events
    //////////////////////////////////////////////////////////////
    double Mt1 = 175.;
    double MW1 = 80.;
    double Mt2 = 175.;
    double MW2 = 80.;
		
    double Mnu1 = 0.;
    double Mnu2 = 0.;
		
    GetLepTopHem(0, Mt1, MW1, 5., 0.005, Mnu1);
    GetLepTopHem(1, Mt2, MW2, 5., 0.005, Mnu2);
		
    double EB1 = (Mt1*Mt1 - MW1*MW1)/(2.*Mt1);
    double EB2 = (Mt2*Mt2 - MW2*MW2)/(2.*Mt2);
		
    double EVIS1 = (Mt1*Mt1 - Mnu1*Mnu1)/(Mt1);
    double EVIS2 = (Mt2*Mt2 - Mnu2*Mnu2)/(Mt2);
		
    double EW1 = (Mt1*Mt1 + MW1*MW1)/(2.*Mt1);
    double EW2 = (Mt2*Mt2 + MW2*MW2)/(2.*Mt2);
		
    double EL1 = (MW1*MW1-Mnu1*Mnu1)/(2.*MW1);
    double EL2 = (MW2*MW2-Mnu2*Mnu2)/(2.*MW2);
		
    double gamma1 = EW1/MW1;
    double gamma2 = EW2/MW2;
    double beta1 = 1.;
    double beta2 = 1.;
		
    //put them in the lab frame
    BoostToLabFrame(gamma);
		
    double g_eff = (P[0]+P[1]).M()/(2.*sqrt(P[0].M()*P[1].M()));
		
    PtBoost(PT);
    //////////////////////////////////////////////////////////////
    // END CODE TO generate toy GEN LEVEL event
    //////////////////////////////////////////////////////////////
        
    //////////////////////////////////////////////////////////////
    // BEGIN CODE TO analyze GEN LEVEL events
    //////////////////////////////////////////////////////////////



    // LL - This is where I can start to replace things


    
    //////////////////////////////////////////////////////////////
    // First, we calculate approximate neutrino 4-vectors
    // in W frames using two different strategies
    //////////////////////////////////////////////////////////////
    
    RJTool->newEvent();

    //RJTool->addTruParticle();

    // give b's to the tool
    RJTool->addVisParticle("b",C_2[0]);
    RJTool->addVisParticle("b",C_2[1]);

    // give leptons to the tool
    RJTool->addVisParticle("l",C_1_1[0]);
    RJTool->addVisParticle("l",C_1_1[1]);

    TVector3 MET = (vMiss[0]+vMiss[1]).Vect();
    MET.SetZ(0.0);

    RJTool->addMET( MET );

    RJTool->setHemisphereMode(0); //top symmetry

    RJTool->guessInvParticles();
    RJTool->getObservables();





    // TLorentzVector NU1_top, NU2_top, NU1_W, NU2_W;
        
    // //first, drawing symmetry between the tops
    // for(;;){
    //   TVector3 MET = (vMiss[0]+vMiss[1]).Vect();
    //   MET.SetZ(0.0);
    //   //get b-quarks in lab frame
    //   TLorentzVector B1, B2;
    //   B1 = C_2[0];
    //   B2 = C_2[1];
    //   //get leptons in lab frame
    //   TLorentzVector L1, L2;
    //   L1 = C_1_1[0];
    //   L2 = C_1_1[1];
		
    //   //We separate the objects into decay hemispheres
    //   TLorentzVector H1 = (B1+L1);
    //   TLorentzVector H2 = (B2+L2);
            
    //   //Now, we must perform longitudinal boost from lab frame to CMz frame
    //   TVector3 BL = H1.Vect()+H2.Vect();
    //   BL.SetX(0.0);
    //   BL.SetY(0.0);
    //   BL = (1./(H1.E()+H2.E()))*BL;
		
    //   B1.Boost(-BL);
    //   L1.Boost(-BL);
    //   B2.Boost(-BL);
    //   L2.Boost(-BL);
    //   H1.Boost(-BL);
    //   H2.Boost(-BL);
		 
    //   //Now, we need to 'guess' the invariant mass of the weakly interacting system
    //   double HM2min = H1.M2();
    //   double HM2max = H2.M2();
    //   if(HM2min > HM2max){
    //   	double temp = HM2max;
    //   	HM2max = HM2min;
    //   	HM2min = temp;
    //   }
	

    //   double Minv2 = (H1+H2).M2() - 4.*HM2min;
      

    //   TVector3 PT = MET+(H1+H2).Vect();
    //   double Einv2 = MET.Mag2() + Minv2;
    //   PT.SetZ(0.0);
    //   TVector3 BT = (1./( (H1+H2).E() + sqrt(Einv2) ))*PT;
    //   B1.Boost(-BT);
    //   L1.Boost(-BT);
    //   B2.Boost(-BT);
    //   L2.Boost(-BT);
    //   H1.Boost(-BT);
    //   H2.Boost(-BT);
            
    //   //Now, in CM approx frame
    //   double c = ( H1.E()+H2.E() + sqrt( (H1.E()+H2.E())*(H1.E()+H2.E()) -(H1+H2).M2() + Minv2 ))/(H1.E()+H2.E())/2.;
            
    //   double Enu1 = (c-1.)*H1.E() + c*H2.E();
    //   double Enu2 = c*H1.E() + (c-1.)*H2.E();
    //   TVector3 pnu1 = (c-1.)*H1.Vect() - c*H2.Vect();
    //   TVector3 pnu2 = (c-1.)*H2.Vect() - c*H1.Vect();
            
    //   //define neutrino 4-vectors
    //   NU1_top.SetPxPyPzE(pnu1.X(),pnu1.Y(),pnu1.Z(),Enu1);
    //   NU2_top.SetPxPyPzE(pnu2.X(),pnu2.Y(),pnu2.Z(),Enu2);

    //   //now, go back to lab frame
    //   NU1_top.Boost(BT);
    //   NU2_top.Boost(BT);
    //   NU1_top.Boost(BL);
    //   NU2_top.Boost(BL);
            
    //   break;
    // }
		
    // //second, drawing symmetry between the Ws
    // for(;;){
    //   TVector3 MET = (vMiss[0]+vMiss[1]).Vect();
    //   MET.SetZ(0.0);
    //   //get b-quarks in lab frame
    //   TLorentzVector B1, B2;
    //   B1 = C_2[0];
    //   B2 = C_2[1];
    //   //get leptons in lab frame
    //   TLorentzVector L1, L2;
    //   L1 = C_1_1[0];
    //   L2 = C_1_1[1];
            
    //   //We separate the objects into decay hemispheres
    //   TLorentzVector H1 = (B1+L1);
    //   TLorentzVector H2 = (B2+L2);
            
    //   //Now, we must perform longitudinal boost from lab frame to CMz frame
    //   TVector3 BL = H1.Vect()+H2.Vect();
    //   BL.SetX(0.0);
    //   BL.SetY(0.0);
    //   BL = (1./(H1.E()+H2.E()))*BL;
            
    //   B1.Boost(-BL);
    //   L1.Boost(-BL);
    //   B2.Boost(-BL);
    //   L2.Boost(-BL);
    //   H1.Boost(-BL);
    //   H2.Boost(-BL);
            
    //   //Now, we need to 'guess' the invariant mass of the weakly interacting system
    //   double Minv2 = (L1+L2).M2();
    //   TVector3 PT = MET+(L1+L2).Vect();
    //   double Einv2 = MET.Mag2() + Minv2;
    //   PT.SetZ(0.0);
    //   TVector3 BT = (1./( (L1+L2).E() + sqrt(Einv2) ))*PT;
    //   L1.Boost(-BT);
    //   L2.Boost(-BT);
            
    //   //Now, in CM approx frame _of WW system_
    //   double c = ( L1.E()+L2.E() + sqrt( (L1.E()+L2.E())*(L1.E()+L2.E()) -(L1+L2).M2() + Minv2 ))/(L1.E()+L2.E())/2.;
            
    //   double Enu1 = (c-1.)*L1.E() + c*L2.E();
    //   double Enu2 = c*L1.E() + (c-1.)*L2.E();
    //   TVector3 pnu1 = (c-1.)*L1.Vect() - c*L2.Vect();
    //   TVector3 pnu2 = (c-1.)*L2.Vect() - c*L1.Vect();
            
    //   //define neutrino 4-vectors
    //   NU1_W.SetPxPyPzE(pnu1.X(),pnu1.Y(),pnu1.Z(),Enu1);
    //   NU2_W.SetPxPyPzE(pnu2.X(),pnu2.Y(),pnu2.Z(),Enu2);
            
    //   //now, go back to lab frame
    //   NU1_W.Boost(BT);
    //   NU2_W.Boost(BT);
    //   NU1_W.Boost(BL);
    //   NU2_W.Boost(BL);
            
    //   break;
    // }
        
    // TLorentzVector B1, B2;
    // B1 = C_2[0];
    // B2 = C_2[1];
    // //get leptons in lab frame
    // TLorentzVector L1, L2;
    // L1 = C_1_1[0];
    // L2 = C_1_1[1];
    // //get truth neutrinos in lab frame
    // TLorentzVector NU1, NU2;
    // NU1 = C_1_2[0];
    // NU2 = C_1_2[1];

    // //Now we have two different sets of neutrino guesses in lab frame - define matching vis particles
    // TLorentzVector B1_top = B1;
    // TLorentzVector B2_top = B2;
    // TLorentzVector L1_top = L1;
    // TLorentzVector L2_top = L2;
    // TLorentzVector B1_W = B1;
    // TLorentzVector B2_W = B2;
    // TLorentzVector L1_W = L1;
    // TLorentzVector L2_W = L2;
		
    // TLorentzVector H1 = L1+B1+NU1;
    // TLorentzVector H2 = L2+B2+NU2;
    // TLorentzVector H1_top = L1+B1+NU1_top;
    // TLorentzVector H2_top = L2+B2+NU2_top;
    // TLorentzVector H1_W = L1+B1+NU1_W;
    // TLorentzVector H2_W = L2+B2+NU2_W;

    // //In the lab frame - let's move now to the ttbar CM frame
    // //boosts
    // TVector3 BltoCM = (H1+H2).BoostVector();
    // TVector3 BltoCM_top = (H1_top+H2_top).BoostVector();
    // TVector3 BltoCM_W = (H1_W+H2_W).BoostVector();

    // H1.Boost(-BltoCM);
    // H2.Boost(-BltoCM);
    // L1.Boost(-BltoCM);
    // L2.Boost(-BltoCM);
    // NU1.Boost(-BltoCM);
    // NU2.Boost(-BltoCM);
    // B1.Boost(-BltoCM);
    // B2.Boost(-BltoCM);
    // H1_top.Boost(-BltoCM_top);
    // H2_top.Boost(-BltoCM_top);
    // L1_top.Boost(-BltoCM_top);
    // L2_top.Boost(-BltoCM_top);
    // B1_top.Boost(-BltoCM_top);
    // B2_top.Boost(-BltoCM_top);
    // NU1_top.Boost(-BltoCM_top);
    // NU2_top.Boost(-BltoCM_top);
    // H1_W.Boost(-BltoCM_W);
    // H2_W.Boost(-BltoCM_W);
    // L1_W.Boost(-BltoCM_W);
    // L2_W.Boost(-BltoCM_W);
    // B1_W.Boost(-BltoCM_W);
    // B2_W.Boost(-BltoCM_W);
    // NU1_W.Boost(-BltoCM_W);
    // NU2_W.Boost(-BltoCM_W);

    // //angles
    // double costhetaTT = fabs( H1.Vect().Unit().Dot( BltoCM.Unit() ));
    // double dphiTT = fabs(H1.Vect().DeltaPhi( BltoCM ));
    // double dphiM  = fabs( (B1+B2+L1+L2).Vect().DeltaPhi( BltoCM ));
    // double costhetaTT_top = fabs( H1_top.Vect().Unit().Dot( BltoCM_top.Unit() ));
    // double dphiTT_top = fabs(H1_top.Vect().DeltaPhi( BltoCM_top ));
    // double dphiM_top  = fabs( (B1_top+B2_top+L1_top+L2_top).Vect().DeltaPhi( BltoCM_top ));
    // double costhetaTT_W = fabs( H1_W.Vect().Unit().Dot( BltoCM_W.Unit() ));
    // double dphiTT_W = fabs(H1_W.Vect().DeltaPhi( BltoCM_W ));
    // double dphiM_W  = fabs( (B1_W+B2_W+L1_W+L2_W).Vect().DeltaPhi( BltoCM_W ));
    // //scale
    // double MTT = (H1+H2).M();
    // double MTT_top = (H1_top+H2_top).M();
    // double MTT_W = (H1_W+H2_W).M();
		
    // //2D plots of angles
    // hist_dphiTT_true_v_top->Fill(dphiTT,dphiTT_top);
    // hist_dphiTT_true_v_W->Fill(dphiTT,dphiTT_W);
    // hist_dphiM_true_v_top->Fill(dphiM,dphiM_top);
    // hist_dphiM_true_v_W->Fill(dphiM,dphiM_W);
    // hist_dphiTT_dphiM_true->Fill(dphiTT,dphiM);
    // hist_dphiTT_dphiM_top->Fill(dphiTT_top,dphiM_top);
    // hist_dphiTT_dphiM_W->Fill(dphiTT_W,dphiM_W);

    // //2D plots of scale and angle
    // hist_MTT_dphiM_top->Fill(MTT_top/MTT, dphiM_top);
    // hist_MTT_dphiM_W->Fill(MTT_W/MTT, dphiM_W);
    
    // //vector normal to decay plane of CM frame
    // TVector3 vNORM_CM_T1 = B1.Vect().Cross((L1+NU1).Vect());
    // TVector3 vNORM_CM_T2 = B2.Vect().Cross((L2+NU2).Vect());
    // double dphi_T1_T2 = vNORM_CM_T1.Angle(vNORM_CM_T2);
    // double dot_dphi_T1_T2 = B1.Vect().Dot(vNORM_CM_T2);
    // if(dot_dphi_T1_T2 < 0.0 && dphi_T1_T2 > 0.0){
    //   dphi_T1_T2 = TMath::Pi()*2. - dphi_T1_T2;
    // }
    // TVector3 vNORM_CM_T1_top = B1_top.Vect().Cross((L1_top+NU1_top).Vect());
    // TVector3 vNORM_CM_T2_top = B2_top.Vect().Cross((L2_top+NU2_top).Vect());
    // double dphi_T1_T2_top = vNORM_CM_T1_top.Angle(vNORM_CM_T2_top);
    // double dot_dphi_T1_T2_top = B1_top.Vect().Dot(vNORM_CM_T2_top);
    // if(dot_dphi_T1_T2_top < 0.0 && dphi_T1_T2_top > 0.0){
    //   dphi_T1_T2_top = TMath::Pi()*2. - dphi_T1_T2_top;
    // }
    // TVector3 vNORM_CM_T1_W = B1_W.Vect().Cross((L1_W+NU1_W).Vect());
    // TVector3 vNORM_CM_T2_W = B2_W.Vect().Cross((L2_W+NU2_W).Vect());
    // double dphi_T1_T2_W = vNORM_CM_T1_W.Angle(vNORM_CM_T2_W);
    // double dot_dphi_T1_T2_W = B1.Vect().Dot(vNORM_CM_T2_W);
    // if(dot_dphi_T1_T2_W < 0.0 && dphi_T1_T2_W > 0.0){
    //   dphi_T1_T2_W = TMath::Pi()*2. - dphi_T1_T2_W;
    // }

    // double ddphi_T1_T2_top = sin(dphi_T1_T2_top-dphi_T1_T2);
    // double ddphi_T1_T2_W = sin(dphi_T1_T2_W-dphi_T1_T2);

    // //To the next frames!!!!
    // //now, to 'top' CM frame approxs
    // TVector3 BCMtoT1 = H1.BoostVector();
    // TVector3 BCMtoT2 = H2.BoostVector();
    // TVector3 BCMtoT1_top = H1_top.BoostVector();
    // TVector3 BCMtoT2_top = H2_top.BoostVector();
    // TVector3 BCMtoT1_W = H1_W.BoostVector();
    // TVector3 BCMtoT2_W = H2_W.BoostVector();

    // B1.Boost(-BCMtoT1);
    // B2.Boost(-BCMtoT2);
    // L1.Boost(-BCMtoT1);
    // L2.Boost(-BCMtoT2);
    // NU1.Boost(-BCMtoT1);
    // NU2.Boost(-BCMtoT2);
    // B1_top.Boost(-BCMtoT1_top);
    // B2_top.Boost(-BCMtoT2_top);
    // L1_top.Boost(-BCMtoT1_top);
    // L2_top.Boost(-BCMtoT2_top);
    // NU1_top.Boost(-BCMtoT1_top);
    // NU2_top.Boost(-BCMtoT2_top);
    // B1_W.Boost(-BCMtoT1_W);
    // B2_W.Boost(-BCMtoT2_W);
    // L1_W.Boost(-BCMtoT1_W);
    // L2_W.Boost(-BCMtoT2_W);
    // NU1_W.Boost(-BCMtoT1_W);
    // NU2_W.Boost(-BCMtoT2_W);
    // //decay angles in top frame
    // double costhetaT1 = B1.Vect().Unit().Dot(BCMtoT1.Unit());
    // double costhetaT2 = B2.Vect().Unit().Dot(BCMtoT2.Unit());
    // double costhetaT1_top = B1_top.Vect().Unit().Dot(BCMtoT1_top.Unit());
    // double costhetaT2_top = B2_top.Vect().Unit().Dot(BCMtoT2_top.Unit());
    // double costhetaT1_W = B1_W.Vect().Unit().Dot(BCMtoT1_W.Unit());
    // double costhetaT2_W = B2_W.Vect().Unit().Dot(BCMtoT2_W.Unit());
    // double dcosthetaT1_top = costhetaT1_top*sqrt(1.-costhetaT1*costhetaT1)-costhetaT1*sqrt(1.-costhetaT1_top*costhetaT1_top);
    // double dcosthetaT2_top = costhetaT2_top*sqrt(1.-costhetaT2*costhetaT2)-costhetaT2*sqrt(1.-costhetaT2_top*costhetaT2_top);
    // double dcosthetaT1_W = costhetaT1_W*sqrt(1.-costhetaT1*costhetaT1)-costhetaT1*sqrt(1.-costhetaT1_W*costhetaT1_W);
    // double dcosthetaT2_W = costhetaT2_W*sqrt(1.-costhetaT2*costhetaT2)-costhetaT2*sqrt(1.-costhetaT2_W*costhetaT2_W);

    // //vectors normal to decay planes of T frames
    // TVector3 vNORM_T1_B = B1.Vect().Cross(BCMtoT1);
    // TVector3 vNORM_T2_B = B2.Vect().Cross(BCMtoT2);
    // TVector3 vNORM_T1_B_top = B1_top.Vect().Cross(BCMtoT1_top);
    // TVector3 vNORM_T2_B_top = B2_top.Vect().Cross(BCMtoT2_top);
    // TVector3 vNORM_T1_B_W = B1_W.Vect().Cross(BCMtoT1_W);
    // TVector3 vNORM_T2_B_W = B2_W.Vect().Cross(BCMtoT2_W);
    // //vectors normal to W decay planes in T frames
    // TVector3 vNORM_T1_W = L1.Vect().Cross(NU1.Vect());
    // TVector3 vNORM_T2_W = L2.Vect().Cross(NU2.Vect());
    // TVector3 vNORM_T1_W_top = L1_top.Vect().Cross(NU1_top.Vect());
    // TVector3 vNORM_T2_W_top = L2_top.Vect().Cross(NU2_top.Vect());
    // TVector3 vNORM_T1_W_W = L1_W.Vect().Cross(NU1_W.Vect());
    // TVector3 vNORM_T2_W_W = L2_W.Vect().Cross(NU2_W.Vect());

    // double dphi_W_T1 = vNORM_T1_W.Angle(vNORM_T1_B);
    // double dphi_W_T2 = vNORM_T2_W.Angle(vNORM_T2_B);
    // double dot_dphi_W_T1 = L1.Vect().Dot(vNORM_T1_B);
    // double dot_dphi_W_T2 = L2.Vect().Dot(vNORM_T2_B);
    // if(dot_dphi_W_T1 < 0.0 && dphi_W_T1 > 0.0){
    //   dphi_W_T1 = TMath::Pi()*2. - dphi_W_T1;
    // }
    // if(dot_dphi_W_T2 < 0.0 && dphi_W_T2 > 0.0){
    //   dphi_W_T2 = TMath::Pi()*2. - dphi_W_T2;
    // }
    // double dphi_W_T1_top = vNORM_T1_W_top.Angle(vNORM_T1_B_top);
    // double dphi_W_T2_top = vNORM_T2_W_top.Angle(vNORM_T2_B_top);
    // double dot_dphi_W_T1_top = L1_top.Vect().Dot(vNORM_T1_B_top);
    // double dot_dphi_W_T2_top = L2_top.Vect().Dot(vNORM_T2_B_top);
    // if(dot_dphi_W_T1_top < 0.0 && dphi_W_T1_top > 0.0){
    //   dphi_W_T1_top = TMath::Pi()*2. - dphi_W_T1_top;
    // }
    // if(dot_dphi_W_T2_top < 0.0 && dphi_W_T2_top > 0.0){
    //   dphi_W_T2_top = TMath::Pi()*2. - dphi_W_T2_top;
    // }
    // double dphi_W_T1_W = vNORM_T1_W_W.Angle(vNORM_T1_B_W);
    // double dphi_W_T2_W = vNORM_T2_W_W.Angle(vNORM_T2_B_W);
    // double dot_dphi_W_T1_W = L1_W.Vect().Dot(vNORM_T1_B_W);
    // double dot_dphi_W_T2_W = L2_W.Vect().Dot(vNORM_T2_B_W);
    // if(dot_dphi_W_T1_W < 0.0 && dphi_W_T1_W > 0.0){
    //   dphi_W_T1_W = TMath::Pi()*2. - dphi_W_T1_W;
    // }
    // if(dot_dphi_W_T2_W < 0.0 && dphi_W_T2_W > 0.0){
    //   dphi_W_T2_W = TMath::Pi()*2. - dphi_W_T2_W;
    // }

    // double ddphi_W_T1_top = sin(dphi_W_T1_top-dphi_W_T1);
    // double ddphi_W_T2_top = sin(dphi_W_T2_top-dphi_W_T2);
    // double ddphi_W_T1_W = sin(dphi_W_T1_W-dphi_W_T1);
    // double ddphi_W_T2_W = sin(dphi_W_T2_W-dphi_W_T2);
    

    // //gamma for asymmetric boost
    // double gammaT = 1./pow( (1.-BCMtoT1.Mag2())*(1.-BCMtoT2.Mag2()),1./4. );
    // double gammaT_top = 1./sqrt(1.-BCMtoT1_top.Mag2());
    // double gammaT_W = 1./pow( (1.-BCMtoT1_W.Mag2())*(1.-BCMtoT2_W.Mag2()),1./4. );
    // //scale variables
    // double MT1 = H1.M();
    // double MT2 = H2.M();
    // double MT_top = H1_top.M();
    // double MT1_W = H1_W.M();
    // double MT2_W = H2_W.M();
  
    // double Eb1 = B1.E();
    // double Eb2 = B2.E();
    // double Eb1_top = B1_top.E();
    // double Eb2_top = B2_top.E();
    // double Eb1_W = B1_W.E();
    // double Eb2_W = B2_W.E();

    // //2D plots of angles
    // hist_costhetaT1_costhetaT2->Fill(costhetaT1,costhetaT2);
    // hist_costhetaT1_costhetaT2_top->Fill(costhetaT1_top,costhetaT2_top);
    // hist_costhetaT1_costhetaT2_W->Fill(costhetaT1_W,costhetaT2_W);
    // hist_dcosthetaT1_costhetaT1_top->Fill(dcosthetaT1_top,costhetaT1_top);
    // hist_dcosthetaT1_costhetaT1_W->Fill(dcosthetaT1_W,costhetaT1_W);
    // hist_dcosthetaT1_dcosthetaT2_top->Fill(dcosthetaT1_top,dcosthetaT2_top);
    // hist_dcosthetaT1_dcosthetaT2_W->Fill(dcosthetaT1_W,dcosthetaT2_W);

    // hist_dphi_W_T1_dphi_W_T2->Fill(dphi_W_T1,dphi_W_T2);
    // hist_dphi_W_T1_dphi_W_T2_top->Fill(dphi_W_T1_top,dphi_W_T2_top);
    // hist_dphi_W_T1_dphi_W_T2_W->Fill(dphi_W_T1_W,dphi_W_T2_W);
    // hist_ddphi_W_T1_ddphi_W_T2_top->Fill(ddphi_W_T1_top,ddphi_W_T2_top);
    // hist_ddphi_W_T1_ddphi_W_T2_W->Fill(ddphi_W_T1_W,ddphi_W_T2_W);

    // hist_dphi_W_T1_costhetaT1_top->Fill(dphi_W_T1_top,costhetaT1_top);
    // hist_dphi_W_T1_costhetaT1_top->Fill(dphi_W_T2_top,costhetaT2_top);
    // hist_dphi_W_T1_costhetaT1_W->Fill(dphi_W_T1_W,costhetaT1_W);
    // hist_dphi_W_T1_costhetaT1_W->Fill(dphi_W_T2_W,costhetaT2_W);

    // hist_dphi_T1_T2_dphi_W_T1->Fill(dphi_T1_T2,dphi_W_T1);
    // hist_dphi_T1_T2_dphi_W_T1_top->Fill(dphi_T1_T2_top,dphi_W_T1_top);
    // hist_dphi_T1_T2_dphi_W_T1_W->Fill(dphi_T1_T2_W,dphi_W_T1_W);
    
    
    // //2D plots scales and angles
    // hist_dcosthetaT1_gammaT_top->Fill(dcosthetaT1_top,gammaT_top);
    // hist_dcosthetaT1_gammaT_top->Fill(dcosthetaT2_top,gammaT_top);
    // hist_dcosthetaT1_gammaT_W->Fill(dcosthetaT1_W,gammaT_W);
    // hist_dcosthetaT1_gammaT_W->Fill(dcosthetaT2_W,gammaT_W);

    // //2D plots of scales
    // hist_gammaTrue_gamma_top->Fill(gammaT,gammaT_top);
    // hist_gammaTrue_gamma_W->Fill(gammaT,gammaT_W);
    // hist_Eb1_gammaT_top->Fill(Eb1_top/Eb1,gammaT_top);
    // hist_Eb1_gammaT_top->Fill(Eb2_top/Eb2,gammaT_top);
    // hist_Eb1_gammaT_W->Fill(Eb1_W/Eb1,gammaT_W);
    // hist_Eb1_gammaT_W->Fill(Eb2_W/Eb2,gammaT_W);
    // hist_Eb1_costhetaT_top->Fill(Eb1_top/Eb1,costhetaT1_top);
    // hist_Eb1_costhetaT_top->Fill(Eb2_top/Eb2,costhetaT2_top);
    // hist_Eb1_costhetaT_W->Fill(Eb1_W/Eb1,costhetaT1_W);
    // hist_Eb1_costhetaT_W->Fill(Eb2_W/Eb2,costhetaT2_W);
    // hist_Eb1_dcosthetaT_top->Fill(Eb1_top/Eb1,dcosthetaT1_top);
    // hist_Eb1_dcosthetaT_top->Fill(Eb2_top/Eb2,dcosthetaT2_top);
    // hist_Eb1_dcosthetaT_W->Fill(Eb1_W/Eb1,dcosthetaT1_W);
    // hist_Eb1_dcosthetaT_W->Fill(Eb2_W/Eb2,dcosthetaT2_W);
    // hist_Eb1_Eb2_top->Fill(Eb1_top/Eb1,Eb2_top/Eb2);
    // hist_Eb1_Eb2_W->Fill(Eb1_W/Eb1,Eb2_W/Eb2);
    // hist_MTT_Eb1_top->Fill(MTT_top/MTT,Eb1_top/Eb1);
    // hist_MTT_Eb1_top->Fill(MTT_top/MTT,Eb2_top/Eb2);
    // hist_MTT_Eb1_W->Fill(MTT_W/MTT,Eb1_W/Eb1);
    // hist_MTT_Eb1_W->Fill(MTT_W/MTT,Eb2_W/Eb2);
    // hist_MT1_MT2_W->Fill(MT1_W/MT1,MT2_W/MT2);

    // //Heading to the last frames!!!!!
    // //from the T rest frames to the W rest frames
    // TVector3 BTtoW1 = (L1+NU1).BoostVector();
    // TVector3 BTtoW2 = (L2+NU2).BoostVector();
    // TVector3 BTtoW1_top = (L1_top+NU1_top).BoostVector();
    // TVector3 BTtoW2_top = (L2_top+NU2_top).BoostVector();
    // TVector3 BTtoW1_W = (L1_W+NU1_W).BoostVector();
    // TVector3 BTtoW2_W = (L2_W+NU2_W).BoostVector();

    // L1.Boost(-BTtoW1);
    // NU1.Boost(-BTtoW1);
    // L2.Boost(-BTtoW2);
    // NU2.Boost(-BTtoW2);
    // L1_top.Boost(-BTtoW1_top);
    // NU1_top.Boost(-BTtoW1_top);
    // L2_top.Boost(-BTtoW2_top);
    // NU2_top.Boost(-BTtoW2_top);
    // L1_W.Boost(-BTtoW1_W);
    // NU1_W.Boost(-BTtoW1_W);
    // L2_W.Boost(-BTtoW2_W);
    // NU2_W.Boost(-BTtoW2_W);

    // double El1 = L1.E();
    // double El2 = L2.E();
    // double El1_top = L1_top.E();
    // double El2_top = L2_top.E();
    // double El1_W = L1_W.E();
    // double El2_W = L2_W.E();
    
    // //calculate some angles
    // double costhetaW1 = L1.Vect().Unit().Dot(BTtoW1.Unit());
    // double costhetaW2 = L2.Vect().Unit().Dot(BTtoW2.Unit());
    // double costhetaW1_top = L1_top.Vect().Unit().Dot(BTtoW1_top.Unit());
    // double costhetaW2_top = L2_top.Vect().Unit().Dot(BTtoW2_top.Unit());
    // double costhetaW1_W = L1_W.Vect().Unit().Dot(BTtoW1_W.Unit());
    // double costhetaW2_W = L2_W.Vect().Unit().Dot(BTtoW2_W.Unit());
    // double dcosthetaW1_top = costhetaW1*sqrt(1.-costhetaW1_top*costhetaW1_top) - costhetaW1_top*sqrt(1.-costhetaW1*costhetaW1);
    // double dcosthetaW2_top = costhetaW2*sqrt(1.-costhetaW2_top*costhetaW2_top) - costhetaW2_top*sqrt(1.-costhetaW2*costhetaW2);
    // double dcosthetaW1_W = costhetaW1*sqrt(1.-costhetaW1_W*costhetaW1_W) - costhetaW1_W*sqrt(1.-costhetaW1*costhetaW1);
    // double dcosthetaW2_W = costhetaW2*sqrt(1.-costhetaW2_W*costhetaW2_W) - costhetaW2_W*sqrt(1.-costhetaW2*costhetaW2);

    // //2D plots of angles
    // hist_costhetaW1_costhetaW2->Fill(costhetaW1,costhetaW2);
    // hist_costhetaW1_costhetaW2_top->Fill(costhetaW1_top,costhetaW2_top);
    // hist_costhetaW1_costhetaW2_W->Fill(costhetaW1_W,costhetaW2_W);
    // hist_dcosthetaW1_costhetaW1_top->Fill(dcosthetaW1_top,costhetaW1_top);
    // hist_dcosthetaW1_costhetaW1_W->Fill(dcosthetaW1_W,costhetaW1_W);
    // hist_dcosthetaW1_dcosthetaW2_top->Fill(dcosthetaW1_top,dcosthetaW2_top);
    // hist_dcosthetaW1_dcosthetaW2_W->Fill(dcosthetaW1_W,dcosthetaW2_W);

    // //2D plots of scale
    // hist_El1_gammaT_top->Fill(El1_top/El1,gammaT_top);
    // hist_El1_gammaT_top->Fill(El2_top/El2,gammaT_top);
    // hist_El1_gammaT_W->Fill(El1_W/El1,gammaT_W);
    // hist_El1_gammaT_W->Fill(El2_W/El2,gammaT_W);
    // hist_El1_costhetaW_top->Fill(El1_top/El1,costhetaW1_top);
    // hist_El1_costhetaW_top->Fill(El2_top/El2,costhetaW2_top);
    // hist_El1_costhetaW_W->Fill(El1_W/El1,costhetaW1_W);
    // hist_El1_costhetaW_W->Fill(El2_W/El2,costhetaW2_W);
    // hist_El1_dcosthetaW_top->Fill(El1_top/El1,dcosthetaW1_top);
    // hist_El1_dcosthetaW_top->Fill(El2_top/El2,dcosthetaW2_top);
    // hist_El1_dcosthetaW_W->Fill(El1_W/El1,dcosthetaW1_W);
    // hist_El1_dcosthetaW_W->Fill(El2_W/El2,dcosthetaW2_W);
    // hist_El1_El2_top->Fill(El1_top/El1,El2_top/El2);
    // hist_El1_El2_W->Fill(El1_W/El1,El2_W/El2);
    // hist_MTT_El1_top->Fill(MTT_top/MTT,El1_top/El1);
    // hist_MTT_El1_top->Fill(MTT_top/MTT,El2_top/El2);
    // hist_MTT_El1_W->Fill(MTT_W/MTT,El1_W/El1);
    // hist_MTT_El1_W->Fill(MTT_W/MTT,El2_W/El2);
    // hist_Eb1_El1_top->Fill(Eb1_top/Eb1,El1_top/El1);
    // hist_Eb1_El1_top->Fill(Eb2_top/Eb2,El2_top/El2);
    // hist_Eb1_El1_W->Fill(Eb1_W/Eb1,El1_W/El1);
    // hist_Eb1_El1_W->Fill(Eb2_W/Eb2,El2_W/El2);
    // hist_Eb1_El2_top->Fill(Eb1_top/Eb1,El2_top/El2);
    // hist_Eb1_El2_top->Fill(Eb2_top/Eb2,El1_top/El1);
    // hist_Eb1_El2_W->Fill(Eb1_W/Eb1,El2_W/El2);
    // hist_Eb1_El2_W->Fill(Eb2_W/Eb2,El1_W/El1);
    // hist_El1_Mt1_top->Fill(El1_top/El1,MT_top/MT1);
    // hist_El1_Mt1_top->Fill(El2_top/El2,MT_top/MT2);
    // hist_El1_Mt1_W->Fill(El1_W/El1,MT1_W/MT1);
    // hist_El1_Mt1_W->Fill(El2_W/El2,MT2_W/MT2);
    // hist_El1_costhetaT1_top->Fill(El1_top/El1,costhetaT1_top);
    // hist_El1_costhetaT1_top->Fill(El2_top/El2,costhetaT2_top);
    // hist_El1_costhetaT1_W->Fill(El1_W/El1,costhetaT1_W);
    // hist_El1_costhetaT1_W->Fill(El2_W/El2,costhetaT2_W);
      
    //   hist_dphi_T1_T2_costhetaT1_top->Fill(dphi_T1_T2_top,costhetaT1_top);
    //   hist_dphi_T1_T2_costhetaT1_top->Fill(dphi_T1_T2_top,costhetaT2_top);
    //   hist_dphi_T1_T2_costhetaT1_W->Fill(dphi_T1_T2_W,costhetaT1_W);
    //   hist_dphi_T1_T2_costhetaT1_W->Fill(dphi_T1_T2_W,costhetaT2_W);
    //   hist_dphi_T1_T2_costhetaW1_top->Fill(dphi_T1_T2_top,costhetaW1_top);
    //   hist_dphi_T1_T2_costhetaW1_top->Fill(dphi_T1_T2_top,costhetaW2_top);
    //   hist_dphi_T1_T2_costhetaW1_W->Fill(dphi_T1_T2_W,costhetaW1_W);
    //   hist_dphi_T1_T2_costhetaW1_W->Fill(dphi_T1_T2_W,costhetaW2_W);
    //   hist_dphi_W_T1_costhetaT2_top->Fill(dphi_W_T1_top,costhetaT2_top);
    //   hist_dphi_W_T1_costhetaT2_top->Fill(dphi_W_T2_top,costhetaT1_top);
    //   hist_dphi_W_T1_costhetaT2_W->Fill(dphi_W_T1_W,costhetaT2_W);
    //   hist_dphi_W_T1_costhetaT2_W->Fill(dphi_W_T2_W,costhetaT1_W);
    //   hist_dphi_W_T1_costhetaW1_top->Fill(dphi_W_T1_top,costhetaW1_top);
    //   hist_dphi_W_T1_costhetaW1_top->Fill(dphi_W_T2_top,costhetaW2_top);
    //   hist_dphi_W_T1_costhetaW1_W->Fill(dphi_W_T1_W,costhetaW1_W);
    //   hist_dphi_W_T1_costhetaW1_W->Fill(dphi_W_T2_W,costhetaW2_W);
    //   hist_dphi_W_T1_costhetaW2_top->Fill(dphi_W_T1_top,costhetaW2_top);
    //   hist_dphi_W_T1_costhetaW2_top->Fill(dphi_W_T2_top,costhetaW1_top);
    //   hist_dphi_W_T1_costhetaW2_W->Fill(dphi_W_T1_W,costhetaW2_W);
    //   hist_dphi_W_T1_costhetaW2_W->Fill(dphi_W_T2_W,costhetaW1_W);

    //   hist_Eb1_dphi_T1_T2_top->Fill(Eb1_top/Eb1,dphi_T1_T2_top);
    // hist_Eb1_dphi_T1_T2_W->Fill(Eb1_W/Eb1,dphi_T1_T2_W);
    // hist_El1_dphi_T1_T2_top->Fill(El1_top/El1,dphi_T1_T2_top);
    // hist_El1_dphi_T1_T2_W->Fill(El1_W/El1,dphi_T1_T2_W);
    // hist_MTT_dphi_T1_T2_top->Fill(MTT_top/MTT,dphi_T1_T2_top);
    // hist_MTT_dphi_T1_T2_W->Fill(MTT_W/MTT,dphi_T1_T2_W);
    // hist_Eb1_dphi_W_T1_top->Fill(Eb1_top/Eb1,dphi_W_T1_top);
    // hist_Eb1_dphi_W_T1_W->Fill(Eb1_W/Eb1,dphi_W_T1_W);
    // hist_El1_dphi_W_T1_top->Fill(El1_top/El1,dphi_W_T1_top);
    // hist_El1_dphi_W_T1_W->Fill(El1_W/El1,dphi_W_T1_W);
    // hist_MTT_dphi_W_T1_top->Fill(MTT_top/MTT,dphi_W_T1_top);
    // hist_MTT_dphi_W_T1_W->Fill(MTT_W/MTT,dphi_W_T1_W);
    // hist_Eb1_dphi_W_T2_top->Fill(Eb1_top/Eb1,dphi_W_T2_top);
    // hist_Eb1_dphi_W_T2_W->Fill(Eb1_W/Eb1,dphi_W_T2_W);
    // hist_El1_dphi_W_T2_top->Fill(El1_top/El1,dphi_W_T2_top);
    // hist_El1_dphi_W_T2_W->Fill(El1_W/El1,dphi_W_T2_W);
    // double El1Eb1_top = El1_top/(El1_top+Eb1_top);
    // double El1Eb1_W = El1_W/(El1_W+Eb1_W);
    // double El1Eb2_top = El1_top/(El1_top+Eb2_top);
    // double El1Eb2_W = El1_W/(El1_W+Eb2_W);
    // hist_El1Eb1_dphi_T1_T2_top->Fill(El1Eb1_top,dphi_T1_T2_top);
    // hist_El1Eb1_dphi_T1_T2_W->Fill(El1Eb1_W,dphi_T1_T2_W);
    // hist_El1Eb2_dphi_T1_T2_top->Fill(El1Eb2_top,dphi_T1_T2_top);
    // hist_El1Eb2_dphi_T1_T2_W->Fill(El1Eb2_W,dphi_T1_T2_W);
    // hist_El1Eb1_costhetaT1_top->Fill(El1Eb1_top,costhetaT1_top);
    // hist_El1Eb1_costhetaT1_W->Fill(El1Eb1_W,costhetaT1_W); 
    // hist_El1Eb2_costhetaT1_top->Fill(El1Eb2_top,costhetaT1_top);
    // hist_El1Eb2_costhetaT1_W->Fill(El1Eb2_W,costhetaT1_W);
    // hist_El1Eb1_costhetaW1_top->Fill(El1Eb1_top,costhetaW1_top);
    // hist_El1Eb1_costhetaW1_W->Fill(El1Eb1_W,costhetaW1_W); 
    // hist_El1Eb2_costhetaW1_top->Fill(El1Eb2_top,costhetaW1_top);
    // hist_El1Eb2_costhetaW1_W->Fill(El1Eb2_W,costhetaW1_W); 
    // hist_El1Eb1_dphi_W_T1_top->Fill(El1Eb1_top,dphi_W_T1_top); 
    // hist_El1Eb1_dphi_W_T1_W->Fill(El1Eb1_W,dphi_W_T1_W); 
    // hist_El1Eb2_dphi_W_T1_top->Fill(El1Eb2_top,dphi_W_T1_top);
    // hist_El1Eb2_dphi_W_T1_W->Fill(El1Eb2_W,dphi_W_T1_W); 
    // hist_El1Eb1_costhetaT2_top->Fill(El1Eb1_top,costhetaT2_top);
    // hist_El1Eb1_costhetaT2_W->Fill(El1Eb1_W,costhetaT2_W); 
    // hist_El1Eb2_costhetaT2_top->Fill(El1Eb2_top,costhetaT2_top);
    // hist_El1Eb2_costhetaT2_W->Fill(El1Eb2_W,costhetaT2_W); 
    // hist_El1Eb1_costhetaW2_top->Fill(El1Eb1_top,costhetaW2_top);
    // hist_El1Eb1_costhetaW2_W->Fill(El1Eb1_W,costhetaW2_W); 
    // hist_El1Eb2_costhetaW2_top->Fill(El1Eb2_top,costhetaW2_top);
    // hist_El1Eb2_costhetaW2_W->Fill(El1Eb2_W,costhetaW2_W); 
    // hist_El1Eb1_dphi_W_T2_top->Fill(El1Eb1_top,dphi_W_T2_top); 
    // hist_El1Eb1_dphi_W_T2_W->Fill(El1Eb1_W,dphi_W_T2_W); 
    // hist_El1Eb2_dphi_W_T2_top->Fill(El1Eb2_top,dphi_W_T2_top);
    // hist_El1Eb2_dphi_W_T2_W->Fill(El1Eb2_W,dphi_W_T2_W); 

    // hist_Eb1_costhetaT2_top->Fill(Eb1_top/Eb1,costhetaT2_top);
    // hist_Eb1_costhetaT2_W->Fill(Eb1_W/Eb1,costhetaT2_W);
    // hist_Eb1_costhetaW1_top->Fill(Eb1_top/Eb1,costhetaW1_top);
    // hist_Eb1_costhetaW1_W->Fill(Eb1_W/Eb1,costhetaW1_W);
    // hist_Eb1_costhetaW2_top->Fill(Eb1_top/Eb1,costhetaW2_top);
    // hist_Eb1_costhetaW2_W->Fill(Eb1_W/Eb1,costhetaW2_W);
    // hist_El1_costhetaT2_top->Fill(El1_top/El1,costhetaT2_top);
    // hist_El1_costhetaT2_W->Fill(El1_W/El1,costhetaT2_W);
    // hist_El1_costhetaW2_top->Fill(El1_top/El1,costhetaW2_top);
    // hist_El1_costhetaW2_W->Fill(El1_W/El1,costhetaW2_W);

    // hist_costhetaT1_costhetaW1_top->Fill(costhetaT1_top,costhetaW1_top);
    // hist_costhetaT1_costhetaW2_top->Fill(costhetaT1_top,costhetaW2_top);
    // hist_costhetaT1_costhetaW1_W->Fill(costhetaT1_W,costhetaW1_W);
    // hist_costhetaT1_costhetaW2_W->Fill(costhetaT1_W,costhetaW2_W);
    // hist_costhetaT1_costhetaW1->Fill(costhetaT1,costhetaW1);
    // hist_costhetaT1_costhetaW2->Fill(costhetaT1,costhetaW2);
    // hist_dphi_T1_W1_costhetaW1->Fill(dphi_W_T1,costhetaW1);
    // hist_dphi_T1_W1_costhetaW2->Fill(dphi_W_T1,costhetaW2);
    // hist_dphi_T1_W1_costhetaT1->Fill(dphi_W_T1,costhetaT1);
    // hist_dphi_T1_W1_costhetaT2->Fill(dphi_W_T1,costhetaT2);
  }
	
  // TCanvas *c1 = Plot_Me_2D("c1", hist_dphiTT_true_v_top, "| #Delta #phi_{TT}^{true} |", "| #Delta #phi_{TT}^{top} |");
  // TCanvas *c2 = Plot_Me_2D("c2", hist_dphiTT_true_v_W, "| #Delta #phi_{TT}^{true} |", "| #Delta #phi_{TT}^{W} |");
  // TCanvas *c3 = Plot_Me_2D("c3", hist_dphiM_true_v_top, "| #Delta #phi_{M}^{true} |", "| #Delta #phi_{M}^{top} |");
  // TCanvas *c4 = Plot_Me_2D("c4", hist_dphiM_true_v_W, "| #Delta #phi_{M}^{true} |", "| #Delta #phi_{M}^{W} |");
  // TCanvas *c5 = Plot_Me_2D("c5", hist_dphiTT_dphiM_true, "| #Delta #phi_{TT}^{true} |", "| #Delta #phi_{M}^{true} |");
  // TCanvas *c6 = Plot_Me_2D("c6", hist_dphiTT_dphiM_top, "| #Delta #phi_{TT}^{top} |", "| #Delta #phi_{M}^{top} |");
  // TCanvas *c7 = Plot_Me_2D("c7", hist_dphiTT_dphiM_W, "| #Delta #phi_{TT}^{W} |", "| #Delta #phi_{M}^{W} |");
  // TCanvas *c8 = Plot_Me_2D("c8", hist_MTT_dphiM_top, "M_{tt}^{top} / M_{tt}^{true}", "| #Delta #phi_{M}^{top} |");
  // TCanvas *c9 = Plot_Me_2D("c9", hist_MTT_dphiM_W, "M_{tt}^{W} / M_{tt}^{true}", "| #Delta #phi_{M}^{W} |");
  // TCanvas *c10 = Plot_Me_2D("c10", hist_costhetaT1_costhetaT2, "cos #theta_{T1}^{true}", "cos #theta_{T2}^{true}");
  // TCanvas *c11 = Plot_Me_2D("c11", hist_costhetaT1_costhetaT2_top, "cos #theta_{T1}^{top}", "cos #theta_{T2}^{top}");
  // TCanvas *c12 = Plot_Me_2D("c12", hist_costhetaT1_costhetaT2_W, "cos #theta_{T1}^{W}", "cos #theta_{T2}^{W}");
  //TCanvas *c13 = Plot_Me_2D("c13", hist_dcosthetaT1_costhetaT1_top, "sin (#theta_{T1}^{top}-#theta_{T1}^{true})", "cos #theta_{T1}^{top}");
  //TCanvas *c14 = Plot_Me_2D("c14", hist_dcosthetaT1_costhetaT1_W, "sin (#theta_{T1}^{W}-#theta_{T1}^{true})", "cos #theta_{T1}^{W}");
  //TCanvas *c15 = Plot_Me_2D("c15", hist_dcosthetaT1_dcosthetaT2_top, "sin (#theta_{T1}^{top}-#theta_{T1}^{true})", "sin (#theta_{T2}^{top}-#theta_{T2}^{true})");
  //TCanvas *c16 = Plot_Me_2D("c16", hist_dcosthetaT1_dcosthetaT2_W, "sin (#theta_{T1}^{W}-#theta_{T1}^{true})", "sin (#theta_{T2}^{W}-#theta_{T2}^{true})");
  TCanvas *c17 = Plot_Me_2D("c17", hist_gammaTrue_gamma_top, "#gamma_{T}^{true}", "#gamma_{T}^{top}");
  TCanvas *c18 = Plot_Me_2D("c18", hist_gammaTrue_gamma_W, "#gamma_{T}^{true}", "#gamma_{T}^{W}");
  //TCanvas *c19 = Plot_Me_2D("c19", hist_dcosthetaT1_gammaT_top, "sin (#theta_{T1}^{top}-#theta_{T1}^{true})", "#gamma_{T}^{top}");
  //TCanvas *c20 = Plot_Me_2D("c20", hist_dcosthetaT1_gammaT_W, "sin (#theta_{T1}^{W}-#theta_{T1}^{true})", "#gamma_{T}^{W}");
  TCanvas *c21 = Plot_Me_2D("c21", hist_Eb1_gammaT_top, "E_{b1}^{top}/E_{b1}^{true}", "#gamma_{T}^{top}");
  TCanvas *c22 = Plot_Me_2D("c22", hist_Eb1_gammaT_W, "E_{b1}^{W}/E_{b1}^{true}", "#gamma_{T}^{W}");
  TCanvas *c23 = Plot_Me_2D("c23", hist_Eb1_costhetaT_top, "E_{b1}^{top}/E_{b1}^{true}", "cos #theta_{T1}^{top}");
  TCanvas *c24 = Plot_Me_2D("c24", hist_Eb1_costhetaT_W, "E_{b1}^{W}/E_{b1}^{true}", "cos #theta_{T1}^{W}");
  TCanvas *c25 = Plot_Me_2D("c25", hist_Eb1_dcosthetaT_top, "E_{b1}^{top}/E_{b1}^{true}", "sin (#theta_{T1}^{top}-#theta_{T1}^{true})");
  TCanvas *c26 = Plot_Me_2D("c26", hist_Eb1_dcosthetaT_W, "E_{b1}^{W}/E_{b1}^{true}", "sin (#theta_{T1}^{W}-#theta_{T1}^{true})");
  TCanvas *c27 = Plot_Me_2D("c27", hist_Eb1_Eb2_top, "E_{b1}^{top}/E_{b1}^{true}", "E_{b2}^{top}/E_{b2}^{true}");
  TCanvas *c28 = Plot_Me_2D("c28", hist_Eb1_Eb2_W, "E_{b1}^{W}/E_{b1}^{true}", "E_{b2}^{W}/E_{b2}^{true}");
  TCanvas *c29 = Plot_Me_2D("c29", hist_MTT_Eb1_top, "M_{tt}^{top}/M_{tt}", "E_{b1}^{top}/E_{b1}^{true}");
  TCanvas *c30 = Plot_Me_2D("c30", hist_MTT_Eb1_W, "M_{tt}^{W}/M_{tt}", "E_{b1}^{W}/E_{b1}^{true}");
  TCanvas *c31 = Plot_Me_2D("c31", hist_MT1_MT2_W, "M_{t1}^{W}/M_{t1}^{true}", "M_{t2}^{W}/M_{t2}^{true}");
  TCanvas *c32 = Plot_Me_2D("c32", hist_costhetaW1_costhetaW2, "cos #theta_{W1}^{true}", "cos #theta_{W2}^{true}");
  TCanvas *c33 = Plot_Me_2D("c33", hist_costhetaW1_costhetaW2_top, "cos #theta_{W1}^{top}", "cos #theta_{W2}^{top}");
  TCanvas *c34 = Plot_Me_2D("c34", hist_costhetaW1_costhetaW2_W, "cos #theta_{W1}^{W}", "cos #theta_{W2}^{W}");
  //TCanvas *c35 = Plot_Me_2D("c35", hist_dcosthetaW1_costhetaW1_top, "sin (#theta_{W1}^{top}-#theta_{W1}^{true})", "cos #theta_{W1}^{top}");
  //TCanvas *c36 = Plot_Me_2D("c36", hist_dcosthetaW1_costhetaW1_W, "sin (#theta_{W1}^{W}-#theta_{W1}^{true})", "cos #theta_{W1}^{W}");
  //TCanvas *c37 = Plot_Me_2D("c37", hist_dcosthetaW1_dcosthetaW2_top, "sin (#theta_{W1}^{top}-#theta_{W1}^{true})", "sin (#theta_{W2}^{top}-#theta_{W2}^{true})");
  //TCanvas *c38 = Plot_Me_2D("c38", hist_dcosthetaW1_dcosthetaW2_W, "sin (#theta_{W1}^{W}-#theta_{W1}^{true})", "sin (#theta_{W2}^{W}-#theta_{W2}^{true})");
  TCanvas *c39 = Plot_Me_2D("c39", hist_El1_gammaT_top, "E_{l1}^{top}/E_{l1}^{true}", "#gamma_{T}^{top}");
  TCanvas *c40 = Plot_Me_2D("c40", hist_El1_gammaT_W, "E_{l1}^{W}/E_{l1}^{true}", "#gamma_{T}^{W}");
  TCanvas *c41 = Plot_Me_2D("c41", hist_El1_costhetaW_top, "E_{l1}^{top}/E_{l1}^{true}", "cos #theta_{W1}^{top}");
  TCanvas *c42 = Plot_Me_2D("c42", hist_El1_costhetaW_W, "E_{l1}^{W}/E_{l1}^{true}", "cos #theta_{W1}^{W}");
  TCanvas *c43 = Plot_Me_2D("c43", hist_El1_dcosthetaW_top, "E_{l1}^{top}/E_{l1}^{true}", "sin (#theta_{W1}^{top}-#theta_{W1}^{true})");
  TCanvas *c44 = Plot_Me_2D("c44", hist_El1_dcosthetaW_W, "E_{l1}^{W}/E_{l1}^{true}", "sin (#theta_{W1}^{W}-#theta_{W1}^{true})");
  TCanvas *c45 = Plot_Me_2D("c45", hist_El1_El2_top, "E_{l1}^{top}/E_{l1}^{true}", "E_{l2}^{top}/E_{l2}^{true}");
  TCanvas *c46 = Plot_Me_2D("c46", hist_El1_El2_W, "E_{l1}^{W}/E_{l1}^{true}", "E_{l2}^{W}/E_{l2}^{true}");
  TCanvas *c47 = Plot_Me_2D("c47", hist_MTT_El1_top, "M_{tt}^{top}/M_{tt}", "E_{l1}^{top}/E_{l1}^{true}");
  TCanvas *c48 = Plot_Me_2D("c48", hist_MTT_El1_W, "M_{tt}^{W}/M_{tt}", "E_{l1}^{W}/E_{l1}^{true}");
  TCanvas *c49 = Plot_Me_2D("c49", hist_Eb1_El1_top, "E_{b1}^{top}/E_{b1}^{true}", "E_{l1}^{top}/E_{l1}^{true}");
  TCanvas *c50 = Plot_Me_2D("c50", hist_Eb1_El1_W, "E_{b1}^{W}/E_{b1}^{true}", "E_{l1}^{W}/E_{l1}^{true}");
  TCanvas *c51 = Plot_Me_2D("c51", hist_Eb1_El2_top, "E_{b1}^{top}/E_{b1}^{true}", "E_{l2}^{top}/E_{l2}^{true}");
  TCanvas *c52 = Plot_Me_2D("c52", hist_Eb1_El2_W, "E_{b1}^{W}/E_{b1}^{true}", "E_{l2}^{W}/E_{l2}^{true}");
  TCanvas *c53 = Plot_Me_2D("c53", hist_El1_Mt1_top, "E_{l1}^{top}/E_{l1}^{true}", "M_{t1}^{top}/M_{t1}^{true}");
  TCanvas *c54 = Plot_Me_2D("c54", hist_El1_Mt1_W, "E_{l1}^{W}/E_{l1}^{true}", "M_{t1}^{W}/M_{t1}^{true}");
  TCanvas *c55 = Plot_Me_2D("c55", hist_El1_costhetaT1_top, "E_{l1}^{top}/E_{l1}^{true}", "cos #theta_{T1}^{top}");
  TCanvas *c56 = Plot_Me_2D("c56", hist_El1_costhetaT1_W, "E_{l1}^{W}/E_{l1}^{true}", "cos #theta_{T1}^{W}");
  //TCanvas *c57 = Plot_Me_2D("c57", hist_dphi_W_T1_dphi_W_T2, "#Delta#phi_{T1,W1}^{true}", "#Delta#phi_{T2,W2}^{true}");
  TCanvas *c58 = Plot_Me_2D("c58", hist_dphi_W_T1_dphi_W_T2_top, "#Delta#phi_{T1,W1}^{top}", "#Delta#phi_{T2,W2}^{top}");
  TCanvas *c59 = Plot_Me_2D("c59", hist_dphi_W_T1_dphi_W_T2_W, "#Delta#phi_{T1,W1}^{W}", "#Delta#phi_{T2,W2}^{W}");
  TCanvas *c60 = Plot_Me_2D("c60", hist_ddphi_W_T1_ddphi_W_T2_top, "sin(#Delta#phi_{T1,W1}^{top}-#Delta#phi_{T1,W1}^{true})", 
			    "sin(#Delta#phi_{T2,W2}^{top}-#Delta#phi_{T2,W2}^{true})");
  TCanvas *c61 = Plot_Me_2D("c61", hist_ddphi_W_T1_ddphi_W_T2_W, "sin(#Delta#phi_{T1,W1}^{W}-#Delta#phi_{T1,W1}^{true})", 
			    "sin(#Delta#phi_{T2,W2}^{W}-#Delta#phi_{T2,W2}^{true})");
  TCanvas *c62 = Plot_Me_2D("c62", hist_dphi_W_T1_costhetaT1_top, "#Delta#phi_{T1,W1}^{top}", "cos #theta_{T1}^{top}");
  TCanvas *c63 = Plot_Me_2D("c63", hist_dphi_W_T1_costhetaT1_W, "#Delta#phi_{T1,W1}^{W}", "cos #theta_{T1}^{W}");
  //TCanvas *c64 = Plot_Me_2D("c64", hist_dphi_T1_T2_dphi_W_T1, "#Delta#phi_{T1,T2}^{true}", "#Delta#phi_{T1,W1}^{true}");
  TCanvas *c65 = Plot_Me_2D("c65", hist_dphi_T1_T2_dphi_W_T1_top, "#Delta#phi_{T1,T2}^{top}", "#Delta#phi_{T1,W1}^{top}");
  TCanvas *c66 = Plot_Me_2D("c66", hist_dphi_T1_T2_dphi_W_T1_W, "#Delta#phi_{T1,T2}^{W}", "#Delta#phi_{T1,W1}^{W}");
    TCanvas *c67 = Plot_Me_2D("c67", hist_dphi_T1_T2_costhetaT1_top, "#Delta#phi_{T1,T2}^{top}", "cos #theta_{T1}^{top}");
    TCanvas *c68 = Plot_Me_2D("c68", hist_dphi_T1_T2_costhetaT1_W, "#Delta#phi_{T1,T2}^{W}", "cos #theta_{T1}^{W}");
    TCanvas *c69 = Plot_Me_2D("c69", hist_dphi_T1_T2_costhetaW1_top, "#Delta#phi_{T1,T2}^{top}", "cos #theta_{W1}^{top}");
    TCanvas *c70 = Plot_Me_2D("c70", hist_dphi_T1_T2_costhetaW1_W, "#Delta#phi_{T1,T2}^{W}", "cos #theta_{W1}^{W}");
    TCanvas *c71 = Plot_Me_2D("c71", hist_dphi_W_T1_costhetaT2_top, "#Delta#phi_{T1,W1}^{top}", "cos #theta_{T2}^{top}");
    TCanvas *c72 = Plot_Me_2D("c72", hist_dphi_W_T1_costhetaT2_W, "#Delta#phi_{T1,W1}^{W}", "cos #theta_{T2}^{W}");
    TCanvas *c73 = Plot_Me_2D("c73", hist_dphi_W_T1_costhetaW1_top, "#Delta#phi_{T1,W1}^{top}", "cos #theta_{W1}^{top}");
    TCanvas *c74 = Plot_Me_2D("c74", hist_dphi_W_T1_costhetaW1_W, "#Delta#phi_{T1,W1}^{W}", "cos #theta_{W1}^{W}");
    TCanvas *c75 = Plot_Me_2D("c75", hist_dphi_W_T1_costhetaW2_top, "#Delta#phi_{T1,W1}^{top}", "cos #theta_{W2}^{top}");
    TCanvas *c76 = Plot_Me_2D("c76", hist_dphi_W_T1_costhetaW2_W, "#Delta#phi_{T1,W1}^{W}", "cos #theta_{W2}^{W}");

    TCanvas *c77 = Plot_Me_2D("c77", hist_Eb1_dphi_T1_T2_top, "E_{b1}^{top}/E_{b1}^{true}", "#Delta#phi_{T1,T2}^{top}");
    TCanvas *c78 = Plot_Me_2D("c78", hist_Eb1_dphi_T1_T2_W, "E_{b1}^{W}/E_{b1}^{true}", "#Delta#phi_{T1,T2}^{W}");
    TCanvas *c79 = Plot_Me_2D("c79", hist_El1_dphi_T1_T2_top, "E_{l1}^{top}/E_{l1}^{true}", "#Delta#phi_{T1,T2}^{top}");
    TCanvas *c80 = Plot_Me_2D("c80", hist_El1_dphi_T1_T2_W, "E_{l1}^{W}/E_{l1}^{true}", "#Delta#phi_{T1,T2}^{W}");
    TCanvas *c81 = Plot_Me_2D("c81", hist_MTT_dphi_T1_T2_top, "M_{tt}^{top}/M_{tt}^{true}", "#Delta#phi_{T1,T2}^{top}");
    TCanvas *c82 = Plot_Me_2D("c82", hist_MTT_dphi_T1_T2_W, "M_{tt}^{W}/M_{tt}^{true}", "#Delta#phi_{T1,T2}^{W}");
    TCanvas *c83 = Plot_Me_2D("c83", hist_Eb1_dphi_W_T1_top, "E_{b1}^{top}/E_{b1}^{true}", "#Delta#phi_{T1,W1}^{top}");
    TCanvas *c84 = Plot_Me_2D("c84", hist_Eb1_dphi_W_T1_W, "E_{b1}^{W}/E_{b1}^{true}", "#Delta#phi_{T1,W1}^{W}");
    TCanvas *c85 = Plot_Me_2D("c85", hist_El1_dphi_W_T1_top, "E_{l1}^{top}/E_{l1}^{true}", "#Delta#phi_{T1,W1}^{top}");
    TCanvas *c86 = Plot_Me_2D("c86", hist_El1_dphi_W_T1_W, "E_{l1}^{W}/E_{l1}^{true}", "#Delta#phi_{T1,W1}^{W}");
    TCanvas *c87 = Plot_Me_2D("c87", hist_MTT_dphi_W_T1_top, "M_{tt}^{top}/M_{tt}^{true}", "#Delta#phi_{T1,W1}^{top}");
    TCanvas *c88 = Plot_Me_2D("c88", hist_MTT_dphi_W_T1_W, "M_{tt}^{W}/M_{tt}^{true}", "#Delta#phi_{T1,W1}^{W}");
    TCanvas *c89 = Plot_Me_2D("c89", hist_Eb1_dphi_W_T2_top, "E_{b1}^{top}/E_{b1}^{true}", "#Delta#phi_{T2,W2}^{top}");
    TCanvas *c90 = Plot_Me_2D("c90", hist_Eb1_dphi_W_T2_W, "E_{b1}^{W}/E_{b1}^{true}", "#Delta#phi_{T2,W2}^{W}");
    TCanvas *c91 = Plot_Me_2D("c91", hist_El1_dphi_W_T2_top, "E_{l1}^{top}/E_{l1}^{true}", "#Delta#phi_{T2,W2}^{top}");
    TCanvas *c92 = Plot_Me_2D("c92", hist_El1_dphi_W_T2_W, "E_{l1}^{W}/E_{l1}^{true}", "#Delta#phi_{T2,W2}^{W}");
    TCanvas *c93 = Plot_Me_2D("c93", hist_El1Eb1_dphi_T1_T2_top, "[E_{l1}/(E_{l1}+E_{b1})]^{top}", "#Delta#phi_{T1,T2}^{top}");
    TCanvas *c94 = Plot_Me_2D("c94", hist_El1Eb1_dphi_T1_T2_W, "[E_{l1}/(E_{l1}+E_{b1})]^{W}", "#Delta#phi_{T1,T2}^{W}");
    TCanvas *c95 = Plot_Me_2D("c95", hist_El1Eb2_dphi_T1_T2_top, "[E_{l1}/(E_{l1}+E_{b2})]^{top}", "#Delta#phi_{T1,T2}^{top}");
    TCanvas *c96 = Plot_Me_2D("c96", hist_El1Eb2_dphi_T1_T2_W , "[E_{l1}/(E_{l1}+E_{b2})]^{W}", "#Delta#phi_{T1,T2}^{W}");
    TCanvas *c97 = Plot_Me_2D("c97", hist_El1Eb1_costhetaT1_top, "[E_{l1}/(E_{l1}+E_{b1})]^{top}", "cos #theta_{T1}^{top}");
    TCanvas *c98 = Plot_Me_2D("c98", hist_El1Eb1_costhetaT1_W , "[E_{l1}/(E_{l1}+E_{b1})]^{W}", "cos #theta_{T1}^{W}");
    TCanvas *c99 = Plot_Me_2D("c99", hist_El1Eb2_costhetaT1_top, "[E_{l1}/(E_{l1}+E_{b2})]^{top}", "cos #theta_{T1}^{top}");
    TCanvas *c100 = Plot_Me_2D("c100", hist_El1Eb2_costhetaT1_W, "[E_{l1}/(E_{l1}+E_{b2})]^{W}", "cos #theta_{T1}^{W}");
    TCanvas *c101 = Plot_Me_2D("c101", hist_El1Eb1_costhetaW1_top, "[E_{l1}/(E_{l1}+E_{b1})]^{top}", "cos #theta_{W1}^{top}");
    TCanvas *c102 = Plot_Me_2D("c102", hist_El1Eb1_costhetaW1_W , "[E_{l1}/(E_{l1}+E_{b1})]^{W}", "cos #theta_{W1}^{W}");
    TCanvas *c103 = Plot_Me_2D("c103", hist_El1Eb2_costhetaW1_top, "[E_{l1}/(E_{l1}+E_{b2})]^{top}", "cos #theta_{W1}^{top}");
    TCanvas *c104 = Plot_Me_2D("c104", hist_El1Eb2_costhetaW1_W , "[E_{l1}/(E_{l1}+E_{b2})]^{W}", "cos #theta_{W1}^{W}");
    TCanvas *c105 = Plot_Me_2D("c105", hist_El1Eb1_dphi_W_T1_top , "[E_{l1}/(E_{l1}+E_{b1})]^{top}", "#Delta#phi_{T1,W1}^{top}");
    TCanvas *c106 = Plot_Me_2D("c106", hist_El1Eb1_dphi_W_T1_W , "[E_{l1}/(E_{l1}+E_{b1})]^{W}", "#Delta#phi_{T1,W1}^{W}");
    TCanvas *c107 = Plot_Me_2D("c107", hist_El1Eb2_dphi_W_T1_top, "[E_{l1}/(E_{l1}+E_{b2})]^{top}", "#Delta#phi_{T1,W1}^{top}");
    TCanvas *c108 = Plot_Me_2D("c108", hist_El1Eb2_dphi_W_T1_W , "[E_{l1}/(E_{l1}+E_{b2})]^{W}", "#Delta#phi_{T1,W1}^{W}");
    TCanvas *c109 = Plot_Me_2D("c109", hist_El1Eb1_costhetaT2_top, "[E_{l1}/(E_{l1}+E_{b1})]^{top}", "cos #theta_{T2}^{top}");
    TCanvas *c110 = Plot_Me_2D("c110", hist_El1Eb1_costhetaT2_W , "[E_{l1}/(E_{l1}+E_{b1})]^{W}", "cos #theta_{T2}^{W}");
    TCanvas *c111 = Plot_Me_2D("c111", hist_El1Eb2_costhetaT2_top, "[E_{l1}/(E_{l1}+E_{b2})]^{top}", "cos #theta_{T2}^{top}");
    TCanvas *c112 = Plot_Me_2D("c112", hist_El1Eb2_costhetaT2_W , "[E_{l1}/(E_{l1}+E_{b2})]^{W}", "cos #theta_{T2}^{W}");
    TCanvas *c113 = Plot_Me_2D("c113", hist_El1Eb1_costhetaW2_top, "[E_{l1}/(E_{l1}+E_{b1})]^{top}", "cos #theta_{W2}^{top}");
    TCanvas *c114 = Plot_Me_2D("c114", hist_El1Eb1_costhetaW2_W , "[E_{l1}/(E_{l1}+E_{b1})]^{W}", "cos #theta_{W2}^{W}");
    TCanvas *c115 = Plot_Me_2D("c115", hist_El1Eb2_costhetaW2_top, "[E_{l1}/(E_{l1}+E_{b2})]^{top}", "cos #theta_{W2}^{top}");
    TCanvas *c116 = Plot_Me_2D("c116", hist_El1Eb2_costhetaW2_W , "[E_{l1}/(E_{l1}+E_{b2})]^{W}", "cos #theta_{W2}^{W}");
    TCanvas *c117 = Plot_Me_2D("c117", hist_El1Eb1_dphi_W_T2_top , "[E_{l1}/(E_{l1}+E_{b1})]^{top}", "#Delta#phi_{T2,W2}^{top}");
    TCanvas *c118 = Plot_Me_2D("c118", hist_El1Eb1_dphi_W_T2_W , "[E_{l1}/(E_{l1}+E_{b1})]^{W}", "#Delta#phi_{T2,W2}^{W}");
    TCanvas *c119 = Plot_Me_2D("c119", hist_El1Eb2_dphi_W_T2_top, "[E_{l1}/(E_{l1}+E_{b2})]^{top}", "#Delta#phi_{T2,W2}^{top}");
    TCanvas *c120 = Plot_Me_2D("c120", hist_El1Eb2_dphi_W_T2_W, "[E_{l1}/(E_{l1}+E_{b2})]^{W}", "#Delta#phi_{T2,W2}^{W}");

    TCanvas *c121 = Plot_Me_2D("c121", hist_Eb1_costhetaT2_top, "E_{b1}^{top}/E_{b1}^{true}", "cos #theta_{T2}^{top}");
    TCanvas *c122 = Plot_Me_2D("c122", hist_Eb1_costhetaT2_W, "E_{b1}^{W}/E_{b1}^{true}", "cos #theta_{T2}^{W}");
    TCanvas *c123 = Plot_Me_2D("c123", hist_Eb1_costhetaW1_top, "E_{b1}^{top}/E_{b1}^{true}", "cos #theta_{W1}^{top}");
    TCanvas *c124 = Plot_Me_2D("c124", hist_Eb1_costhetaW1_W, "E_{b1}^{W}/E_{b1}^{true}", "cos #theta_{W1}^{W}");
    TCanvas *c125 = Plot_Me_2D("c125", hist_Eb1_costhetaW2_top, "E_{b1}^{top}/E_{b1}^{true}", "cos #theta_{W2}^{top}");
    TCanvas *c126 = Plot_Me_2D("c126", hist_Eb1_costhetaW2_W, "E_{b1}^{W}/E_{b1}^{true}", "cos #theta_{W2}^{W}");
    TCanvas *c127 = Plot_Me_2D("c127", hist_El1_costhetaT2_top, "E_{l1}^{top}/E_{l1}^{true}", "cos #theta_{T2}^{top}");
    TCanvas *c128 = Plot_Me_2D("c128", hist_El1_costhetaT2_W, "E_{l1}^{W}/E_{l1}^{true}", "cos #theta_{T2}^{W}");
    TCanvas *c129 = Plot_Me_2D("c129", hist_El1_costhetaW2_top, "E_{l1}^{top}/E_{l1}^{true}", "cos #theta_{W2}^{top}");
    TCanvas *c130 = Plot_Me_2D("c130", hist_El1_costhetaW2_W, "E_{l1}^{W}/E_{l1}^{true}", "cos #theta_{W2}^{W}");

    TCanvas *c131 = Plot_Me_2D("c131", hist_costhetaT1_costhetaW1_top, "cos #theta_{T1}^{top}", "cos #theta_{W1}^{top}");
    TCanvas *c132 = Plot_Me_2D("c132", hist_costhetaT1_costhetaW2_top , "cos #theta_{T1}^{top}", "cos #theta_{W2}^{top}");
    TCanvas *c133 = Plot_Me_2D("c133", hist_costhetaT1_costhetaW1_W, "cos #theta_{T1}^{W}", "cos #theta_{W1}^{W}");
    TCanvas *c134 = Plot_Me_2D("c134", hist_costhetaT1_costhetaW2_W, "cos #theta_{T1}^{W}", "cos #theta_{W2}^{W}");
    TCanvas *c135 = Plot_Me_2D("c135", hist_costhetaT1_costhetaW1, "cos #theta_{T1}^{true}", "cos #theta_{W1}^{true}");
    TCanvas *c136 = Plot_Me_2D("c136", hist_costhetaT1_costhetaW2, "cos #theta_{T1}^{true}", "cos #theta_{W2}^{true}");
    TCanvas *c137 = Plot_Me_2D("c137", hist_dphi_T1_W1_costhetaW1, "#Delta#phi_{T1,W1}^{true}", "cos #theta_{W1}^{true}");
    TCanvas *c138 = Plot_Me_2D("c138", hist_dphi_T1_W1_costhetaW2, "#Delta#phi_{T1,W1}^{true}", "cos #theta_{W2}^{true}");
    TCanvas *c139 = Plot_Me_2D("c139", hist_dphi_T1_W1_costhetaT1, "#Delta#phi_{T1,W1}^{true}", "cos #theta_{T1}^{true}");
    TCanvas *c140 = Plot_Me_2D("c140", hist_dphi_T1_W1_costhetaT2, "#Delta#phi_{T1,W1}^{true}", "cos #theta_{T2}^{true}");
}

void BoostToLabFrame(double gamma_eff){
  double M1 = P[0].M();
  double M2 = P[1].M();
  double shat = 4.*M1*M2*gamma_eff*gamma_eff + (M1-M2)*(M1-M2);
  double gamma_1 = (sqrt(shat)/(2.*M1))*(1.+ (M1*M1-M2*M2)/shat);
  double gamma_2 = (sqrt(shat)/(2.*M2))*(1.+ (M2*M2-M1*M1)/shat);
  double beta_1 = sqrt(1. - 1./(gamma_1*gamma_1));
  double beta_2 = sqrt(1. - 1./(gamma_2*gamma_2));
  double c_S = 1.-2.*gRandom->Rndm();
  double s_S = sqrt(1.-c_S*c_S);
  double phi_S = TMath::Pi()*2.*gRandom->Rndm();
  TVector3 vBETA1, vBETA2;
  vBETA1.SetXYZ(beta_1*cos(phi_S)*s_S,beta_1*sin(phi_S)*s_S,beta_1*c_S);
  vBETA2.SetXYZ(-beta_2*cos(phi_S)*s_S,-beta_2*sin(phi_S)*s_S,-beta_2*c_S);
	
  BoostHem(0, vBETA1);
  BoostHem(1, vBETA2);
	
  return;
	
}

void PtBoost(double PT){
  double phi_pt = TMath::Pi()*2.*gRandom->Rndm();
  TVector3 B_pt;
  B_pt.SetPtEtaPhi(PT/sqrt(PT*PT+1.),0.0,phi_pt);
			 
  BoostHem(0,B_pt);
  BoostHem(1,B_pt);
}

void BoostHem(int iHem, TVector3 vBETA){
  P[iHem].Boost(vBETA);
  vMiss[iHem].Boost(vBETA);
  C_1[iHem].Boost(vBETA);
  C_2[iHem].Boost(vBETA);
  C_1_1[iHem].Boost(vBETA);
  C_1_2[iHem].Boost(vBETA);
  C_1_1_1[iHem].Boost(vBETA);
  C_1_1_2[iHem].Boost(vBETA);
	
}
void GetLepTopHem(int iHem, double Mtop, double MW, double Mb, double Mlep, double Mnu){
  //last decay
  double PL = sqrt((MW*MW-(Mnu-Mlep)*(Mnu-Mlep))*(MW*MW-(Mnu+Mlep)*(Mnu+Mlep)))/(2.*MW);
  double EL = (MW*MW-Mnu*Mnu+Mlep*Mlep)/(2.*MW);
  double Enu = (MW*MW+Mnu*Mnu-Mlep*Mlep)/(2.*MW);
	
  //second to last
  double Pb = sqrt((Mtop*Mtop-(Mb-MW)*(Mb-MW))*(Mtop*Mtop-(Mb+MW)*(Mb+MW)))/(2.*Mtop);
  double Eb = (Mtop*Mtop-MW*MW+Mb*Mb)/(2.*Mtop);
  double EW = (Mtop*Mtop+MW*MW-Mb*Mb)/(2.*Mtop);
	
  //first, we do leptons/neutrinos in W frames:
  double c1_L = 1.-2.*gRandom->Rndm();
  double s1_L = sqrt(1.-c1_L*c1_L);
  double phi1_L = TMath::Pi()*2.*gRandom->Rndm();
	
  //lepton
  C_1_1[iHem].SetPxPyPzE(PL*cos(phi1_L)*s1_L,PL*sin(phi1_L)*s1_L,PL*c1_L,EL);
  //neutrino
  C_1_2[iHem].SetPxPyPzE(-PL*cos(phi1_L)*s1_L,-PL*sin(phi1_L)*s1_L,-PL*c1_L,Enu);
	
	
  //Now, we do b's and W's in T frames
  double c1_B = 1.-2.*gRandom->Rndm();
  double s1_B = sqrt(1.-c1_B*c1_B);
  double phi1_B = TMath::Pi()*2.*gRandom->Rndm();

  //b-quark
  C_2[iHem].SetPxPyPzE(Pb*cos(phi1_B)*s1_B,Pb*sin(phi1_B)*s1_B,Pb*c1_B,Eb);
  //W
  C_1[iHem].SetPxPyPzE(-Pb*cos(phi1_B)*s1_B,-Pb*sin(phi1_B)*s1_B,-Pb*c1_B,EW);
	
  //Now, put the leptons and nu's into T-frames
  TVector3 W1boost = C_1[iHem].BoostVector();
	
  C_1_1[iHem].Boost(W1boost);
  C_1_2[iHem].Boost(W1boost);
	
  //system of weakly interacting particles from this hemisphere
  vMiss[iHem] = C_1_2[iHem];
	
  //total system of particles = parent
  P[iHem] = C_1[iHem]+C_2[iHem];
	
  return;
	
}


void setstyle() {
  
  // For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(300); //Height of canvas
  gStyle->SetCanvasDefW(600); //Width of canvas
  gStyle->SetCanvasDefX(0);   //POsition on screen
  gStyle->SetCanvasDefY(0);

  // For the Pad:
  gStyle->SetPadBorderMode(0);
  // gStyle->SetPadBorderSize(Width_t size = 1);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  // For the frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);

  // set the paper & margin sizes
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.065);
  gStyle->SetPadRightMargin(0.22);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.17);

  // use large Times-Roman fonts
  gStyle->SetTitleFont(132,"xyz");  // set the all 3 axes title font
  gStyle->SetTitleFont(132," ");    // set the pad title font
  gStyle->SetTitleSize(0.06,"xyz"); // set the 3 axes title size
  gStyle->SetTitleSize(0.06," ");   // set the pad title size
  gStyle->SetLabelFont(132,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelColor(1,"xyz");
  gStyle->SetTextFont(132);
  gStyle->SetTextSize(0.08);
  gStyle->SetStatFont(132);

  // use bold lines and markers
  gStyle->SetMarkerStyle(8);
  gStyle->SetHistLineWidth(1.85);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  //..Get rid of X error bars
  gStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  const Int_t NRGBs = 5;
  const Int_t NCont = 25;
	
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };   
  Double_t red[NRGBs] =   { 1.0,   0.95,  0.95,  0.65,   0.15 };
  Double_t green[NRGBs] = { 1.0,  0.85, 0.7, 0.5,  0.3 };
  Double_t blue[NRGBs] =  { 0.95, 0.6 , 0.3,  0.45, 0.65 };
  double mr[NRGBs],mg[NRGBs],mb[NRGBs];
  for(int i = 0; i < NRGBs; i++){
    mr[i] = red[NRGBs-1-i];
    mg[i] = green[NRGBs-1-i];
    mb[i] = blue[NRGBs-1-i];
  }
  TColor::CreateGradientColorTable(NRGBs, stops, mr, mg, mb, NCont);
	
  gStyle->SetNumberContours(NCont);
	
  gStyle->cd();

}

double PDF_v(double x){
	
  double ret;
	
  ret = (1/x)*pow(x,eta_1)*pow(1-x,eta_2)*(1.+eps_u*sqrt(x)+g_u*x);
	
  return ret;
}
double PDF_S(double x){
	
  double ret;
	
  ret = (1/x)*pow(x,del_S)*pow(1-x,eta_S)*(1.+eps_S*sqrt(x)+g_S*x);
	
  return ret;
}
double PDF_g(double x){
	
  double ret;
	
  ret = (1/x)*(A_g*pow(x,del_g)*pow(1-x,eta_g)*(1.+eps_g*sqrt(x)+g_g*x)+A_g1*pow(x,del_g1)*pow(1-x,eta_g1));
	
  return ret;
}
double PDF_TOT(double x){
	
  double ret;
	
  if(type == 0){ // q qbar
    ret = (1/x)*(PDF_v(x)*PDF_S(C/x)+PDF_S(x)*PDF_v(C/x));
  } else { // gg
    ret = (1/x)*(PDF_g(x)*PDF_g(C/x));
  }
	
  return ret;
}

double Calc_dsigma_dgamma(double gamma){
	
  C = 4.*gamma*gamma*M*M/(rootS*rootS);
  TF1 *func = new TF1("func","PDF_TOT(x)",0,10);
	
  double ret;
	
  //ret = (1./pow(gamma,5.))*sqrt(1.-1./(gamma*gamma))*func->Integral(C,1.);
  ret = (1./gamma)*sqrt(1.-1./(gamma*gamma))*func->Integral(C,1.);
	
  return ret;
}

TCanvas* Plot_Me_2D(char *titlecan, TH2D* histo, char *titleX, char *titleY){
  TCanvas *c1 = new TCanvas(titlecan,titlecan,600,500);
  c1->Draw();
  c1->SetGridx();
  c1->SetGridy();
  histo->Scale(1./histo->Integral());
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetXaxis()->CenterTitle();
  histo->GetYaxis()->CenterTitle();
  histo->Draw("COLZ");
  char *sname = new char[200];
  sprintf(sname,"fig_%s.pdf",histo->GetTitle());
  //c1->Print(sname);	
  return c1;
	
}


