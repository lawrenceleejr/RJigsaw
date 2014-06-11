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

#include "TRJigsaw.h"

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

int main(){
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
	int N = 1000;


	////////////////////////////////////////////////////////////
	// Create and initialize RJigsaw Tool
	////////////////////////////////////////////////////////////

	Root::TRJigsaw* RJTool = new Root::TRJigsaw();
	RJTool->initialize("./RJigsawConfig/hemisphere1","./RJigsawConfig/hemisphere2");
	RJTool->resetHists();

	// Here are the variables I can use. The indicies are hemisphere (0,1,2), generation, hemibalancemode (0,1)
	// // angles
	// "cosTheta_%d_%d_%d"
	// "dPhi_%d_%d_%d" 
	// "dPhiVis_%d_%d_%d" 
	// // scale
	// "M_%d_%d_%d" 
	// "dPhiDecay_%d_%d_%d" 
	// "cosThetaDecay_%d_%d_%d" 
	// "cosThetaDecay_%d_%d_%d" 
	// "dPhiDecay_%d_%d_%d" 
	// "dPhiDecay_%d_%d_%d" 
	// "dPhiDecay_%d_%ßd_%d" 
	// "gamma_%d_%d_%d"

	////////////////////////////////////////////////////////////
	// Book histograms. Nothing more to do after this.
	// Tool takes in expression that matches key in observables
	// map
	////////////////////////////////////////////////////////////

	RJTool->bookHist1D("M_0_0_0",100.,0.,1000.);
	RJTool->bookHist1D("M_0_0_1",100.,0.,1000.);
	RJTool->bookHist1D("M_0_0_-1",100.,0.,1000.);

	RJTool->bookHist1D("M_1_0_0",100.,0.,1000.);
	RJTool->bookHist1D("M_1_0_1",100.,0.,1000.);
	RJTool->bookHist1D("M_1_0_-1",100.,0.,1000.);
	RJTool->bookHist1D("M_2_0_0",100.,0.,1000.);
	RJTool->bookHist1D("M_2_0_1",100.,0.,1000.);
	RJTool->bookHist1D("M_2_0_-1",100.,0.,1000.);

	RJTool->bookHist1D("M_1_1_0",100.,0.,1000.);
	RJTool->bookHist1D("M_1_1_1",100.,0.,1000.);
	RJTool->bookHist1D("M_1_1_-1",100.,0.,1000.);
	RJTool->bookHist1D("M_2_1_0",100.,0.,1000.);
	RJTool->bookHist1D("M_2_1_1",100.,0.,1000.);
	RJTool->bookHist1D("M_2_1_-1",100.,0.,1000.);


	RJTool->bookHist2D("M_0_0_-1","M_0_0_0", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_0_0_-1","M_0_0_1", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_0_0_0" ,"M_0_0_1", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 

	RJTool->bookHist2D("M_1_0_-1","M_1_0_0", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_1_0_-1","M_1_0_1", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_1_0_0" ,"M_1_0_1", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_2_0_-1","M_2_0_0", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_2_0_-1","M_2_0_1", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_2_0_0" ,"M_2_0_1", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 


	RJTool->bookHist2D("M_1_1_-1","M_1_1_0", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_1_1_-1","M_1_1_1", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_1_1_0" ,"M_1_1_1", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_2_1_-1","M_2_1_0", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_2_1_-1","M_2_1_1", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_2_1_0" ,"M_2_1_1", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 

	RJTool->bookHist2D("M_1_0_-1","M_2_0_-1", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_1_0_0","M_2_0_0", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_1_0_1","M_2_0_1", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 

	RJTool->bookHist2D("M_1_1_-1","M_2_1_-1", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_1_1_0" ,"M_2_1_0", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_1_1_1" ,"M_2_1_1", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 

	RJTool->bookHist2D("M_1_0_-1","M_1_1_-1", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_1_0_0" ,"M_1_1_0", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 
	RJTool->bookHist2D("M_1_0_1" ,"M_1_1_1", 100. , 0. , 1000. , 100. , 0. , 1000.   ); 

	// Replicating original code

	RJTool->bookHist2D("dPhi_0_0_-1","dPhi_0_0_0", 50. , 0. , TMath::Pi()   , 50. , 0. , TMath::Pi()   ); 
	RJTool->bookHist2D("dPhi_0_0_-1","dPhi_0_0_1", 50. , 0. , TMath::Pi()   , 50. , 0. , TMath::Pi()   ); 

	RJTool->bookHist2D("dPhiVis_0_0_-1","dPhiVis_0_0_0", 50. , 0. , TMath::Pi()   , 50. , 0. , TMath::Pi()   ); 
	RJTool->bookHist2D("dPhiVis_0_0_-1","dPhiVis_0_0_1", 50. , 0. , TMath::Pi()   , 50. , 0. , TMath::Pi()   ); 

	RJTool->bookHist2D("dPhi_0_0_-1","dPhiVis_0_0_-1", 50. , 0. , TMath::Pi()   , 50. , 0. , TMath::Pi()   ); 
	RJTool->bookHist2D("dPhi_0_0_0","dPhiVis_0_0_0", 50. , 0. , TMath::Pi()   , 50. , 0. , TMath::Pi()   ); 
	RJTool->bookHist2D("dPhi_0_0_1","dPhiVis_0_0_1", 50. , 0. , TMath::Pi()   , 50. , 0. , TMath::Pi()   ); 

	RJTool->bookHist2D("M_0_0_-1","dPhiVis_0_0_-1", 50. , 0. , 1000.  , 50. , 0. , TMath::Pi()   ); 
	RJTool->bookHist2D("M_0_0_0","dPhiVis_0_0_0", 50. , 0. , 1000.  , 50. , 0. , TMath::Pi()   ); 
	RJTool->bookHist2D("M_0_0_1","dPhiVis_0_0_1", 50. , 0. , 1000.  , 50. , 0. , TMath::Pi()   ); 

	RJTool->bookHist2D("gamma_0_0_-1","gamma_0_0_0", 50. , 0. , 10. , 50. , 0. , 10. ); 
	RJTool->bookHist2D("gamma_0_0_-1","gamma_0_0_1", 50. , 0. , 10. , 50. , 0. , 10. ); 






	////////////////////////////////////////////////////////////
	// Generate fake events taking flat ME's for all decay angles
	////////////////////////////////////////////////////////////


	//here, we set up the stuff for dynamic gamma
	TH1D *h_gamma = (TH1D*) new TH1D("newgamma","newgamma",500,1.0,6.0);
	if(b_gamma){
		cout << "generating gamma distribution" << endl;
		for(int ibin = 1; ibin <= 500; ibin++){
			double g = h_gamma->GetBinCenter(ibin);
			double entry = Calc_dsigma_dgamma(g);

			if(entry > 0.) h_gamma->SetBinContent(ibin, entry);
		}
		cout << "done" << endl;
	}

	for(int i = 0; i < N; i++){

		if(i%100==0){
			cout << i << endl;
		}

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


		////////////////////////////////////////////////////////////
		// The event has been generated. Let's hand it to the tool
		////////////////////////////////////////////////////////////

		// Reset tool from last event

		RJTool->newEvent();

		// Give b's to the tool

		RJTool->addVisParticle("b",C_2[0],1);
		RJTool->addVisParticle("b",C_2[1],2);

		// Give leptons to the tool

		RJTool->addVisParticle("l",C_1_1[0],1);
		RJTool->addVisParticle("l",C_1_1[1],2);


		// Give truth nus to the tool

		RJTool->addTruParticle("nu",C_1_2[0],1);
		RJTool->addTruParticle("nu",C_1_2[1],2);

		// Let's construct a fake MET vector

		TVector3 MET = (C_1_2[0]+C_1_2[1]).Vect();
		MET.SetZ(0.0);

		// And then tell the tool about it

		RJTool->addMET( MET );

		// Set the hemisphere balance mode to use the top symmetry and guess neutrino momenta

		RJTool->setHemisphereMode(0); //top symmetry
		RJTool->guessInvParticles();

		// Calculate everything and fill the observables map with all variables
		// this is where the magic happens

		RJTool->getObservables();

		// If you want access to the variables to do whatever with...

		std::map< TString, double > observables = RJTool->getObservablesMap();

		// std::cout << observables["M_0_0_0"] << std::endl;

		// Try again with other Hemisphere mode

		RJTool->setHemisphereMode(1); //W symmetry
		RJTool->guessInvParticles();
		RJTool->getObservables();

		// Try again using true neutrinos

		RJTool->setHemisphereMode(-1); //W symmetry
		RJTool->guessInvParticles();
		RJTool->getObservables();

		// Fill the hists with weight 1.

		RJTool->fillHists(1.);

	}

	// Write out all the histograms to a root file

	TFile* f = new TFile("./output.root", "RECREATE");
	RJTool->writeHists(f);

	return 0;

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
	//TF1 *func = new TF1("func","PDF_TOT(x)",0,10);
	TF1 *func = new TF1("func","1",0,10);

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


