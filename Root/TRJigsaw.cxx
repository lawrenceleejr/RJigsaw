/******************************************************************************
Name:        TRJigsaw

Author:      
Created:     

Description: 
******************************************************************************/

// Preprocessor magic for debugging
#define XXX std::cout << " I am here: " << __FILE__ << ":" << __LINE__ << std::endl;

// This class' header
#include "RJigsaw/TRJigsaw.h"

// include math
#include <math.h>

// ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TLeaf.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TAxis.h>
#include <TString.h>
#include <TRandom3.h>
#include <assert.h>

#include <sstream>


ClassImp(Root::TRJigsaw)


////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////

Root::TRJigsaw::TRJigsaw(const char* name) :
TNamed(name,"notitle")
{
	visParticles.clear();
	 //METVector;
	hemiBalanceMode = 0;
}

////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////

Root::TRJigsaw::~TRJigsaw()
{
}


////////////////////////////////////////////////////////////
// Initialization
////////////////////////////////////////////////////////////

void Root::TRJigsaw::initialize(std::string filename1, std::string filename2){

	// Reads in filenameN and stores particle listings into vectors

	std::cout << "###################################################" << std::endl;
	std::cout << "###################################################" << std::endl;
	std::cout << "##" << std::endl;
	std::cout << "## TRJigsaw::initialize()" << std::endl;
	std::cout << "##" << std::endl;
	std::cout << "## Initializing with:" << std::endl;
	std::cout << "## Hemisphere 1: " << filename1 << std::endl;
	std::cout << "## Hemisphere 2: " << filename2 << std::endl;
	std::cout << "###################################################" << std::endl;
	std::cout << "###################################################" << std::endl;

	std::cout <<  std::endl;
	std::cout <<  std::endl;

	std::cout << "###################################################" << std::endl;
	std::cout << "##" << std::endl;
	std::cout << "## Config below " << std::endl;
	std::cout << "##" << std::endl;


	hemisphereConfig[0].clear();
	hemisphereConfig[1].clear();
	hemisphereConfig[2].clear();

	std::string STRING;
	std::ifstream infile;

	infile.open(filename1.data());
	assert(!infile.fail() );
	while(!infile.eof())
	{
		getline(infile,STRING);
		std::stringstream ss(STRING);
		//std::cout << STRING << std::endl;
		std::istream_iterator<std::string> begin(ss);
		std::istream_iterator<std::string> end;
		std::vector<std::string> vstrings(begin,end);
		std::copy(vstrings.begin(), vstrings.end(), std::ostream_iterator<std::string>(std::cout, " ") );
		std::vector< TString > tStringVector( vstrings.begin(), vstrings.end() ) ;
		if( !tStringVector.empty()) hemisphereConfig[1].push_back(  tStringVector );
		std::cout << std::endl;
		std::cout << std::endl;
	}
	infile.close();

	infile.open(filename2.data());
	assert(!infile.fail() );
	while(!infile.eof())
	{
		getline(infile,STRING);
		std::stringstream ss(STRING);		
		//std::cout << STRING << std::endl;
		std::istream_iterator<std::string> begin(ss);
		std::istream_iterator<std::string> end;
		std::vector<std::string> vstrings(begin,end);
		std::copy(vstrings.begin(), vstrings.end(), std::ostream_iterator<std::string>(std::cout, " ") );
		std::vector< TString > tStringVector( vstrings.begin(), vstrings.end() ) ;
		if(!tStringVector.empty()) hemisphereConfig[2].push_back(  tStringVector );
		std::cout << std::endl;
		std::cout << std::endl;
	}
	infile.close();

	// This just makes things easier for the first vertex
	std::vector<TString> dummyVector;
	dummyVector.push_back("system");

	hemisphereConfig[0].push_back(dummyVector);


	std::cout << "##" << std::endl;
	std::cout << "###################################################" << std::endl;


	return;


}

void Root::TRJigsaw::addVisParticle(
	TString particleType,
	TLorentzVector particleMomentum,
	int hemisphere = 0
	){

	particleClass particle;
	particle.particleType = particleType;
	particle.particleMomentum = particleMomentum;
	particle.particleMomentumForBoosting = particleMomentum;
	particle.hemisphere = hemisphere;
	visParticles.push_back( particle );

	return;

}

void Root::TRJigsaw::addTruParticle(
	TString particleType,
	TLorentzVector particleMomentum,
	int hemisphere = 0
	){

	particleClass particle;
	particle.particleType = particleType;
	particle.particleMomentum = particleMomentum;
	particle.particleMomentumForBoosting = particleMomentum;
	particle.hemisphere = hemisphere;

	truParticles.push_back( particle );

	return;

}


void Root::TRJigsaw::guessInvParticles(){

	invParticles.clear();

	TLorentzVector Inv1, Inv2;

	//We separate the objects into decay hemispheres

	TLorentzVector H1 = sumParticles(false, 0, -1, 0, 1);
	TLorentzVector H2 = sumParticles(false, 0, -1, 0, 2);

	//Now, we must perform longitudinal boost from lab frame to CMz frame

	TVector3 BL = H1.Vect()+H2.Vect();
	BL.SetX(0.0);
	BL.SetY(0.0);
	BL = (1./(H1.E()+H2.E()))*BL;

	H1.Boost(-BL);
	H2.Boost(-BL);

	if(hemiBalanceMode==-1){ // get inv particles from truth

		// Set InvN to truth nu momenta
		Inv1 = findTruParticleMomentum("nu",1);
		Inv2 = findTruParticleMomentum("nu",2);
		
	} else if(hemiBalanceMode==0){ // for ttbar, this is a top symmetry mode

		//Now, we need to 'guess' the invariant mass of the weakly interacting system

		double HM2min = std::min(H1.M2(), H2.M2() );

		double Minv2 = (H1+H2).M2() - 4.*HM2min;

		TVector3 PT = *METVector+(H1+H2).Vect();
		double Einv2 = METVector->Mag2() + Minv2;
		PT.SetZ(0.0);
		TVector3 BT = (1./( (H1+H2).E() + sqrt(Einv2) ))*PT;

		H1.Boost(-BT);
		H2.Boost(-BT);

		//Now, in CM approx frame

		double c = ( H1.E()+H2.E() + sqrt( (H1.E()+H2.E())*(H1.E()+H2.E()) -(H1+H2).M2() + Minv2 ))/(H1.E()+H2.E())/2.;

		double Enu1 = (c-1.)*H1.E() + c*H2.E();
		double Enu2 = c*H1.E() + (c-1.)*H2.E();
		TVector3 pnu1 = (c-1.)*H1.Vect() - c*H2.Vect();
		TVector3 pnu2 = (c-1.)*H2.Vect() - c*H1.Vect();

		//define neutrino 4-vectors

		Inv1.SetPxPyPzE(pnu1.X(),pnu1.Y(),pnu1.Z(),Enu1);
		Inv2.SetPxPyPzE(pnu2.X(),pnu2.Y(),pnu2.Z(),Enu2);

		//now, go back to lab frame

		Inv1.Boost(BT);
		Inv2.Boost(BT);
		Inv1.Boost(BL);
		Inv2.Boost(BL);

	} else if(hemiBalanceMode==1){ // for ttbar, this is a W symmetry mode

		TLorentzVector L1, L2;
		L1 = sumParticles(false, 0, -1, 1, 1);
		L2 = sumParticles(false, 0, -1, 1, 2);

		//Now, we need to 'guess' the invariant mass of the weakly interacting system

		double Minv2 = (L1+L2).M2();
		TVector3 PT = *METVector+(L1+L2).Vect();
		double Einv2 = METVector->Mag2() + Minv2;
		PT.SetZ(0.0);
		TVector3 BT = (1./( (L1+L2).E() + sqrt(Einv2) ))*PT;
		L1.Boost(-BT);
		L2.Boost(-BT);

		//Now, in CM approx frame _of WW system_

		double c = ( L1.E()+L2.E() + sqrt( (L1.E()+L2.E())*(L1.E()+L2.E()) -(L1+L2).M2() + Minv2 ))/(L1.E()+L2.E())/2.;

		double Enu1 = (c-1.)*L1.E() + c*L2.E();
		double Enu2 = c*L1.E() + (c-1.)*L2.E();
		TVector3 pnu1 = (c-1.)*L1.Vect() - c*L2.Vect();
		TVector3 pnu2 = (c-1.)*L2.Vect() - c*L1.Vect();

		//define neutrino 4-vectors

		Inv1.SetPxPyPzE(pnu1.X(),pnu1.Y(),pnu1.Z(),Enu1);
		Inv2.SetPxPyPzE(pnu2.X(),pnu2.Y(),pnu2.Z(),Enu2);

		//now, go back to lab frame

		Inv1.Boost(BT);
		Inv2.Boost(BT);
		Inv1.Boost(BL);
		Inv2.Boost(BL);

	} else { 
		std::cout << "Unrecognized hemiBalanceMode Value" << std::endl;
		assert(0); 
	}

	particleClass invParticle1;
	invParticle1.particleType = TString("inv");
	invParticle1.particleMomentum = Inv1;
	invParticle1.particleMomentumForBoosting = Inv1;
	invParticle1.hemisphere = 1;

	particleClass invParticle2;
	invParticle2.particleType = TString("inv");
	invParticle2.particleMomentum = Inv2;
	invParticle2.particleMomentumForBoosting = Inv2;
	invParticle2.hemisphere = 2;

	invParticles.push_back(invParticle1);
	invParticles.push_back(invParticle2);

	return;

}


std::pair<TLorentzVector,TLorentzVector> Root::TRJigsaw::calcHemispheres(){

	// ATH_MSG_DEBUG("size of m_selectedInputJetContainer " << nJets );

	/////////// Starting Hemisphere Reconstruction ////////////////////////
	/////////// Adapted from Chris Rogan //////////////////////////////////

	double M_min = -1.;
	int j_count;

	int itemp = 0;
	int count = 0;

	double M_temp=0;

	TLorentzVector h1,h2;
	TLorentzVector h1temp,h2temp;

	// Brute force loop over 2^N-1 jet combinations //////////////////////
	// Not great, we know. But it does the job for now... ////////////////

	int N_comb = pow(2.,visParticles.size() );

	for(int i = 1; i < N_comb-1; i++){

		itemp = i;
		j_count = N_comb/2.;
		count = 0;
		h1temp *= 0; h2temp *= 0;
		while(j_count > 0){
		    if(itemp/j_count == 1){
		        h1temp += visParticles.at(count).particleMomentum;
		    } else {
		        h2temp += visParticles.at(count).particleMomentum;
		    }
		    itemp -= j_count*(itemp/j_count);
		    j_count /= 2.;
		    count++;
		}

		if(h1temp.M2()==0. || h2temp.M2()==0.) continue;

		M_temp = h1temp.M2() + h2temp.M2();

		// Trying to minimize combined hemisphere mass //////////

		if(M_temp < M_min || M_min < 0){
		    M_min = M_temp;
		    h1 = h1temp;
		    h2 = h2temp;
		}
	}

	if(h2.Pt() > h1.Pt() ){
	  TLorentzVector temp = h1;
	  h1 = h2;
	  h2 = temp;
	}

	return std::pair<TLorentzVector,TLorentzVector>(h1,h2);

}

void Root::TRJigsaw::getObservables(TLorentzVector h1, TLorentzVector h2){

	TVector3 vBETA_z = (1./(h1.E()+h2.E()))*(h1+h2).Vect();
	vBETA_z.SetX(0.0);
	vBETA_z.SetY(0.0);

	//transformation from lab frame to approximate rest frame along beam-axis
	h1.Boost(-vBETA_z);
	h2.Boost(-vBETA_z);

	TVector3 pT_CM = (h1+h2).Vect() + *METVector;
	pT_CM.SetZ(0.0); //should be redundant...

	float m_Minv2 = (h1+h2).M2() - 4.*h1.M()*h2.M();
	float m_Einv = sqrt((*METVector).Mag2()+m_Minv2);

	//////////////////////
	// definition of m_shatR
	//////////////////////
	float m_shatR = sqrt( ((h1+h2).E()+m_Einv)*((h1+h2).E()+m_Einv) - pT_CM.Mag2() );

	TVector3 vBeta_R = (1./sqrt(pT_CM.Mag2() + m_shatR*m_shatR))*pT_CM;
    float m_gammainv_R = sqrt(1.-vBeta_R.Mag2());

	//transformation from lab frame to R frame
	h1.Boost(-vBeta_R);
	h2.Boost(-vBeta_R);

	/////////////
	//
	// R-frame
	//
	/////////////


	TLorentzVector h1h2 = h1+h2;

	float m_dphi_Beta_R = (h1h2.Vect()).DeltaPhi(vBeta_R);
	float m_dphi_leg1_leg2 = h1.Vect().DeltaPhi(h2.Vect());
	float m_costheta_R =  fabs(h1h2.Vect().Dot(vBeta_R)/(h1h2.Vect().Mag()*vBeta_R.Mag()));

	TVector3 vBeta_Rp1 = (1./(h1.E()+h2.E()))*(h1.Vect() - h2.Vect());

	////////////////////////
	// definition of m_gaminvR
	////////////////////////

	float m_gammainv_Rp1 = sqrt(1.-vBeta_Rp1.Mag2());
	float m_dphi_Beta_Rp1_Beta_R = vBeta_Rp1.DeltaPhi(vBeta_R);


	h1.Boost(-vBeta_Rp1);
	h2.Boost(vBeta_Rp1);

	float m_mdelta_R = 2.*h1.E();
	float m_costheta_Rp1 = fabs(h1.Vect().Dot(vBeta_Rp1)/(h1.Vect().Mag()*vBeta_Rp1.Mag()));

	observables[ TString::Format("sHatR") ]  = m_shatR;
	observables[ TString::Format("gammainv_R") ]  = m_gammainv_R;
	observables[ TString::Format("dphi_Beta_R") ]  = m_dphi_Beta_R;
	observables[ TString::Format("dphi_leg1_leg2") ]  = m_dphi_leg1_leg2;
	observables[ TString::Format("costheta_R") ]  = m_costheta_R;
	observables[ TString::Format("gammainv_Rp1") ]  = m_gammainv_Rp1;
	observables[ TString::Format("dphi_Beta_Rp1_Beta_R") ]  = m_dphi_Beta_Rp1_Beta_R;
	observables[ TString::Format("mdelta_R") ]  = m_mdelta_R;
	observables[ TString::Format("costheta_Rp1") ]  = m_costheta_Rp1;

}


void Root::TRJigsaw::getObservables(){

	// Ok I have to clean this stuff up now

	// Let's declare some stuff ahead of the loop

	TVector3 tmpBoostVector;
	TLorentzVector leg1, leg1Vis;
	TLorentzVector leg2, leg2Vis;
	TLorentzVector leg1_1, leg1_2;
	TLorentzVector leg2_1, leg2_2;
	TVector3 decayPlaneVector1;
	TVector3 decayPlaneVector2;

	///////////////////////////////////////////////////////////
	// Looping over hemispheres and generations for each
	// hemisphere. 
	// iHemisphere = 0: H1+H2 Frame 
	///////////////////////////////////////////////////////////

	for( int iHemisphere = 0; iHemisphere < 3; iHemisphere++){

			for( unsigned int iGeneration = 0; iGeneration < hemisphereConfig[iHemisphere].size(); iGeneration++ ){ // iGeneration is index in config file line

				// Really can't boost the first vertex more than once

				if(iHemisphere==0) assert(hemisphereConfig[0].size()==1) ;

				////////////////////////////////////////////////////////////
				// Let's start by boosting into this frame
				////////////////////////////////////////////////////////////

				// Calculate total momenta for children at vertex
				// Also momenta of the N+1 step for decay angles

				if(iHemisphere==0){
					leg1    = sumParticles(true, 0, -1, 0, 1); 
					leg1Vis = sumParticles(false, 0, -1, 0, 1);
					leg2    = sumParticles(true, 0, -1, 0, 2);
					leg2Vis = sumParticles(false, 0, -1, 0, 2);
					leg1_1  = sumParticles(true, 0, -1, 1, 1); // Child of legN
					leg2_1  = sumParticles(true, 0, -1, 1, 2);

				} else {
					leg1    = sumParticles(true, 0, 1, iGeneration, iHemisphere);  //b - leg 1 assumed to be stable - ORDER CONFIG FILE APPROPRIATELY
					leg1Vis = sumParticles(false, 0, 1, iGeneration, iHemisphere);
					leg2    = sumParticles(true, 1, -1, iGeneration, iHemisphere); 
					leg2Vis = sumParticles(false, 1, -1, iGeneration, iHemisphere);
					leg1_1  = leg1; // Assumed to be stable
					leg2_1  = sumParticles(true, 0, 1, iGeneration+1, iHemisphere); 
				}

				// Assuming 4mom conservation to get the other child of legN

				leg1_2 = leg1 - leg1_1;
				leg2_2 = leg2 - leg2_1; // eg p_nu = p_W - p_lep

				// Calculate the boost of this whole system

				tmpBoostVector = (leg1Vis+leg2Vis).BoostVector();
				tmpBoostVector.SetX(0.0);
				tmpBoostVector.SetY(0.0);

				// Do the boosting for all relevant particles

				boostParticles(-tmpBoostVector, true, 0, -1, iGeneration, iHemisphere); 

				// and momentum sums

				leg1.Boost(-tmpBoostVector);
				leg1Vis.Boost(-tmpBoostVector);
				leg2.Boost(-tmpBoostVector);
				leg2Vis.Boost(-tmpBoostVector);
				leg1_1.Boost(-tmpBoostVector);
				leg1_2.Boost(-tmpBoostVector);
				leg2_1.Boost(-tmpBoostVector);
				leg2_2.Boost(-tmpBoostVector);


				TVector3 pT_CM = (leg1Vis+leg2Vis).Vect() + *METVector;
				pT_CM.SetZ(0.0);

				

				TLorentzVector leg1leg2 = leg1Vis+leg2Vis;

			    float SHATR = sqrt( 2.*(leg1leg2.E()*leg1leg2.E() - leg1leg2.Vect().Dot(pT_CM) + leg1leg2.E()*sqrt( leg1leg2.E()*leg1leg2.E() + pT_CM.Mag2() - 2.*leg1leg2.Vect().Dot(pT_CM) )));
			    float Eleg1leg2 = leg1leg2.E();

		        TVector3 vBeta_R = (1./sqrt(pT_CM.Mag2() + SHATR*SHATR))*pT_CM;
	            float gamma_R = 1./sqrt(1.-vBeta_R.Mag2());


				boostParticles(-vBeta_R, true, 0, -1, iGeneration, iHemisphere); 

				// and momentum sums

				leg1.Boost(-vBeta_R);
				leg1Vis.Boost(-vBeta_R);
				leg2.Boost(-vBeta_R);
				leg2Vis.Boost(-vBeta_R);
				leg1_1.Boost(-vBeta_R);
				leg1_2.Boost(-vBeta_R);
				leg2_1.Boost(-vBeta_R);
				leg2_2.Boost(-vBeta_R);


				float dphi_Beta_R = (leg1leg2.Vect()).DeltaPhi(vBeta_R);
				float dphi_leg1_leg2 = leg1Vis.Vect().DeltaPhi(leg2Vis.Vect());

				float costhetaR =  fabs(leg1leg2.Vect().Dot(vBeta_R)/(leg1leg2.Vect().Mag()*vBeta_R.Mag()));

			    TVector3 vBeta_Rp1 = (1./(leg1Vis.E()+leg2Vis.E()))*(leg1Vis.Vect() - leg2Vis.Vect());
			    float gamma_Rp1 = 1./sqrt(1.-vBeta_Rp1.Mag2());

				float dphi_Beta_Rp1_Beta_R = vBeta_Rp1.DeltaPhi(vBeta_R);

				leg1Vis.Boost(-vBeta_Rp1);
				leg2Vis.Boost(vBeta_Rp1);

				float MDeltaR = 2.*leg1Vis.E();

				float Eleg1 = leg1Vis.E();
				float Eleg2 = leg2Vis.E();

				float costhetaRp1 = fabs(leg1Vis.Vect().Dot(vBeta_Rp1)/(leg1Vis.Vect().Mag()*vBeta_Rp1.Mag()));


				observables[ TString::Format("sHatR_%d_%d_%d", iHemisphere, iGeneration, hemiBalanceMode) ]  = SHATR;
				observables[ TString::Format("E12_%d_%d_%d", iHemisphere, iGeneration, hemiBalanceMode) ]  = Eleg1leg2;
				observables[ TString::Format("gamma_R_%d_%d_%d", iHemisphere, iGeneration, hemiBalanceMode) ]  = gamma_R;
				observables[ TString::Format("dphi_Beta_R_%d_%d_%d", iHemisphere, iGeneration, hemiBalanceMode) ]  = dphi_Beta_R;
				observables[ TString::Format("dphi_leg1_leg2_%d_%d_%d", iHemisphere, iGeneration, hemiBalanceMode) ]  = dphi_leg1_leg2;
				observables[ TString::Format("gamma_Rp1_%d_%d_%d", iHemisphere, iGeneration, hemiBalanceMode) ]  = gamma_Rp1;
				observables[ TString::Format("dphi_Beta_Rp1_Beta_R_%d_%d_%d", iHemisphere, iGeneration, hemiBalanceMode) ]  = dphi_Beta_Rp1_Beta_R;
				observables[ TString::Format("MDeltaR_%d_%d_%d", iHemisphere, iGeneration, hemiBalanceMode) ]  = MDeltaR;
				observables[ TString::Format("Eleg1_%d_%d_%d", iHemisphere, iGeneration, hemiBalanceMode) ]  = Eleg1;
				observables[ TString::Format("Eleg2_%d_%d_%d", iHemisphere, iGeneration, hemiBalanceMode) ]  = Eleg2;
				observables[ TString::Format("costhetaRp1_%d_%d_%d", iHemisphere, iGeneration, hemiBalanceMode) ]  = costhetaRp1;
				observables[ TString::Format("costhetaR_%d_%d_%d", iHemisphere, iGeneration, hemiBalanceMode) ]  = costhetaR;


				////////////////////////////////////////////////////////////
				// Let's calculate some observables in this rest frame
				////////////////////////////////////////////////////////////

				// Angles

				// observables[ TString::Format("cosTheta_%d_%d_%d", iHemisphere, iGeneration, hemiBalanceMode) ]     = fabs( leg1.Vect().Unit().Dot( tmpBoostVector.Unit() ));
				// observables[ TString::Format("dPhi_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode) ]        = fabs( leg1.Vect().DeltaPhi( tmpBoostVector ));
				// observables[ TString::Format("dPhiVis_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode)  ]    = fabs( (leg1Vis+leg2Vis).Vect().DeltaPhi( tmpBoostVector ));

				// // Scale

				// observables[ TString::Format("M_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode) ]  = (leg1+leg2).M();

				// // Acoplanarity-type angles

				// decayPlaneVector1 = leg1_1.Vect().Cross(leg1_2.Vect());
				// decayPlaneVector2 = leg2_1.Vect().Cross(leg2_2.Vect());

				// observables[ TString::Format("dPhiDecay_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode)  ]  = decayPlaneVector1.Angle(decayPlaneVector2);
				// observables[ TString::Format("cosThetaDecay_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode)  ]  = leg1_1.Vect().Dot(decayPlaneVector2);
				// if( observables[ TString::Format("cosThetaDecay_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode)  ] < 0.0 && 
				// 	observables[ TString::Format("dPhiDecay_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode)      ] > 0.0 ){
				// 	observables[ TString::Format("dPhiDecay_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode) ] *= -1.0;
				// 	observables[ TString::Format("dPhiDecay_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode) ] +=  TMath::Pi()*2. ; 
				// }

				// // Gamma observable
				// observables[ TString::Format("gamma_%d_%d_%d",iHemisphere, iGeneration, hemiBalanceMode) ] = 
				// 1./pow( (1.-leg1.BoostVector().Mag2())*(1.-leg2.BoostVector().Mag2()),1./4. );

				// // Let's get the final mass scale
				// observables[ TString::Format("MDecay_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode) ]  = (leg1).M();



			}
	}

	unboostParticles(true, 0, -1, 0, 0); 
	return;

}




TLorentzVector Root::TRJigsaw::sumParticles(bool includeInvisible, int startingParticle = 0, int endingParticle = -1, int generation = 0, int hemisphere = 0){

	TLorentzVector tmpSum(0,0,0,0);

	if( generation >= int(hemisphereConfig[hemisphere].size())  ) return tmpSum;

	if(endingParticle == -1) endingParticle = int(hemisphereConfig[hemisphere][generation].size());

	for( int iparticle=startingParticle; iparticle<endingParticle; iparticle++){
		if(hemisphereConfig[hemisphere][generation][iparticle]=="nu" && includeInvisible){
			for( unsigned int jparticle = 0; jparticle < invParticles.size(); jparticle++){
				// //std::cout << jparticle << "/" << invParticles.size() << std::endl;
				if(invParticles.at(jparticle).hemisphere != hemisphere ) continue; // IS THIS RIGHT??????????????????????
				tmpSum = tmpSum + invParticles.at(jparticle).particleMomentumForBoosting;
			}
			continue;
		}
		for( unsigned int jparticle = 0; jparticle < visParticles.size(); jparticle++){
			// //std::cout << jparticle << "/" << visParticles.size() << std::endl;
			if(visParticles.at(jparticle).particleType != hemisphereConfig[hemisphere][generation][iparticle] ) continue;
			if(visParticles.at(jparticle).hemisphere != hemisphere ) continue;
			tmpSum = tmpSum + visParticles.at(jparticle).particleMomentumForBoosting;
		}
	}
	
	return tmpSum;

}

void Root::TRJigsaw::boostParticles(TVector3 boost, bool includeInvisible = true, int startingParticle = 0, int endingParticle = -1, int generation = 0, int hemisphere = 0){

	TLorentzVector tmpSum(0,0,0,0);
	if(hemisphere){
	  if(endingParticle == -1) endingParticle = int(hemisphereConfig[hemisphere][generation].size());
		for(int iparticle=startingParticle; iparticle<endingParticle; iparticle++){
			if(hemisphereConfig[hemisphere][generation][iparticle]=="nu" && includeInvisible){
				for( unsigned int jparticle = 0; jparticle < invParticles.size(); jparticle++){
					if(invParticles.at(jparticle).hemisphere != hemisphere ) continue;
					(invParticles.at(jparticle).particleMomentumForBoosting).Boost(boost);
				}
			}
			for( unsigned int jparticle = 0; jparticle < visParticles.size(); jparticle++){
				if(visParticles.at(jparticle).particleType != hemisphereConfig[hemisphere][generation][iparticle] ) continue;
				if(visParticles.at(jparticle).hemisphere != hemisphere ) continue;
				(visParticles.at(jparticle).particleMomentumForBoosting).Boost(boost);
			}
		}
	} else {
		boostParticles(boost,includeInvisible,startingParticle,endingParticle,generation,1);
		boostParticles(boost,includeInvisible,startingParticle,endingParticle,generation,2);
	}

	return;
}


void Root::TRJigsaw::unboostParticles(bool includeInvisible = true, int startingParticle = 0, int endingParticle = -1, int generation = 0, int hemisphere = 0){

	if(hemisphere){
	  if(endingParticle == -1) endingParticle = int(hemisphereConfig[hemisphere][generation].size());
		for( int iparticle=startingParticle; iparticle<endingParticle; iparticle++){
			if(hemisphereConfig[hemisphere][generation][iparticle]=="nu" && includeInvisible){
				for( unsigned int jparticle = 0; jparticle < invParticles.size(); jparticle++){
					if(invParticles.at(jparticle).hemisphere != hemisphere ) continue;
					invParticles.at(jparticle).particleMomentumForBoosting = invParticles.at(jparticle).particleMomentum;
				}
			}
			for( unsigned int jparticle = 0; jparticle < visParticles.size(); jparticle++){
				if(visParticles.at(jparticle).particleType != hemisphereConfig[hemisphere][generation][iparticle] ) continue;
				if(visParticles.at(jparticle).hemisphere != hemisphere ) continue;
				visParticles.at(jparticle).particleMomentumForBoosting = visParticles.at(jparticle).particleMomentum;
			}
		}
	} else {
		unboostParticles(includeInvisible,startingParticle,endingParticle,generation,1);
		unboostParticles(includeInvisible,startingParticle,endingParticle,generation,2);
	}

	return;
}

TLorentzVector Root::TRJigsaw::findTruParticleMomentum(TString particleType, int hemisphere = 0){

	// Tool to get momentum of a given truth particle

	for( unsigned int jparticle = 0; jparticle < truParticles.size(); jparticle++){
			if(truParticles.at(jparticle).particleType != particleType ) continue;
			if(truParticles.at(jparticle).hemisphere != hemisphere ) continue;
			return truParticles.at(jparticle).particleMomentum;
	}

	return TLorentzVector(0.,0.,0.,0.);

}



void Root::TRJigsaw::bookHist1D(TString expression, 
	float xbin=0, float xlow=0, float xhigh=0){

	 // Function to declare a histogram from the expression and place it 
	 // either into the 1d or 2d map of hists for use later

	Hists1D[expression] = new TH1D(expression,expression,xbin,xlow,xhigh);
	return;
}


void Root::TRJigsaw::bookHist2D(TString expression1, TString expression2,
	float xbin=0, float xlow=0, float xhigh=0,
	float ybin=0, float ylow=0, float yhigh=0 ){

	TString name = TString::Format("%s_vs_%s",expression2.Data(),expression1.Data() );
	Hists2D[name] = new TH2D(name,name,xbin,xlow,xhigh,ybin,ylow,yhigh);
	return;

}


void Root::TRJigsaw::fillHists(float weight=1.){

   // Function to read all the histogram objects in the histogram 
   // stores and to fill them using the expression(==key)
   // The expression should probably look like "observable[blah]_vs_observable[blah]"

	for(std::map<TString, TH1D* >::const_iterator i = Hists1D.begin(); i != Hists1D.end(); ++i)
	{
	    TString k = i->first;
	    TH1D* hist = i->second;
	    hist->Fill( observables[k] , weight );
	}


    std::string s, var1, var2;

	for(std::map<TString, TH2D* >::const_iterator i = Hists2D.begin(); i != Hists2D.end(); ++i)
	{
	    TString k = i->first;

	    s = k.Data();
	    var1 = s.substr(0, s.find("_vs_")  );
	    s.erase(0, s.find("_vs_") + 4 );
	    var2 = s.substr(0, s.find("_vs_")  );

	    // std::cout << var1 << " " << var2 << std::endl;

	    TH2D* hist = i->second;
	    hist->Fill( observables[var2], observables[var1] , weight );
	}

	return;

}


void Root::TRJigsaw::writeHists(TFile* f){

   // Function to read all the histogram objects in the histogram 
   // stores and to fill them using the expression(==key)
   // The expression should probably look like "observable[blah]_vs_observable[blah]"

	f->cd();

	std::cout << "Writing Hists:" << std::endl;

	for(std::map<TString, TH1D* >::const_iterator i = Hists1D.begin(); i != Hists1D.end(); ++i)
	{
	    TH1D* hist = i->second;
	    std::cout << i->first << std::endl;
	    hist->Write();
	}

	for(std::map<TString, TH2D* >::const_iterator i = Hists2D.begin(); i != Hists2D.end(); ++i)
	{
	    TH2D* hist = i->second;
	    std::cout << i->first << std::endl;
	    hist->Write();
	}


	f->Close();

	return;

}

