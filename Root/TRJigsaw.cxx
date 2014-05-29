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


ClassImp(Root::TRJigsaw)


//=============================================================================
// Constructor
//=============================================================================
Root::TRJigsaw::TRJigsaw(const char* name) :
TNamed(name,"notitle")
{
	visParticles.clear();
	 //METVector;
	hemiBalanceMode = 0;
}

Root::TRJigsaw::~TRJigsaw()
{
}




// void Root::TRJigsaw::initialize(string filename){

//    pairProd = true; // read in from conf
//    nDecays = 2; // read in from conf

//    decayModel.clear()

//    std::vector<TString> decayLine;

//    for(int iDecay=0; iDecay<nDecays; iDecay++){
//       decayLine.clear();
//       decayLine.push_back("t");
//       decayLine.push_back("b");
//       decayLine.push_back("w");

//       decayModel.push_back(decayLine);
//    }

// }

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

	// boostParticles(-BL, false, 0, -1, 0, 0); 
	H1.Boost(-BL);
	H2.Boost(-BL);

	if(hemiBalanceMode==0){ // for ttbar, this is a top symmetry mode


		//Now, we need to 'guess' the invariant mass of the weakly interacting system
		double HM2min = std::min(H1.M2(), H2.M2() );
		double HM2max = std::max(H1.M2(), H2.M2() );

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

	} else if(hemiBalanceMode==1){

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

	} else { assert(0); }

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


void Root::TRJigsaw::getObservables(){

		///////////////////////////////////////////////////////////
		// This thing is supposed to be recursive, so let's do it
		///////////////////////////////////////////////////////////

		// Let's declare some stuff ahead of the loop

	TVector3 tmpBoostVector;
	TLorentzVector leg1, leg1Vis;
	TLorentzVector leg2, leg2Vis;
	TLorentzVector leg1_1, leg1_2;
	TLorentzVector leg2_1, leg2_2;
	TVector3 decayPlaneVector1;
	TVector3 decayPlaneVector2;

		// This just makes things easier for the first vertex
	std::vector<TString> dummyVector;
	dummyVector.push_back("system");
	hemisphereConfig[0].push_back(dummyVector);


		///////////////////////////////////////////////////////////
		// Looping over hemispheres and generations for each
		// hemisphere. 
		// iHemisphere = 0: H1+H2 Frame 
		///////////////////////////////////////////////////////////

	for( int iHemisphere = 0; iHemisphere < 3; iHemisphere++){
			for( int iGeneration = 0; iGeneration < hemisphereConfig[iHemisphere].size(); iGeneration++ ){ // iGeneration is index in config file line

				// Really can't boost the first vertex more than once

				if(iHemisphere==0) assert(hemisphereConfig[0].size()==1) ;

				////////////////////////////////////////////////////////////
				// Let's start by boosting into this frame
				////////////////////////////////////////////////////////////

				// Calculate total momenta for children at vertex
				// Also momenta of the N+1 step for decay angles

				if(iHemisphere==0){
					leg1    = sumParticles(true, 0, -1, 0, 1); // t
					leg1Vis = sumParticles(false, 0, -1, 0, 1);
					leg2    = sumParticles(true, 0, -1, 0, 2);
					leg2Vis = sumParticles(false, 0, -1, 0, 2);
					leg1_1  = sumParticles(true, 0, -1, 1, 1); // W
					leg2_1  = sumParticles(true, 0, -1, 1, 2);

				} else {
					leg1    = sumParticles(true, 0, 1, iGeneration, iHemisphere);  //b - leg 1 assumed to be stable - ORDER CONFIG FILE APPROPRIATELY
					leg1Vis = sumParticles(false, 0, 1, iGeneration, iHemisphere);
					leg2    = sumParticles(true, 1, -1, iGeneration, iHemisphere); // W
					leg2Vis = sumParticles(false, 1, -1, iGeneration, iHemisphere);
					leg1_1  = leg1;
					leg2_1  = sumParticles(true, 0, 1, iGeneration+1, iHemisphere); // decay of W
				}

				// Assuming 4mom conservation to get the other leg

				leg1_2 = leg1 - leg1_1;
				leg2_2 = leg2 - leg2_1; // eg p_nu = p_W - p_lep

				// Calculate the boost of this whole system

				tmpBoostVector = (leg1+leg2).BoostVector();

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

				// Now everything is in the rest frame of the vertex!

				////////////////////////////////////////////////////////////
				// Let's calculate some observables in this rest frame
				////////////////////////////////////////////////////////////

				// angles
				observables[ TString::Format("cosTheta_%d_%d_%d", iHemisphere, iGeneration, hemiBalanceMode) ]     = fabs( leg1.Vect().Unit().Dot( tmpBoostVector.Unit() ));
				observables[ TString::Format("dPhi_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode) ]        = fabs( leg1.Vect().DeltaPhi( tmpBoostVector ));
				observables[ TString::Format("dPhiVis_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode)  ]    = fabs( (leg1Vis+leg2Vis).Vect().DeltaPhi( tmpBoostVector ));
				// scale
				observables[ TString::Format("M_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode) ]  = (leg1+leg2).M();

				// acoplanarity-type angles
				decayPlaneVector1 = leg1_1.Vect().Cross(leg1_2.Vect());
				decayPlaneVector2 = leg2_1.Vect().Cross(leg2_2.Vect());

				observables[ TString::Format("dPhiDecay_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode)  ]  = decayPlaneVector1.Angle(decayPlaneVector2);
				observables[ TString::Format("codThetaDecay_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode)  ]  = leg1_1.Vect().Dot(decayPlaneVector2);
				if( observables[ TString::Format("codThetaDecay_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode)  ] < 0.0 && 
					observables[ TString::Format("dPhiDecay_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode)      ] > 0.0 ){
					observables[ TString::Format("dPhiDecay_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode) ] *= 1.0;
				observables[ TString::Format("dPhiDecay_%d_%d_%d" , iHemisphere, iGeneration, hemiBalanceMode) ] +=  TMath::Pi()*2. ; 
			}

				// gamma observable
			observables[ TString::Format("gamma_%d_%d_%d",iHemisphere, iGeneration, hemiBalanceMode) ] = 
			1./pow( (1.-leg1.BoostVector().Mag2())*(1.-leg2.BoostVector().Mag2()),1./4. );

		}
	}

	return;

}



void Root::TRJigsaw::bookHist(int dimension, TString expression, 
	float xbin=0, float xlow=0, float xhigh=0,
	float ybin=0, float ylow=0, float yhigh=0 ){

	 // Function to declare a histogram from the expression and place it 
	 // either into the 1d or 2d map of hists for use later

	if(dimension==1){
		Hists1D[expression] = new TH1D(expression,expression,xbin,xlow,xhigh);
	}
	else if(dimension==2){
		Hists2D[expression] = new TH2D(expression,expression,xbin,xlow,xhigh,ybin,ylow,yhigh);
	}

	return;

}

TLorentzVector Root::TRJigsaw::sumParticles(bool includeInvisible, int startingParticle = 0, int endingParticle = -1, int generation = 0, int hemisphere = 0){

	TLorentzVector tmpSum(0,0,0,0);
	if(endingParticle == -1) endingParticle = hemisphereConfig[hemisphere][generation].size();
	for(int iparticle=startingParticle; iparticle<endingParticle; iparticle++){
		if(hemisphereConfig[hemisphere][generation][iparticle]=="nu" && includeInvisible){
			for( int jparticle = 0; jparticle < invParticles.size(); jparticle++){
				if(invParticles.at(jparticle).hemisphere != 1 ) continue;
				tmpSum = tmpSum + invParticles.at(jparticle).particleMomentumForBoosting;
			}
		}
		for( int jparticle = 0; jparticle < visParticles.size(); jparticle++){
			if(visParticles.at(jparticle).particleType != hemisphereConfig[hemisphere][generation][iparticle] ) continue;
			if(visParticles.at(jparticle).hemisphere != 1 ) continue;
			tmpSum = tmpSum + invParticles.at(jparticle).particleMomentumForBoosting;
		}
	}
	
	return tmpSum;

}

void Root::TRJigsaw::boostParticles(TVector3 boost, bool includeInvisible = true, int startingParticle = 0, int endingParticle = -1, int generation = 0, int hemisphere = 0){

	TLorentzVector tmpSum(0,0,0,0);
	if(hemisphere){
		if(endingParticle == -1) endingParticle = hemisphereConfig[hemisphere][generation].size()
			for(int iparticle=startingParticle; iparticle<endingParticle; iparticle++){
				if(hemisphereConfig[hemisphere][generation][iparticle]=="nu" && includeInvisible){
					for( int jparticle = 0; jparticle < invParticles.size(); jparticle++){
						if(invParticles.at(jparticle).hemisphere != 1 ) continue;
						(invParticles.at(jparticle).particleMomentumForBoosting).Boost(boost);
					}
				}
				for( int jparticle = 0; jparticle < visParticles.size(); jparticle++){
					if(visParticles.at(jparticle).particleType != hemisphereConfig[hemisphere][generation][iparticle] ) continue;
					if(visParticles.at(jparticle).hemisphere != 1 ) continue;
					(invParticles.at(jparticle).particleMomentumForBoosting).Boost(boost);
				}
			}
		} else {
			boostParticles(boost,includeInvisible,startingParticle,endingParticle,generation,1)
			boostParticles(boost,includeInvisible,startingParticle,endingParticle,generation,2)
		}

		return;
	}




// void Root::TRJigsaw::fillHists(){

//    // Function to read all the histogram objects in the histogram 
//    // stores and to fill them using the expression(==key)
//    // The expression should probably look like "observable[blah]_vs_observable[blah]"

//    for( ){ // Loop through keys in 1DHists



//    }  


// }






// //combine run number's histograms, periods, everything, then remove trace of second run number
// void Root::TRJigsaw::MergeMCRunNumbers(Int_t mcRunNumber1, Int_t mcRunNumber2, Int_t newMCRunNumber) {
//    if(m_isInitialized) {
//       Error("MergeMCRunNumbers","You cannot MergeMCRunNumbers after initializing the tool. Reorder your code!");
//       throw std::runtime_error("Throwing 5: You cannot MergeMCRunNumbers after initializing the tool. Reorder your code!");
//    }
//    if(newMCRunNumber==0) newMCRunNumber=mcRunNumber1;
//    m_mcMerges[newMCRunNumber][mcRunNumber1]=mcRunNumber2;

// }

