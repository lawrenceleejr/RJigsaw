// Dear emacs, this is -*-c++-*-

#ifndef __TRJIGSAW__
#define __TRJIGSAW__

#define RJIGSAWVERSION 1 //useful maybe for other packages to detect with. Each release I'll try to remember to increment this

/**
@class TRJigsaw
@brief Tools to calculate Jigsaw variables

@author  Lawrence Lee
Based on the original work by Chris Rogan and Paul Jackson
*/

#include "TNamed.h"
#include <TFile.h>
#include "TVectorD.h"

#include <TObject.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include <fstream>


// STL includes
#include <iostream>
#include <stdexcept>

class TH1;
class TH1D;
class TH2D;
class TTree;
class TFile;

namespace Root {

	class TRJigsaw : public TNamed {

	public: 
		/** Standard constructor */
		TRJigsaw(const char* name="TRJigsaw");

		/** Standard destructor */
		~TRJigsaw();

		/////////////////////////////////////////////////////////////////////////////////////////
		// initialize() needs two text files formatted as in the top examples in the package

		void initialize(std::string filename1,std::string filename2);

		// Resets the particle vectors

		void newEvent(){
			truParticles.clear();
			visParticles.clear();
			invParticles.clear();
			observables.clear();
		};

		// Add a truth particle to that vector

		void addTruParticle(
			TString particleType,
			TLorentzVector particleMomentum,
			int hemisphere
		);

		// Add any visible particle momentum to the visParticle vector

		void addVisParticle(
			TString particleType,
			TLorentzVector particleMomentum,
			int hemisphere
			);

		// The MET vector access functions

		void addMET(TVector3 InputMETVector){ METVector = (TVector3*) InputMETVector.Clone(); return; };
		TVector3* getMET(){ return METVector; };

		// Guess the momenta of the two escape particles

		void guessInvParticles(); 

		// Do all of the jigsaw boosting and dump the observables into a map

		void getObservables();

		// mode is the same as the generation of decay, let's say. So for ttbar, 
		// using the top symmetry will be mode 0 
		// and using the W symmetry will be mode 1

		void setHemisphereMode(int mode){ hemiBalanceMode = mode; return; };
		int getHemisphereMode(){return hemiBalanceMode;}


		/////////////////////////////////////////////
		/////////////////////////////////////////////
		// Histogram Handling
		//

		void resetHists(){
			Hists1D.clear();
			Hists2D.clear();
		}

		void bookHist1D(TString expression, 
			float xbin, float xlow, float xhigh);

		void bookHist2D(TString expression1, TString expression2,
			float xbin, float xlow, float xhigh,
			float ybin, float ylow, float yhigh );

		// Loops through all histograms and fills or writes to a file

		void fillHists(float weight);
		void writeHists(TFile* f);

		// Access to the raw observables map if needed

		std::map< TString, double > getObservablesMap(){ return observables;}

		// A handy particle class

		class particleClass {
		public:
			TString particleType;
			TLorentzVector particleMomentum;
			TLorentzVector particleMomentumForBoosting;
			int hemisphere = 0;
		};

	protected:

		// Get the momentum sum for all particles matching the critera below (hemisphere, generation, etc)
		// Uses the initialization files as a reference

		TLorentzVector sumParticles(bool includeInvisible, int startingParticle , int endingParticle , int generation , int hemisphere );

		// Boost the desired particles by -boost

		void boostParticles(TVector3 boost, bool includeInvisible, int startingParticle, int endingParticle, int generation, int hemisphere);

		TLorentzVector findTruParticleMomentum(TString particleType, int hemisphere);

		std::vector< particleClass > truParticles;
		std::vector< particleClass > visParticles;
		std::vector< particleClass > invParticles;
		TVector3* METVector;

		std::map< TString, double > observables;

		std::map< TString, TH1D* > Hists1D;
		std::map< TString, TH2D* > Hists2D;

		int hemiBalanceMode = 0;

		std::vector< std::vector<TString> > hemisphereConfig[3];


	ClassDef(TRJigsaw,1)


	}; // End: class definition

} // End: namespace Root

#endif
