// Dear emacs, this is -*-c++-*-

#ifndef __TRJIGSAW__
#define __TRJIGSAW__

#define RJIGSAWVERSION 3 //useful maybe for other packages to detect with. Each release I'll try to remember to increment this

/**
   @class TRJigsaw
   @brief Tools to calculate Jigsaw variables

   @author  Lawrence Lee
            Based on the original work by Chris Rogan and Paul Jackson
*/

//#include "particleClass.h"

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

      void initialize(std::string filename);

      // initializeRJigsaw /////////////
      // Read in config file that will look like...
      // ...something

      void newEvent(){
        truParticles.clear();
        visParticles.clear();
        invParticles.clear();
        //METVector.Reset();
        observables.clear();
      };

      void addTruParticle(
        TString particleType,
        TLorentzVector particleMomentum
        );

      void addVisParticle(
        TString particleType,
        TLorentzVector particleMomentum
        );

      void guessInvParticles(); // These should be run together 
      void getObservables();


      void addMET(TVector3 InputMETVector){ METVector = (TVector3*) InputMETVector.Clone(); return; };
      TVector3* getMET(){ return METVector; };

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

      void bookHist(int dimension, TString expression, // If it's 2d, separate by a "_vs_"
              float xbin, float xlow, float xhigh,
              float ybin, float ylow, float yhigh );

      // // void fillHists()
      // // void writeHists()


      class particleClass {
      public:
        TString particleType;
        TLorentzVector particleMomentum;
        TLorentzVector particleMomentumForBoosting;
        int hemisphere = 0;
      };

  protected:

      // Class member variables

      TLorentzVector sumParticles(bool includeInvisible = true, int startingPoint = 0, int hemisphere = 0);

      std::vector< particleClass > truParticles;
      std::vector< particleClass > visParticles;
      std::vector< particleClass > invParticles;
      TVector3* METVector;

      std::map< TString, double > observables;

      std::map< TString, TH1D* > Hists1D;
      std::map< TString, TH2D* > Hists2D;


      int hemiBalanceMode = 0;

      std::vector< std::vector<TString> > hemisphereConfig[3];


  public:
      /** Initialize this class once before the event loop starts 
          If distribution information is provided, it is assumed to be 
          for the standard pileup reweighting */
      // Int_t initialize();

      // Int_t Initialize();


      ClassDef(TRJigsaw,1)


  }; // End: class definition

} // End: namespace Root

#endif
