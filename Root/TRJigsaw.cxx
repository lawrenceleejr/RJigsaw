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
      TLorentzVector particleMomentum
      ){

   particleClass particle;
   particle.particleType = particleType;
   particle.particleMomentum = particleMomentum;
   particle.particleMomentumForBoosting = particleMomentum;
   visParticles.push_back( particle );

   return;

}

void Root::TRJigsaw::addTruParticle(
      TString particleType,
      TLorentzVector particleMomentum
      ){

   particleClass particle;
   particle.particleType = particleType;
   particle.particleMomentum = particleMomentum;
   particle.particleMomentumForBoosting = particleMomentum;
   truParticles.push_back( particle );

   return;

}


void Root::TRJigsaw::guessInvParticles(){

   invParticles.clear();

   TLorentzVector Inv1, Inv2;
   bool Inv1Empty = 1, Inv2Empty = 1;

   if(hemiBalanceMode==-1){ // truth

      for( int iparticle = 0; iparticle < truParticles.size(); iparticle++){
         if(truParticles.at(iparticle).particleType != "inv"){ continue;}
         
         if(Inv1Empty){Inv1 = truParticles.at(iparticle).particleMomentumForBoosting; Inv1Empty=0;}
         else if(Inv2Empty) {Inv2 = truParticles.at(iparticle).particleMomentumForBoosting; Inv2Empty=0;}
         else {assert(0); }
      }

      particleClass invParticle1;
      invParticle1.particleType = TString("inv");
      invParticle1.particleMomentum = Inv1;
      invParticle1.particleMomentumForBoosting = Inv1;
      particleClass invParticle2;
      invParticle2.particleType = TString("inv");
      invParticle2.particleMomentum = Inv2;
      invParticle2.particleMomentumForBoosting = Inv2;

      invParticles.push_back(invParticle1);
      invParticles.push_back(invParticle2);

      return;

   }

   //get b-quarks in lab frame

   TLorentzVector B1,B2;   
   bool B1Empty=1, B2Empty=1;

   for( int iparticle = 0; iparticle < visParticles.size(); iparticle++){
      if(visParticles.at(iparticle).particleType != "b"){ continue;}
      
      if(B1Empty){B1 = visParticles.at(iparticle).particleMomentumForBoosting; B1Empty = 0;}
      else if(B2Empty) {B2 = visParticles.at(iparticle).particleMomentumForBoosting; B2Empty = 0;}
      //else {assert(0); }
   }

   //get leptons in lab frame

   TLorentzVector L1,L2;
   bool L1Empty=1, L2Empty=1;

   for( int iparticle = 0; iparticle < visParticles.size(); iparticle++){
      if(visParticles.at(iparticle).particleType != "l"){ continue;}
      
      if(L1Empty){L1 = visParticles.at(iparticle).particleMomentumForBoosting; L1Empty=0;}
      else if(L2Empty) {L2 = visParticles.at(iparticle).particleMomentumForBoosting; L2Empty=0;}
      //else {assert(0); }
   }

   //We separate the objects into decay hemispheres
   TLorentzVector H1 = (B1+L1);
   TLorentzVector H2 = (B2+L2);

   //Now, we must perform longitudinal boost from lab frame to CMz frame
   TVector3 BL = H1.Vect()+H2.Vect();
   BL.SetX(0.0);
   BL.SetY(0.0);
   BL = (1./(H1.E()+H2.E()))*BL;


   B1.Boost(-BL);
   L1.Boost(-BL);
   B2.Boost(-BL);
   L2.Boost(-BL);
   H1.Boost(-BL);
   H2.Boost(-BL);

   if(hemiBalanceMode==0){ // for ttbar, this is a top symmetry mode

      //Now, we need to 'guess' the invariant mass of the weakly interacting system
      double HM2min = H1.M2();
      double HM2max = H2.M2();
      if(HM2min > HM2max){
         double temp = HM2max;
         HM2max = HM2min;
         HM2min = temp;
      }

      double Minv2 = (H1+H2).M2() - 4.*HM2min;

      TVector3 PT = *METVector+(H1+H2).Vect();
      double Einv2 = METVector->Mag2() + Minv2;
      PT.SetZ(0.0);
      TVector3 BT = (1./( (H1+H2).E() + sqrt(Einv2) ))*PT;
      B1.Boost(-BT);
      L1.Boost(-BT);
      B2.Boost(-BT);
      L2.Boost(-BT);
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

   }// else { assert(0); }

   particleClass invParticle1;
   invParticle1.particleType = TString("inv");
   invParticle1.particleMomentum = Inv1;
   invParticle1.particleMomentumForBoosting = Inv1;
   particleClass invParticle2;
   invParticle2.particleType = TString("inv");
   invParticle2.particleMomentum = Inv2;
   invParticle2.particleMomentumForBoosting = Inv2;

   invParticles.push_back(invParticle1);
   invParticles.push_back(invParticle2);

   return;

}

// void Root::TRJigsaw::getObservables(){
//    // just get all visible particles and organize by place in decay chain


//    std::map<TString, std::vector<TLorentzVector> > LabFrameParticles;
//    for( int iparticle = 0; iparticle < visParticles.size(); iparticle++){
//       LabFrameParticles[visParticles.at(iparticle).particleType].push_back(visParticles.at(iparticle).particleMomentum);
//    }
//    LabFrameParticles["inv"].push_back(invParticles[0].particleMomentum);
//    LabFrameParticles["inv"].push_back(invParticles[1].particleMomentum);

// }



void Root::TRJigsaw::getObservables(){



   // //get b-quarks in lab frame

   // TLorentzVector B1,B2;
   // bool B1Empty=1, B2Empty=1;

   // for( int iparticle = 0; iparticle < visParticles.size(); iparticle++){
   //    if(visParticles.at(iparticle).particleType != "b"){ continue;}
      
   //    if(B1Empty){B1 = visParticles.at(iparticle).particleMomentum; B1Empty=1;}
   //    else if(B2Empty) {B2 = visParticles.at(iparticle).particleMomentum; B2Empty=1;}
   //    else {assert(0); }
   // }

   // //get leptons in lab frame

   // TLorentzVector L1,L2;
   // bool L1Empty=1, L2Empty=1;

   // for( int iparticle = 0; iparticle < visParticles.size(); iparticle++){
   //    if(visParticles.at(iparticle).particleType != "l"){ continue;}
      
   //    if(L1Empty){L1 = visParticles.at(iparticle).particleMomentum; L1Empty=1;}
   //    else if(L2Empty) {L2 = visParticles.at(iparticle).particleMomentum;L2Empty=1; }
   //    else {assert(0); }
   // }

   // assert(invParticles.size()==2);
   // TLorentzVector NU1 = invParticles[0].particleMomentumForBoosting;
   // TLorentzVector NU2 = invParticles[1].particleMomentumForBoosting;


   // let's assume I have all the particles I need and they have hemisphere associations


    ////////////////////////////////////////////////////////////
    // Contruct our hemisphere momenta
    ////////////////////////////////////////////////////////////

    // for each hemisphere, take the first line
    // TLorentzVector H_N = sum of these particles
    TLorentzVector H1    = sumParticles(true, 0, -1, 0, 1);
    TLorentzVector H1Vis = sumParticles(false, 0, -1, 0, 1);
    TLorentzVector H2    = sumParticles(true, 0, -1, 0, 2);
    TLorentzVector H2Vis = sumParticles(false, 0, -1, 0, 2);

    ////////////////////////////////////////////////////////////
    // Let's move from the lab frame to the first rest frame
    ////////////////////////////////////////////////////////////

    // Calculate the boost vector
    TVector3 BltoCM = (H1+H2).BoostVector();

    // Do the boosting
    for( int jparticle = 0; jparticle < invParticles.size(); jparticle++){ invParticles.at(jparticle).particleMomentumForBoosting.Boost(-BltoCM);  }
    for( int jparticle = 0; jparticle < visParticles.size(); jparticle++){ visParticles.at(jparticle).particleMomentumForBoosting.Boost(-BltoCM);  }
    H1.Boost(-BltoCM);
    H2.Boost(-BltoCM);

    // Now everything is in the rest frame of H1+H2!

    ////////////////////////////////////////////////////////////
    // Let's calculate some observables in this frame
    ////////////////////////////////////////////////////////////

    //angles
    observables[ "costheta_"+TString(hemiBalanceMode)+TString(0) ]     = fabs( H1.Vect().Unit().Dot( BltoCM.Unit() ));
    observables[ "dphi_"+TString(hemiBalanceMode)+TString(0) ]        = fabs( H1.Vect().DeltaPhi( BltoCM ));
    observables[ "dphiM_"+TString(hemiBalanceMode)+TString(0) ]        = fabs( (H1Vis+H2Vis).Vect().DeltaPhi( BltoCM ));
    //scale
    observables[ "M_"+TString(hemiBalanceMode)+TString(0) ]  = (H1+H2).M();



    // LL - maybe this needs to wait to be the first thing in the generation loop...

    for(int iGeneration = 1; iGeneration < hemisphere1Config.size()+1; iGeneration++){


    }

    // It's really important that the particles from later decays show up later in the config line
    TLorentzVector tmpIntermediateParticle1 = sumParticles(true, 1, -1, 0, 1);
    TLorentzVector tmpIntermediateParticle2 = sumParticles(true, 1, -1, 0, 2);

    TLorentzVector tmpSurvivorParticle1 = sumParticles(true, 0, 1, 0, 1)
    TLorentzVector tmpSurvivorParticle2 = sumParticles(true, 0, 1, 0, 2)

    //vector normal to decay plane of CM frame
    TVector3 vNORM_CM_T1 = tmpSurvivorParticle1.Vect().Cross(tmpIntermediateParticle1.Vect());
    TVector3 vNORM_CM_T2 = tmpSurvivorParticle2.Vect().Cross(tmpIntermediateParticle2.Vect());

    observables["dphi_1_2_"+TString(hemiBalanceMode)+TString(0) ]  = vNORM_CM_T1.Angle(vNORM_CM_T2);
    observables["dot_dphi_1_2_"+TString(hemiBalanceMode)+TString(0) ]  = tmpSurvivorParticle1.Vect().Dot(vNORM_CM_T2);
    if(observables["dot_dphi_1_2_"+TString(hemiBalanceMode)+TString(0) ]  < 0.0 && observables["dphi_1_2_"+TString(hemiBalanceMode)+TString(0) ]  > 0.0){
      observables["dphi_1_2_"+TString(hemiBalanceMode)+TString(0) ]  = TMath::Pi()*2. - observables["dphi_1_2"];
    }

    ///////////////////////////////////////////
   // for each hemisphere
   // for each line in hemisphere config
   // TLorentzVector myBoost = sum of each of the momenta of those particles in the previous frame
   // boost all particles in that hemisphere to this frame (top rest frame)
   // get observables
    /////////////////////////////////////////

    //To the next frames!!!!
    //now, to 'top' CM frame approxs
    TVector3 tmpBoostVector1 = H1.BoostVector();
    TVector3 tmpBoostVector2 = H2.BoostVector();

    B1.Boost(-tmpBoostVector1);
    B2.Boost(-tmpBoostVector2);
    L1.Boost(-tmpBoostVector1);
    L2.Boost(-tmpBoostVector2);
    NU1.Boost(-tmpBoostVector1);
    NU2.Boost(-tmpBoostVector2);

    observables["costhetaT1_"+TString(hemiBalanceMode) ]  = B1.Vect().Unit().Dot(tmpBoostVector1.Unit());
    observables["costhetaT2_"+TString(hemiBalanceMode) ]  = B2.Vect().Unit().Dot(tmpBoostVector2.Unit());

    //vectors normal to decay planes of T frames
    TVector3 vNORM_T1_B = B1.Vect().Cross(tmpBoostVector1);
    TVector3 vNORM_T2_B = B2.Vect().Cross(tmpBoostVector2);
    //vectors normal to W decay planes in T frames
    TVector3 vNORM_T1_W = L1.Vect().Cross(NU1.Vect());
    TVector3 vNORM_T2_W = L2.Vect().Cross(NU2.Vect());


    observables["dphi_W_T1_"+TString(hemiBalanceMode) ]  = vNORM_T1_W.Angle(vNORM_T1_B);
    observables["dphi_W_T2_"+TString(hemiBalanceMode) ]  = vNORM_T2_W.Angle(vNORM_T2_B);
    observables["dot_dphi_W_T1_"+TString(hemiBalanceMode) ]  = L1.Vect().Dot(vNORM_T1_B);
    observables["dot_dphi_W_T2_"+TString(hemiBalanceMode) ]  = L2.Vect().Dot(vNORM_T2_B);
    if(observables["dot_dphi_W_T1_"+TString(hemiBalanceMode) ]  < 0.0 && observables["dphi_W_T1_"+TString(hemiBalanceMode) ]  > 0.0){
      observables["dphi_W_T1_"+TString(hemiBalanceMode) ]  = TMath::Pi()*2. - observables["dphi_W_T1_"+TString(hemiBalanceMode) ] ;
    }
    if(observables["dot_dphi_W_T2_"+TString(hemiBalanceMode) ]  < 0.0 && observables["dphi_W_T2_"+TString(hemiBalanceMode) ]  > 0.0){
      observables["dphi_W_T2_"+TString(hemiBalanceMode) ]  = TMath::Pi()*2. - observables["dphi_W_T2_"+TString(hemiBalanceMode) ] ;
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

TLorentzVector Root::TRJigsaw::sumParticles(bool includeInvisible = true, int startingParticle = 0, int endingParticle = -1, int generation = 0, int hemisphere = 0){

    TLorentzVector tmpSum(0,0,0,0);
    if(hemisphere==1){
      if(endingParticle == -1) endingParticle = hemisphere1Config[generation].size()
      for(int iparticle=startingParticle; iparticle<endingParticle; iparticle++){
        if(hemisphere1Config[generation][iparticle]=="nu" && includeInvisible){
          for( int jparticle = 0; jparticle < invParticles.size(); jparticle++){
            if(invParticles.at(jparticle).hemisphere != 1 ) continue;
            tmpSum = tmpSum + invParticles.at(jparticle).particleMomentumForBoosting;
          }
        }
        for( int jparticle = 0; jparticle < visParticles.size(); jparticle++){
          if(visParticles.at(jparticle).particleType != hemisphere1Config[generation][iparticle] ) continue;
          if(visParticles.at(jparticle).hemisphere != 1 ) continue;
          tmpSum = tmpSum + invParticles.at(jparticle).particleMomentumForBoosting;
        }
      }
    }
    else if(hemisphere==2){
      if(endingParticle == -1) endingParticle = hemisphere2Config[generation].size()
      for(int iparticle=startingParticle; iparticle<endingParticle; iparticle++){
        if(hemisphere2Config[generation][iparticle]=="nu" && includeInvisible){
          for( int jparticle = 0; jparticle < invParticles.size(); jparticle++){
            if(invParticles.at(jparticle).hemisphere != 2 ) continue;
            tmpSum = tmpSum + invParticles.at(jparticle).particleMomentumForBoosting;
          }
        }
        for( int jparticle = 0; jparticle < visParticles.size(); jparticle++){
          if(visParticles.at(jparticle).particleType != hemisphere2Config[generation][iparticle] ) continue;
          if(visParticles.at(jparticle).hemisphere != 2 ) continue;
          tmpSum = tmpSum + invParticles.at(jparticle).particleMomentumForBoosting;
        }
      }
    }

    return tmpSum;

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

