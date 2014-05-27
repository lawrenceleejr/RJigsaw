
#include <TObject.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TString.h>



class particleClass : public TObject {
public:
	TString particleType;
	TLorentzVector particleMomentum;
	int hemisphere = 0;

	ClassDef(particleClass,1)
};