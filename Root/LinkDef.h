#include <map>
#include <string>
#include "TString.h"
#include <RJigsaw/TRJigsaw.h>

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace Root ;

#pragma link C++ class Root::TRJigsaw+ ; 

#pragma link C++ class map<int,map<int,pair<unsigned int,unsigned int > > >+ ;
#pragma link C++ class map<int,pair<unsigned int,unsigned int > >+ ;
#pragma link C++ class pair<unsigned int,unsigned int>+ ;
#pragma link C++ class map<TString,TH1*>+ ;
#pragma link C++ class map<int,map<int,int> >+ ;
#pragma link C++ class map<TString,map<int, double> >+ ;
#pragma link C++ class map<int,double>+ ;
#pragma link C++ class map<TString,map<int,map<int,double> > >+;
#pragma link C++ class map<int,map<int,double> >+;
#pragma link C++ class map<TString,map<int, map<int,TH1D*> > >+;
#pragma link C++ class map<int,map<int,TH1D*> >+;
#pragma link C++ class map<int,TH1D*>+;
#pragma link C++ class map<TString,map<int, map<int,TH2D*> > >+;
#pragma link C++ class map<int,map<int,TH2D*> >+;
#pragma link C++ class map<int,TH2D*>+;
#pragma link C++ class map<TString,map<int, map<int,TH1*> > >+;
#pragma link C++ class map<int,map<int,TH1*> >+;
#pragma link C++ class map<int,TH1*>+;
//#pragma link C++ class map<TString,map<unsigned int,map<int,map<unsigned int,double> > > >+; 
//#pragma link C++ class map<unsigned int,map<int,map<unsigned int,double> > >+; 
#pragma link C++ class map<TString,map<int,map<unsigned int,double> > >+; 
#pragma link C++ class map<int,map<unsigned int,double> >+; 
#pragma link C++ class map<unsigned int,double>+; 
#pragma link C++ class map<TString,map<int,map<int,bool> > >+;
#pragma link C++ class map<int,map<int,bool> >+;
#pragma link C++ class map<int,bool>+;
#pragma link C++ class map<TString,TString>+;

#pragma link C++ class pair<int,map<int,pair<unsigned int,unsigned int > > >+ ;
#pragma link C++ class pair<int,pair<unsigned int,unsigned int > >+ ;
#pragma link C++ class pair<TString,TH1*>+ ;
#pragma link C++ class pair<int,map<int,int> >+ ;
#pragma link C++ class pair<TString,map<int, double> >+ ;
#pragma link C++ class pair<TString,map<int,map<int,double> > >+;
#pragma link C++ class pair<int,map<int,double> >+;
#pragma link C++ class pair<TString,map<int, map<int,TH1D*> > >+;
#pragma link C++ class pair<int,map<int,TH1D*> >+;
#pragma link C++ class pair<int,TH1D*>+;
#pragma link C++ class pair<TString,map<int, map<int,TH2D*> > >+;
#pragma link C++ class pair<int,map<int,TH2D*> >+;
#pragma link C++ class pair<int,TH2D*>+;
#pragma link C++ class pair<TString,map<int, map<int,TH1*> > >+;
#pragma link C++ class pair<int,map<int,TH1*> >+;
#pragma link C++ class pair<int,TH1*>+;
//#pragma link C++ class pair<TString,map<unsigned int,map<int,map<unsigned int,double> > > >+; 
//#pragma link C++ class pair<unsigned int,map<int,map<unsigned int,double> > >+; 
#pragma link C++ class pair<TString,map<int,map<unsigned int,double> > >+; 
#pragma link C++ class pair<int,map<unsigned int,double> >+; 
#pragma link C++ class pair<unsigned int,double>+; 
#pragma link C++ class pair<TString,map<int,map<int,bool> > >+;
#pragma link C++ class pair<int,map<int,bool> >+;
#pragma link C++ class pair<int,bool>+;


#endif
