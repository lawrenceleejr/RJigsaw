// Dear emacs, this is -*-c++-*-

#ifndef PILEUPREWEIGHTING_PILEUPREWEIGHTINGALG_H
#define PILEUPREWEIGHTING_PILEUPREWEIGHTINGALG_H

/**
   @class PileupReweightingAlg
   @brief Algorithm to get the calculated MC pilup weight and attach it to the event as UserData.

   @author Karsten Koeneke

*/

/******************************************************************************
Name:        PileupReweightingAlg

Author:      Karsten Koeneke (CERN)
Created:     April 2011

Description: Algorithm to get the calculated MC pilup weight and attach it to the event as UserData.
******************************************************************************/


#include <string>

#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ServiceHandle.h"

class IUserDataSvc;
namespace Root{
  class TPileupReweighting;
}

class PileupReweightingAlg : public AthAlgorithm {

public: 
  /** Standard constructor */
  PileupReweightingAlg(const std::string& name, ISvcLocator* pSvcLocator);

  /** Standard destructor */
  virtual ~PileupReweightingAlg();
  
public:
  /** Gaudi Service Interface method implementations: initialize */
  StatusCode initialize();

  /** Gaudi Service Interface method implementations: execute */
  StatusCode execute();

  /** Gaudi Service Interface method implementations: finalize */
  StatusCode finalize();



  // Private members
private:
  /** Get a handle on the UserDataSvc */
  ServiceHandle<IUserDataSvc> m_userDataSvc;
  
  
  /** Get a handle on the TPileupReweighting class */
  Root::TPileupReweighting* m_tPileUp;
  

  /** The name for the UserData event pileup weight variable */
  std::string m_userDataWeightName;


  /** The name for the EventInfo container */
  std::string m_eventInfoName;


  /** The name of the data root file */
  std::string m_dataRootFileName;

  /** The name of the data histogram */
  std::string m_dataRootHistName;


  /** The name of the MC root file */
  std::string m_mcRootFileName;

  /** The name of the MC histogram */
  std::string m_mcRootHistName;


  /** Event counter */
  unsigned long m_nEventsProcessed;

};

#endif
