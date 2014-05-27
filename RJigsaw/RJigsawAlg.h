// Dear emacs, this is -*-c++-*-

#ifndef RJIGSAW_RJIGSAWALG_H
#define RJIGSAW_RJIGSAWALG_H

/**
   @class RJigsawAlg
   @brief 

   @author 

*/

/******************************************************************************
Name:        RJigsawAlg

Author:      
Created:     

Description: 
******************************************************************************/


#include <string>

#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ServiceHandle.h"

class IUserDataSvc;
namespace Root{
  class TRJigsaw;
}

class RJigsawAlg : public AthAlgorithm {

public: 
  /** Standard constructor */
  RJigsawAlg(const std::string& name, ISvcLocator* pSvcLocator);

  /** Standard destructor */
  virtual ~RJigsawAlg();
  
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
  
  
  /** Get a handle on the TRJigsaw class */
  Root::TRJigsaw* m_tRJigsaw;
  

  /** Event counter */
  unsigned long m_nEventsProcessed;

};

#endif
