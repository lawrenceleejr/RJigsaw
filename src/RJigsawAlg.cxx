// Dear emacs, this is -*-c++-*-
/******************************************************************************
Name:        RJigsawAlg

Author:      Lawrence Lee
Created:     April 2011

Description: Algorithm to get the calculated MC pilup weight and attach it to the event as UserData.
******************************************************************************/

// Preprocessor magic for debugging
#define XXX std::cout << "I am here: " << __FILE__ << ":" << __LINE__ << std::endl;

// This class' header
#include "RJigsaw/RJigsawAlg.h"

// Include the ROOT class from this package
#include "RJigsaw/TRJigsaw.h"


// STL includes
#include <string>

// ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TAxis.h>

// Gaudi inlcudes
#include "GaudiKernel/ServiceHandle.h"

// Include the UserDataSvc
#include "AthenaKernel/IUserDataSvc.h"
#include "AthenaKernel/errorcheck.h"





//=============================================================================
// Constructor
//=============================================================================
RJigsawAlg::RJigsawAlg(const std::string& name, ISvcLocator* pSvcLocator) :
  AthAlgorithm(name, pSvcLocator),
  m_userDataSvc( "UserDataSvc", name ),
  m_tRJigsaw(0)
{

}



//=============================================================================
// Destructor
//=============================================================================
RJigsawAlg::~RJigsawAlg()
{
}

StatusCode RJigsawAlg::initialize()
{
  // Print the used configuration
  if (msgLvl(MSG::INFO))
    {
      msg(MSG::INFO) 
        << "==> initialize " << name() << "..."
        << endreq;
      
      // // Print out the used configuration
      // msg(MSG::INFO) << " using userDataSvc                   = " << m_userDataSvc << endreq;
      // msg(MSG::INFO) << " using userDataEventRJigsawWeightName = " << m_userDataWeightName << endreq;

    }

  // The standard status code
  StatusCode sc(StatusCode::SUCCESS);

  
  // Get the UserData service
  CHECK( m_userDataSvc.retrieve() );


  // Initialize the counters to zero
  m_nEventsProcessed = 0;


  //-------------------------------------------------------
  // Get the ROOT class and initialize it
  //-------------------------------------------------------
  m_tRJigsaw = new Root::TRJigsaw( (((std::string)"T")+name()).c_str() );
  int isGood = m_tRJigsaw->initialize();

  if ( isGood == 0 )
    {
      REPORT_MESSAGE(MSG::DEBUG) << "Initialization of TRJigsaw successful";
    }
  else 
    {
      REPORT_MESSAGE(MSG::ERROR) << "Unrecognized return code! Exiting!";
      sc = StatusCode::FAILURE;
    }

  return sc;
}





//=============================================================================
// Athena execute event method
//=============================================================================
StatusCode RJigsawAlg::execute() 
{
  // Declare the simple StatusCode
  StatusCode sc(StatusCode::SUCCESS);


  // Simple status message at the beginning of each event execute
  if (msgLvl(MSG::DEBUG))
    {
      msg(MSG::DEBUG) 
        << "==> execute " << name() << " on " << m_nEventsProcessed << ". event..."
        << endreq;
    }
  ++m_nEventsProcessed;
  
  return sc;
}




//=============================================================================
// Athena finalize method
//=============================================================================
StatusCode RJigsawAlg::finalize() 
{
  // Declare the simple StatusCode
  StatusCode sc(StatusCode::SUCCESS);

  // Print info messages
  if (msgLvl(MSG::INFO))
    {
      msg(MSG::INFO) 
        << "==> finalize " << name() << "...\n"
        << "***************************************************************\n"
        << endreq;
      msg(MSG::INFO) 
        << " Number of processed events:  " << m_nEventsProcessed
        << endreq;
    }

  // Delete the pointers
  if ( m_tRJigsaw ) delete m_tRJigsaw;

  return sc;
}

