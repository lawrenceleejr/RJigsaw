// Dear emacs, this is -*-c++-*-
/******************************************************************************
Name:        PileupReweightingAlg

Author:      Karsten Koeneke (CERN)
Created:     April 2011

Description: Algorithm to get the calculated MC pilup weight and attach it to the event as UserData.
******************************************************************************/

// Preprocessor magic for debugging
#define XXX std::cout << "I am here: " << __FILE__ << ":" << __LINE__ << std::endl;

// This class' header
#include "PileupReweighting/PileupReweightingAlg.h"

// Include the ROOT class from this package
#include "PileupReweighting/TPileupReweighting.h"


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

// Event information
#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"





//=============================================================================
// Constructor
//=============================================================================
PileupReweightingAlg::PileupReweightingAlg(const std::string& name, ISvcLocator* pSvcLocator) :
  AthAlgorithm(name, pSvcLocator),
  m_userDataSvc( "UserDataSvc", name ),
  m_tPileUp(0)
{
  // Properties go here
  declareProperty("userDataSvc",       m_userDataSvc,         "The UserData service" );
  declareProperty("userDataEventPileupWeightName",  m_userDataWeightName="pileupWeight",
                  "Name of the UserData event pileup weight variable" );

  declareProperty("EventInfoName",     m_eventInfoName="MyEvent", "The name of the EventInfo container" );

  declareProperty("dataROOTFileName",  m_dataRootFileName="", "The name of the data root file" );
  declareProperty("dataROOTHistName",  m_dataRootHistName="", "The name of the data histogram" );
  declareProperty("mcROOTFileName",    m_mcRootFileName="",   "The name of the MC root file" );
  declareProperty("mcROOTHistName",    m_mcRootHistName="",   "The name of the MC histogram" );
}



//=============================================================================
// Destructor
//=============================================================================
PileupReweightingAlg::~PileupReweightingAlg()
{
}

StatusCode PileupReweightingAlg::initialize()
{
  // Print the used configuration
  if (msgLvl(MSG::INFO))
    {
      msg(MSG::INFO) 
        << "==> initialize " << name() << "..."
        << endreq;
      
      // Print out the used configuration
      msg(MSG::INFO) << " using userDataSvc                   = " << m_userDataSvc << endreq;
      msg(MSG::INFO) << " using userDataEventPileupWeightName = " << m_userDataWeightName << endreq;

      msg(MSG::INFO) << " using EventInfoName                 = " << m_eventInfoName << endreq;

      msg(MSG::INFO) << " using dataROOTFileName              = " << m_dataRootFileName << endreq;
      msg(MSG::INFO) << " using dataROOTHistName              = " << m_dataRootHistName << endreq;
      msg(MSG::INFO) << " using mcROOTFileName                = " << m_mcRootFileName << endreq;
      msg(MSG::INFO) << " using mcROOTHistName                = " << m_mcRootHistName << endreq;
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
  m_tPileUp = new Root::TPileupReweighting( (((std::string)"T")+name()).c_str() );
  int isGood = m_tPileUp->initialize( m_dataRootFileName,
                                      m_dataRootHistName,
                                      m_mcRootFileName,
                                      m_mcRootHistName );
  if ( isGood == 0 )
    {
      REPORT_MESSAGE(MSG::DEBUG) << "Initialization of TPileupReweighting successful";
    }
  else if ( isGood == -1 )
    {
      REPORT_MESSAGE(MSG::ERROR) << "Could not load the data ROOT file with name " << m_dataRootFileName;
      sc = StatusCode::FAILURE;
    }
  else if ( isGood == -2 )
    {
      REPORT_MESSAGE(MSG::ERROR) << "Could not load the data ROOT histogram with name " << m_dataRootHistName;
      sc = StatusCode::FAILURE;
    }
  else if ( isGood == -3 )
    {
      REPORT_MESSAGE(MSG::ERROR) << "Could not load the MC ROOT file with name " << m_mcRootFileName;
      sc = StatusCode::FAILURE;
    }
  else if ( isGood == -4 )
    {
      REPORT_MESSAGE(MSG::ERROR) << "Could not load the MC ROOT histogram with name " << m_mcRootHistName;
      sc = StatusCode::FAILURE;
    }
  else if ( isGood == -15 )
    {
      REPORT_MESSAGE(MSG::ERROR) << "Could not normalize the data ROOT histogram with name " << m_dataRootHistName << " ...invalid pointer";
      sc = StatusCode::FAILURE;
    }
  else if ( isGood == -25 )
    {
      REPORT_MESSAGE(MSG::ERROR) << "Could not normalize the data ROOT histogram with name " << m_dataRootHistName << " ...empty histogram";
      sc = StatusCode::FAILURE;
    }
  else if ( isGood == -16 )
    {
      REPORT_MESSAGE(MSG::ERROR) << "Could not normalize the MC ROOT histogram with name " << m_mcRootHistName << " ...invalid pointer";
      sc = StatusCode::FAILURE;
    }
  else if ( isGood == -26 )
    {
      REPORT_MESSAGE(MSG::ERROR) << "Could not normalize the MC ROOT histogram with name " << m_mcRootHistName << " ...empty histogram";
      sc = StatusCode::FAILURE;
    }
  else if ( isGood == -10 )
    {
      REPORT_MESSAGE(MSG::ERROR) << "Found that the MC histogram is integer-binned with bin edges at integer numbers!"
                                 << " But in that case, the bin boundaries have to be at half-integer values!!! Fix it!";
      sc = StatusCode::FAILURE;
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
StatusCode PileupReweightingAlg::execute() 
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


  // Declare the default event weight
  float pileupEventWeight(-1.0);


  // retrive the EventInfo container
  const EventInfo* eventInfo;
  sc = evtStore()->retrieve(eventInfo, m_eventInfoName);
  if ( sc.isFailure() || !eventInfo ) 
    {
      if ( m_nEventsProcessed < 5 )
        {
          REPORT_MESSAGE(MSG::INFO) << "Container '" << m_eventInfoName 
                                    << "' could not be retrieved from StoreGate !" 
                                    << " Assuming that a data file is loaded. Nothing will be done!";
        }
      else
        {
          REPORT_MESSAGE(MSG::DEBUG) << "Container '" << m_eventInfoName 
                                     << "' could not be retrieved from StoreGate !" 
                                     << " Assuming that a data file is loaded. Nothing will be done!";
        }
      sc = StatusCode::SUCCESS;
      return sc;
    } 
  else
    {
      REPORT_MESSAGE(MSG::DEBUG) << "Container '" << m_eventInfoName 
                                 << "' retrieved from StoreGate";
    }


  // Only do something when we have a valid EventInfo object
  if ( eventInfo )
    {
      // Get the mu value for this event
      float mu = (float)eventInfo->event_ID()->lumi_block();

      // Get the pileup weight for this event
      pileupEventWeight = m_tPileUp->getPileupWeight( mu );
      

      //---------------------------------------
      // Now, save the result to UserData
      //---------------------------------------
      CHECK( m_userDataSvc->decorateEvent( m_userDataWeightName, pileupEventWeight ) );

    } // End: if ( eventInfo )

  
  return sc;
}




//=============================================================================
// Athena finalize method
//=============================================================================
StatusCode PileupReweightingAlg::finalize() 
{
  // Declare the simple StatusCode
  StatusCode sc(StatusCode::SUCCESS);

  // Print info messages
  if (msgLvl(MSG::INFO))
    {
      msg(MSG::INFO) 
        << "==> finalize " << name() << "...\n"
        << "***************************************************************\n"
        << "Results of " << name() << " pileup reweighting algorithm:\n"
        << "-------------"
        << endreq;
      msg(MSG::INFO) 
        << " Number of processed events:  " << m_nEventsProcessed
        << endreq;
    }

  // Delete the pointers
  if ( m_tPileUp ) delete m_tPileUp;

  return sc;
}

