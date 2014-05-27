#include "RJigsaw/RJigsawUnitTest.h"
#include "RJigsaw/TRJigsaw.h"
#include "GaudiKernel/ServiceHandle.h"

#include "SGTools/BuiltinsClids.h" 
#include "SGTools/StlVectorClids.h"

///Unit test to make a prw file and then use it

RJigsawUnitTest::RJigsawUnitTest(const std::string& name, ISvcLocator* pSvcLocator) :
  AthAlgorithm(name, pSvcLocator),
  m_RJigsawTool(0),m_RJigsawTool2(0)
{
   declareProperty("Mode",m_mode=0);
}

RJigsawUnitTest::~RJigsawUnitTest() {
   delete m_RJigsawTool;delete m_RJigsawTool2;
}

StatusCode RJigsawUnitTest::initialize() {

   m_RJigsawTool = new Root::TRJigsaw;m_RJigsawTool2 = new Root::TRJigsaw;
   m_RJigsawTool->Initialize();

   return StatusCode::SUCCESS;
}


StatusCode RJigsawUnitTest::execute() {
   //CamEvent event;
   const unsigned int* runNumber = 0;
   if(evtStore()->retrieve(runNumber, "RunNumber").isFailure()) { ATH_MSG_ERROR("no runnumber"); }
   const unsigned int* chanNumber = 0;
   if(evtStore()->retrieve(chanNumber, "mc_channel_number").isFailure()) { ATH_MSG_ERROR("no channel number"); }
   const float* avgintperbx = 0;
   if(evtStore()->retrieve(avgintperbx, "averageIntPerXing").isFailure()) { ATH_MSG_ERROR("no mu"); }
   const std::vector<std::vector<double> >* weights = 0;
   if(evtStore()->retrieve(weights, "mcevt_weight").isFailure()) { ATH_MSG_ERROR("no weight"); }

   return StatusCode::SUCCESS;
}

StatusCode RJigsawUnitTest::finalize() {



   return StatusCode::SUCCESS;
}

