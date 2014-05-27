#include "GaudiKernel/DeclareFactoryEntries.h"

#include "PileupReweighting/PileupReweightingAlg.h"
#include "PileupReweighting/PileupReweightingUnitTest.h"

//#include "PileupReweighting/FastLumiSvc.h"

DECLARE_ALGORITHM_FACTORY( PileupReweightingAlg )
DECLARE_ALGORITHM_FACTORY( PileupReweightingUnitTest )

//DECLARE_SERVICE_FACTORY(FastLumiSvc)

DECLARE_FACTORY_ENTRIES( PileupReweighting ) 
{
  DECLARE_ALGORITHM( PileupReweightingAlg );
  DECLARE_ALGORITHM( PileupReweightingUnitTest );
  //DECLARE_SERVICE(FastLumiSvc);
}


