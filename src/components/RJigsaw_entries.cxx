#include "GaudiKernel/DeclareFactoryEntries.h"

#include "RJigsaw/RJigsawAlg.h"
#include "RJigsaw/RJigsawUnitTest.h"

//#include "RJigsaw/FastLumiSvc.h"

DECLARE_ALGORITHM_FACTORY( RJigsawAlg )
DECLARE_ALGORITHM_FACTORY( RJigsawUnitTest )

//DECLARE_SERVICE_FACTORY(FastLumiSvc)

DECLARE_FACTORY_ENTRIES( RJigsaw ) 
{
  DECLARE_ALGORITHM( RJigsawAlg );
  DECLARE_ALGORITHM( RJigsawUnitTest );
  //DECLARE_SERVICE(FastLumiSvc);
}


