// Dear emacs, this is -*-c++-*-

#ifndef PILEUPREWEIGHTING_UNITTEST_H
#define PILEUPREWEIGHTING_UNITTEST_H


#include "AthenaBaseComps/AthAlgorithm.h"

#include "TH1D.h"

namespace Root{
  class TPileupReweighting;
}

class PileupReweightingUnitTest : public AthAlgorithm {

public: 
  /** Standard constructor */
  PileupReweightingUnitTest(const std::string& name, ISvcLocator* pSvcLocator);

  /** Standard destructor */
  virtual ~PileupReweightingUnitTest();

public:
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

private:
   bool m_mode;

  Root::TPileupReweighting* m_pileupTool;
   Root::TPileupReweighting* m_pileupTool2;

   TH1D* m_hist; TH1D* m_hist2; TH1D* m_dataHist;
   TH1D* m_hist_w; TH1D* m_hist2_w; 
};

#endif
