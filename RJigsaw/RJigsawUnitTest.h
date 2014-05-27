// Dear emacs, this is -*-c++-*-

#ifndef RJIGSAW_UNITTEST_H
#define RJIGSAW_UNITTEST_H


#include "AthenaBaseComps/AthAlgorithm.h"

#include "TH1D.h"

namespace Root{
  class TRJigsaw;
}

class RJigsawUnitTest : public AthAlgorithm {

public: 
  /** Standard constructor */
  RJigsawUnitTest(const std::string& name, ISvcLocator* pSvcLocator);

  /** Standard destructor */
  virtual ~RJigsawUnitTest();

public:
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

private:
   bool m_mode;

  Root::TRJigsaw* m_pileupTool;
   Root::TRJigsaw* m_pileupTool2;

   TH1D* m_hist; TH1D* m_hist2; TH1D* m_dataHist;
   TH1D* m_hist_w; TH1D* m_hist2_w; 
};

#endif
