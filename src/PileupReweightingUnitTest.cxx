#include "PileupReweighting/PileupReweightingUnitTest.h"
#include "PileupReweighting/TPileupReweighting.h"
#include "GaudiKernel/ServiceHandle.h"

#include "SGTools/BuiltinsClids.h" 
#include "SGTools/StlVectorClids.h"

///Unit test to make a prw file and then use it

PileupReweightingUnitTest::PileupReweightingUnitTest(const std::string& name, ISvcLocator* pSvcLocator) :
  AthAlgorithm(name, pSvcLocator),
  m_pileupTool(0),m_pileupTool2(0)
{
   declareProperty("Mode",m_mode=0);
}

PileupReweightingUnitTest::~PileupReweightingUnitTest() {
   delete m_pileupTool;delete m_pileupTool2;
}

StatusCode PileupReweightingUnitTest::initialize() {

   m_pileupTool = new Root::TPileupReweighting;m_pileupTool2 = new Root::TPileupReweighting;
   if(m_mode==0) {
      //make the prw file
      m_pileupTool->UsePeriodConfig("MC11b");
   } else {
      //use the prw file
      m_pileupTool->AddConfigFile("/tmp/myfile.prw.root");
      m_pileupTool->AddLumiCalcFile("unittest.lumicalc.root");
      m_pileupTool->SetUnrepresentedDataAction(2);

      m_pileupTool2->AddConfigFile("/tmp/myfile.prw.root");
      m_pileupTool2->AddLumiCalcFile("unittest.lumicalc.root");
      m_pileupTool2->SetUnrepresentedDataAction(1);

      m_hist = new TH1D("hist1","",50,-0.5,49.5);m_hist2 = new TH1D("hist2","",50,-0.5,49.5);
      m_dataHist = new TH1D("histd","",50,-0.5,49.5);
      TFile f1("unittest.lumicalc.root");
      TH1* h = (TH1*)f1.Get("avgintperbx");
      for(int i=0;i<h->GetNbinsX();i++) {
         m_dataHist->Fill(h->GetBinCenter(i),h->GetBinContent(i)/1E6);
      }
      f1.Close();
      m_pileupTool2->Initialize();
   }
   m_pileupTool->Initialize();

   return StatusCode::SUCCESS;
}


StatusCode PileupReweightingUnitTest::execute() {
   //CamEvent event;
   const unsigned int* runNumber = 0;
   if(evtStore()->retrieve(runNumber, "RunNumber").isFailure()) { ATH_MSG_ERROR("no runnumber"); }
   const unsigned int* chanNumber = 0;
   if(evtStore()->retrieve(chanNumber, "mc_channel_number").isFailure()) { ATH_MSG_ERROR("no channel number"); }
   const float* avgintperbx = 0;
   if(evtStore()->retrieve(avgintperbx, "averageIntPerXing").isFailure()) { ATH_MSG_ERROR("no mu"); }
   const std::vector<std::vector<double> >* weights = 0;
   if(evtStore()->retrieve(weights, "mcevt_weight").isFailure()) { ATH_MSG_ERROR("no weight"); }


   if(m_mode==0) {
      m_pileupTool->Fill(*runNumber,*chanNumber,weights->at(0)[0],*avgintperbx);
   } else {
      double w1 = m_pileupTool->GetCombinedWeight(*runNumber,*chanNumber,*avgintperbx);
      double w2 = m_pileupTool2->GetCombinedWeight(*runNumber,*chanNumber,*avgintperbx);

      m_hist->Fill(*avgintperbx,w1*(weights->at(0)[0]));
      m_hist2->Fill(*avgintperbx,w2*(weights->at(0)[0]));

   }
   return StatusCode::SUCCESS;
}

StatusCode PileupReweightingUnitTest::finalize() {

   if(m_mode==0) m_pileupTool->WriteToFile("/tmp/myfile.prw.root");
   else {
      //scale the histograms by Lumi/sumOfPileupWeights
      m_hist->Scale(m_pileupTool->GetIntegratedLumi()/m_hist->GetSumOfWeights());
      m_hist2->Scale(m_pileupTool2->GetIntegratedLumi()/m_hist2->GetSumOfWeights());
      m_dataHist->Scale(m_pileupTool->GetIntegratedLumi()/m_dataHist->GetSumOfWeights());

      //iterate over bins and print the deviations bigger than 0.0000001
      ATH_MSG_INFO("hist 1");
      for(int i=0;i<m_hist->GetNbinsX();i++) {
         if(fabs(m_hist->GetBinContent(i)-m_dataHist->GetBinContent(i))>0.000001) {
            ATH_MSG_INFO("Difference @ mu=" << m_hist->GetBinCenter(i) << "  " << m_hist->GetBinContent(i) << " vs " << m_dataHist->GetBinContent(i));
         }
      }

      m_dataHist->Scale(m_pileupTool2->GetIntegratedLumi()/m_dataHist->GetSumOfWeights());

      ATH_MSG_INFO("hist 2");
      for(int i=0;i<m_hist2->GetNbinsX();i++) {
         if(fabs(m_hist2->GetBinContent(i)-m_dataHist->GetBinContent(i))>0.000001) {
            ATH_MSG_INFO("Difference @ mu=" << m_hist2->GetBinCenter(i) << "  " << m_hist2->GetBinContent(i) << " vs " << m_dataHist->GetBinContent(i));
         }
      }

      TFile f2("unittest.results.root","RECREATE");
      m_hist->Write("unrep2_mc");
      m_hist2->Write("unrep1_mc");
      m_dataHist->Write("data");
      f2.Close();

   }

   return StatusCode::SUCCESS;
}

