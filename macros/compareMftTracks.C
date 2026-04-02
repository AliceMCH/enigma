#if !defined(__CLING__) || defined(__ROOTCLING__)
/*
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <filesystem>
#include <vector>

#include <Rtypes.h>
#include <ROOT/RDataFrame.hxx>
#include <TChain.h>

#include <fairlogger/Logger.h>
#include "Framework/Logger.h"

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

#include "MFTAlignment/TracksToRecords.h"
*/
#endif

//_________________________________
// alienv setenv O2Physics/latest -c root -l
// .L ~/cernbox/alice/enigma/macros/runTracksToRecords.C++
// runTracksToRecords()
//_________________________________
void compareMftTracks(std::string treeFileName1 = "mfttracks.root",
                      std::string treeFileName2 = "mfttracks-realign.root",
                      const int minClusters = 5)
{
  //return;
  TFile treeFile1(treeFileName1.data());
  TFile treeFile2(treeFileName2.data());

  TTree* tracksTrees[2] = { (TTree*)treeFile1.Get("o2sim"), (TTree*)treeFile2.Get("o2sim") };
  std::cout << "Tracks trees: " << tracksTrees[0] << ", " << tracksTrees[1] << std::endl;
  std::array<std::vector<o2::mft::TrackMFT>*, 2> mftTracks = { nullptr, nullptr };

  std::array<TH1F*, 3> hTrackChi2 = {
      new TH1F("mftTrackChi2Ref", "MFT track #chi^{2} (ref)", 100, 0, 100),
      new TH1F("mftTrackChi2", "MFT track #chi^{2}", 100, 0, 100),
      new TH1F("mftTrackChi2Ratio", "MFT track #chi^{2} (ratio", 100, 0, 100)
  };
  std::array<TH1F*, 3> hTrackNClus = {
      new TH1F("mftTrackNClusRef", "MFT track # of clusters (ref)", 10, 0.5, 10.5),
      new TH1F("mftTrackNClus", "MFT track # of clusters", 10, 0.5, 10.5),
      new TH1F("mftTrackNClusRatio", "MFT track # of clusters (ratio)", 10, 0.5, 10.5)
  };
  std::array<TH2F*, 3> hTrackXY = {
      new TH2F("mftTrackXYRef", "MFT track Y vs. X (ref)", 150, -15, 15, 150, -15, 15),
      new TH2F("mftTrackXY", "MFT track Y vs. X", 150, -15, 15, 150, -15, 15),
      new TH2F("mftTrackXYRatio", "MFT track Y vs. X (ratio)", 150, -15, 15, 150, -15, 15)
  };
  std::array<TH1F*, 3> hPhiMFT = {
      new TH1F("phiMFTRef", "#phi pf MFT track (ref);#phi (rad)", 100, -TMath::Pi(), TMath::Pi()),
      new TH1F("phiMFT", "#phi pf MFT track;#phi (rad)", 100, -TMath::Pi(), TMath::Pi()),
      new TH1F("phiMFTRatio", "#phi pf MFT track (ratio);#phi (rad)", 100, -TMath::Pi(), TMath::Pi())
  };

  for (int i = 0; i < 2; i++) {
    tracksTrees[i]->SetBranchAddress("MFTTrack", &(mftTracks[i]));
    std::cout << "  tracksTrees #" << i << ": " << mftTracks[i] << std::endl;

    int nTracks = 0;
    for (int j = 0; j < tracksTrees[i]->GetEntries(); j++) {
      tracksTrees[i]->GetEntry(j);
      for (const auto& mftTrack : (*(mftTracks[i]))) {
        hTrackXY[i]->Fill(mftTrack.getX(), mftTrack.getY());
        hTrackChi2[i]->Fill(mftTrack.getChi2OverNDF());
        hPhiMFT[i]->Fill(mftTrack.getPhi());
        hTrackNClus[i]->Fill(mftTrack.getNumberOfPoints());

        if (mftTrack.getNumberOfPoints() < minClusters) {
          continue;
        }
        nTracks += 1;
      }
    }
    std::cout << nTracks << std::endl;
  }

  TCanvas c("c", "c", 1200, 800);

  hTrackXY[0]->Draw("colz");
  c.SaveAs("mftTracks.pdf(");
  hTrackXY[1]->Draw("colz");
  c.SaveAs("mftTracks.pdf");
  hTrackXY[2]->Add(hTrackXY[1]);
  hTrackXY[2]->Divide(hTrackXY[0]);
  hTrackXY[2]->SetMinimum(0);
  hTrackXY[2]->SetMaximum(2);
  hTrackXY[2]->Draw("colz");
  c.SaveAs("mftTracks.pdf");

  hTrackChi2[0]->Draw();
  c.SaveAs("mftTracks.pdf");
  hTrackChi2[1]->Draw();
  c.SaveAs("mftTracks.pdf");
  hTrackChi2[2]->Add(hTrackChi2[1]);
  hTrackChi2[2]->Divide(hTrackChi2[0]);
  hTrackChi2[2]->SetMinimum(0);
  hTrackChi2[2]->SetMaximum(2);
  hTrackChi2[2]->Draw("colz");
  c.SaveAs("mftTracks.pdf");

  hTrackNClus[0]->Draw();
  c.SaveAs("mftTracks.pdf");
  hTrackNClus[1]->Draw();
  c.SaveAs("mftTracks.pdf");
  hTrackNClus[2]->Add(hTrackNClus[1]);
  hTrackNClus[2]->Divide(hTrackNClus[0]);
  hTrackNClus[2]->SetMinimum(0);
  hTrackNClus[2]->SetMaximum(2);
  hTrackNClus[2]->Draw("colz");
  c.SaveAs("mftTracks.pdf");

  hPhiMFT[0]->Draw();
  c.SaveAs("mftTracks.pdf");
  hPhiMFT[1]->Draw();
  c.SaveAs("mftTracks.pdf");
  hPhiMFT[2]->Add(hPhiMFT[1]);
  hPhiMFT[2]->Divide(hPhiMFT[0]);
  hPhiMFT[2]->SetMinimum(0);
  hPhiMFT[2]->SetMaximum(2);
  hPhiMFT[2]->Draw("colz");
  c.SaveAs("mftTracks.pdf)");
}
