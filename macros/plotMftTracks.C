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
void plotMftTracks(std::string treeFileName = "mfttracks.root",
                    const int minClusters = 5)
{
  //return;
  TFile treeFile(treeFileName.data());
  TTree* tracksTree = (TTree*)treeFile.Get("o2sim");
  std::vector<o2::mft::TrackMFT>* mftTracks = nullptr;
  tracksTree->SetBranchAddress("MFTTrack", &mftTracks);
  std::vector<o2::itsmft::ROFRecord>* mftTracksROF = nullptr;
  tracksTree->SetBranchAddress("MFTTracksROF", &mftTracksROF);

  TH1F hTrackNClus("mftTrackNClus", "MFT track # of clusters", 10, 0.5, 10.5);
  TH2F hTrackXY("mftTrackXY", "MFT track Y vs. X", 150, -15, 15, 150, -15, 15);
  TH1F hPhiMFT("phiMFT", "#phi pf MFT track;#phi (rad)", 100, -TMath::Pi(), TMath::Pi());

  int nTracks = 0;
  for (int i = 0; i < tracksTree->GetEntries(); i++) {
    tracksTree->GetEntry(i);

    //std::cout << "Entry #" << i << " mftTracks.size(): " << mftTracks->size() << std::endl;

    for (const auto& mftTrack : (*mftTracks)) {

    //std::cout << "Entry #" << i << " mftTracksROF.size(): " << mftTracksROF->size() << std::endl;

    //______________________________________________________
    //for (int iMftRof = 0; iMftRof < mftTracksROF->size(); iMftRof++) {

      //const auto& mftRofRec = (*mftTracksROF)[iMftRof];
      //std::cout << "mftRofRec.getNEntries(): " << mftRofRec.getNEntries() << std::endl;

      //for (int itrk = 0; itrk < mftRofRec.getNEntries(); itrk++) {
        //int trkEntry = mftRofRec.getFirstEntry() + itrk; // entry of icl-th cluster of this ROF in the vector of clusters
        //const auto& mftTrack = (*mftTracks)[trkEntry];

        hTrackXY.Fill(mftTrack.getX(), mftTrack.getY());
        hPhiMFT.Fill(mftTrack.getPhi());
        hTrackNClus.Fill(mftTrack.getNumberOfPoints());

        //for (const auto& mftTrack : (*mftTracks)) {

        //std::cout << std::format("Entry #{} MFT track nClus={} ROF: orbit={} bc={}", i, mftTrack.getNumberOfPoints(), mftRofRec.getBCData().orbit,  mftRofRec.getBCData().bc) << std::endl;
        //std::cout << std::format("    x={:+0.3f} y={:+0.3f} z={:+0.3f} phi={:+0.3f} tanl={:+0.3f} ",
        //    mftTrack.getX(), mftTrack.getY(), mftTrack.getZ(), mftTrack.getPhi(), mftTrack.getTanl()) << std::endl;

        if (mftTrack.getNumberOfPoints() < minClusters) {
          continue;
        }
        nTracks += 1;
      //}
    }
  }
  std::cout << nTracks << std::endl;

  TCanvas c("c", "c", 1200, 800);

  hTrackXY.Draw("colz");
  c.SaveAs("mftTracks.pdf(");

  hTrackNClus.Draw();
  c.SaveAs("mftTracks.pdf");

  hPhiMFT.Draw();
  c.SaveAs("mftTracks.pdf)");
}
