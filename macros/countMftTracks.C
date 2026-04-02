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
void countMftTracks(std::string treeFileName = "mfttracks.root",
                    const int minClusters = 5)
{
  //return;
  TFile treeFile(treeFileName.data());
  TTree* tracksTree = (TTree*)treeFile.Get("o2sim");
  std::vector<o2::mft::TrackMFT>* mftTracks = nullptr;
  tracksTree->SetBranchAddress("MFTTrack", &mftTracks);
  std::vector<o2::itsmft::ROFRecord>* mftTracksROF = nullptr;
  tracksTree->SetBranchAddress("MFTTracksROF", &mftTracksROF);

  int nTracks = 0;
  for (int i = 0; i < tracksTree->GetEntries(); i++) {
    tracksTree->GetEntry(i);

    //std::cout << "Entry #" << i << " mftTracks.size(): " << mftTracks->size() << std::endl;

    //______________________________________________________
    for (int iMftRof = 0; iMftRof < mftTracksROF->size(); iMftRof++) {

      const auto& mftRofRec = (*mftTracksROF)[iMftRof];
      //std::cout << "mftRofRec.getNEntries(): " << mftRofRec.getNEntries() << std::endl;

      for (int itrk = 0; itrk < mftRofRec.getNEntries(); itrk++) {
        int trkEntry = mftRofRec.getFirstEntry() + itrk; // entry of icl-th cluster of this ROF in the vector of clusters
        const auto& mftTrack = (*mftTracks)[trkEntry];

        //for (const auto& mftTrack : (*mftTracks)) {

        //std::cout << std::format("Entry #{} MFT track nClus={} ROF: orbit={} bc={}", i, mftTrack.getNumberOfPoints(), mftRofRec.getBCData().orbit,  mftRofRec.getBCData().bc) << std::endl;
        //std::cout << std::format("    x={:+0.3f} y={:+0.3f} z={:+0.3f} phi={:+0.3f} tanl={:+0.3f} ",
        //    mftTrack.getX(), mftTrack.getY(), mftTrack.getZ(), mftTrack.getPhi(), mftTrack.getTanl()) << std::endl;

        if (mftTrack.getNumberOfPoints() < minClusters) {
          continue;
        }
        nTracks += 1;
      }
    }
  }
  std::cout << nTracks << std::endl;
}
