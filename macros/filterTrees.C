#if !defined(__CLING__) || defined(__ROOTCLING__)

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

#endif

namespace fs = std::filesystem;

//_________________________________
// alienv setenv O2Physics/latest -c root -l
// .L ~/cernbox/alice/enigma/macros/runTracksToRecords.C++
// runTracksToRecords()
//_________________________________
void _filterTrees(std::string mftTracksFileName = "mfttracks.root",
                  std::string mftClustersFileName = "mftclusters.root",
                  std::string itsClustersFileName = "o2clus_its.root",
                  const int minClusters = 5)
{
  // MFT tracks
  TFile mftTracksFile(mftTracksFileName.data());
  TTree* mftTracksTree = (TTree*)mftTracksFile.Get("o2sim");
  std::vector<o2::mft::TrackMFT>* mftTracks = nullptr;
  std::vector<o2::mft::TrackMFT> mftTracksFiltered;
  mftTracksTree->SetBranchAddress("MFTTrack", &mftTracks);
  std::vector<o2::itsmft::ROFRecord>* mftTracksROF = nullptr;
  std::vector<o2::itsmft::ROFRecord> mftTracksROFFiltered;
  mftTracksTree->SetBranchAddress("MFTTracksROF", &mftTracksROF);
  std::vector<int>* mftTrackClusIdx = nullptr;
  std::vector<int> mftTrackClusIdxFiltered;
  mftTracksTree->SetBranchAddress("MFTTrackClusIdx", &mftTrackClusIdx);

  // MFT clusters
  TFile mftClustersFile(mftClustersFileName.data());
  TTree* mftClustersTree = (TTree*)mftClustersFile.Get("o2sim");
  std::vector<o2::itsmft::CompClusterExt>* mftClusters = nullptr;
  std::vector<o2::itsmft::CompClusterExt> mftClustersFiltered;
  mftClustersTree->SetBranchAddress("MFTClusterComp", &mftClusters);
  std::vector<o2::itsmft::ROFRecord>* mftClustersROF = nullptr;
  std::vector<o2::itsmft::ROFRecord> mftClustersROFFiltered;
  mftClustersTree->SetBranchAddress("MFTClustersROF", &mftClustersROF);
  std::vector<unsigned char>* mftClusterPatterns = nullptr;
  std::vector<unsigned char> mftClusterPatternsFiltered;
  mftClustersTree->SetBranchAddress("MFTClusterPatt", &mftClusterPatterns);

  // ITS clusters
  TFile itsClustersFile(itsClustersFileName.data());
  TTree* itsClustersTree = (TTree*)itsClustersFile.Get("o2sim");
  std::vector<o2::itsmft::CompClusterExt>* itsClusters = nullptr;
  std::vector<o2::itsmft::CompClusterExt> itsClustersFiltered;
  itsClustersTree->SetBranchAddress("ITSClusterComp", &itsClusters);
  std::vector<o2::itsmft::ROFRecord>* itsClustersROF = nullptr;
  std::vector<o2::itsmft::ROFRecord> itsClustersROFFiltered;
  itsClustersTree->SetBranchAddress("ITSClustersROF", &itsClustersROF);
  std::vector<unsigned char>* itsPatterns = nullptr;
  std::vector<unsigned char> itsPatternsFiltered;
  itsClustersTree->SetBranchAddress("ITSClusterPatt", &itsPatterns);

  // Filtered outputs
  TFile* mftTracksFileFiltered = nullptr;
  TTree* mftTracksTreeFiltered = nullptr;

  TFile* mftClustersFileFiltered = nullptr;
  TTree* mftClustersTreeFiltered = nullptr;

  TFile* itsClustersFileFiltered = nullptr;
  TTree* itsClustersTreeFiltered = nullptr;


  int nTracks = 0;
  for (int i = 0; i < mftTracksTree->GetEntries(); i++) {
    mftTracksTree->GetEntry(i);
    mftClustersTree->GetEntry(i);
    itsClustersTree->GetEntry(i);

    for (const auto& mftTrack : (*mftTracks)) {

      if (mftTrack.getNumberOfPoints() < minClusters) {
        continue;
      }

      mftTracksFiltered = *mftTracks;
      mftTracksROFFiltered = *mftTracksROF;
      mftTrackClusIdxFiltered = *mftTrackClusIdx;

      mftClustersFiltered = *mftClusters;
      mftClustersROFFiltered = *mftClustersROF;
      mftClusterPatternsFiltered = *mftClusterPatterns;

      itsClustersFiltered = *itsClusters;
      itsClustersROFFiltered = *itsClustersROF;
      itsPatternsFiltered = *itsPatterns;

      if (!mftTracksFileFiltered) {
        mftTracksFileFiltered = new TFile("mfttracks-filtered.root", "RECREATE");
        mftClustersFileFiltered = new TFile("mftclusters-filtered.root", "RECREATE");
        itsClustersFileFiltered = new TFile("o2clus_its-filtered.root", "RECREATE");

        mftTracksFileFiltered->cd();
        mftTracksTreeFiltered = new TTree("o2sim", "MFT tracks");
        mftTracksTreeFiltered->Branch("MFTTrack", &mftTracksFiltered);
        mftTracksTreeFiltered->Branch("MFTTracksROF", &mftTracksROFFiltered);
        mftTracksTreeFiltered->Branch("MFTTrackClusIdx", &mftTrackClusIdxFiltered);

        mftClustersFileFiltered->cd();
        mftClustersTreeFiltered = new TTree("o2sim", "MFT clusters");
        mftClustersTreeFiltered->Branch("MFTClusterComp", &mftClustersFiltered);
        mftClustersTreeFiltered->Branch("MFTClustersROF", &mftClustersROFFiltered);
        mftClustersTreeFiltered->Branch("MFTClusterPatt", &mftClusterPatternsFiltered);

        itsClustersFileFiltered->cd();
        itsClustersTreeFiltered = new TTree("o2sim", "ITS clusters");
        itsClustersTreeFiltered->Branch("ITSClusterComp", &itsClustersFiltered);
        itsClustersTreeFiltered->Branch("ITSClustersROF", &itsClustersROFFiltered);
        itsClustersTreeFiltered->Branch("ITSClusterPatt", &itsPatternsFiltered);
      }

      if (mftTracksTreeFiltered) {
        mftTracksTreeFiltered->Fill();
      }
      if (mftClustersFileFiltered) {
        mftClustersTreeFiltered->Fill();
      }
      if (itsClustersFileFiltered) {
        itsClustersTreeFiltered->Fill();
      }

      nTracks += 1;
    }
  }

  if (mftTracksFileFiltered) {
    mftTracksFileFiltered->cd();
    mftTracksTreeFiltered->Write();
  }
  if (mftClustersFileFiltered) {
    mftClustersFileFiltered->cd();
    mftClustersTreeFiltered->Write();
  }
  if (itsClustersFileFiltered) {
    itsClustersFileFiltered->cd();
    itsClustersTreeFiltered->Write();
  }

  std::cout << nTracks << std::endl;
}

size_t countLines(std::string fileName)
{
  size_t numLines = 0;
  ifstream in(fileName);
  std::string unused;
  while ( std::getline(in, unused) )
     ++numLines;
  return numLines;
}

void filterTrees(std::string folderListFile = "folder-list.txt",
                 const int minClusters = 5)
{
  fs::remove("stop");

  std::ifstream folderList(folderListFile);
  std::string folder;
  std::string runFolderPrev, tenMinutesFolderPrev;
  auto nFolders = countLines(folderListFile);
  int iFolder = 0;
  while (std::getline(folderList, folder)) {
    std::string tempFolder = folder;
    std::string ctfFolder = fs::path(tempFolder).filename();

    tempFolder = fs::path(tempFolder).parent_path();
    std::string tenMinutesFolder = fs::path(tempFolder).filename();

    tempFolder = fs::path(tempFolder).parent_path();
    tempFolder = fs::path(tempFolder).parent_path();
    std::string runFolder = std::string("input/") + fs::path(tempFolder).filename().native();

    fs::create_directories(runFolder);

    std::string outFolder = std::format("{}/{}/{}", runFolder, tenMinutesFolder, ctfFolder);
    std::cout << std::format("{}/{} -> {}", iFolder, nFolders, outFolder) << std::endl;

    iFolder += 1;

    if (std::filesystem::exists("stop")) {
      break;;
    }

    if (std::filesystem::exists(Form("%s/filtering.done", outFolder.c_str()))) {
      if (!fs::exists(Form("%s/o2_tfidinfo.root", outFolder.c_str())) && fs::exists(Form("%s/mfttracks-filtered.root", outFolder.c_str()))) {
        std::string copyCommand = std::format("alien_cp alien://{}/o2_tfidinfo.root file://./{}/o2_tfidinfo.root", folder, outFolder);
        std::cout << copyCommand << std::endl;
        std::system(copyCommand.c_str());
      }
      continue;
    }

    if (runFolderPrev != runFolder || tenMinutesFolderPrev != tenMinutesFolder) {
      //std::string copyCommand = std::string("alien_cp -parent 1 alien://") + fs::path(folder).parent_path() + "/*/mfttracks.root file://./" + outFolder;
      std::string copyCommand = std::format("alien_cp -parent 1 alien://{}/*/mfttracks.root file://./{}", fs::path(folder).parent_path().native(), runFolder);
      //std::string copyCommand = std::string("alien_cp -parent 1 alien://") + outFolder;
      std::cout << copyCommand << std::endl;
      std::system(copyCommand.c_str());
      runFolderPrev = runFolder;
      tenMinutesFolderPrev = tenMinutesFolder;
    }

    if (std::filesystem::exists("stop")) {
      break;;
    }

    std::string mftTracksFileName = Form("%s/mfttracks.root", outFolder.c_str());

    if (!fs::exists(mftTracksFileName)) {
      std::cout << mftTracksFileName << " not found" << std::endl;
      continue;
    }

    // check if we have at least one MFT tracks
    TFile mftTracksFile(mftTracksFileName.data());
    TTree* mftTracksTree = (TTree*)mftTracksFile.Get("o2sim");

    if (!mftTracksTree) {
      std::cout << mftTracksFileName << " does not contain a valid TTree" << std::endl;
      continue;
    }

    std::vector<o2::mft::TrackMFT>* mftTracks = nullptr;
    mftTracksTree->SetBranchAddress("MFTTrack", &mftTracks);

    //std::cout << "Inspecting " << mftTracksFileName << std::endl;
    bool mftTrackFound = false;
    for (int i = 0; i < mftTracksTree->GetEntries(); i++) {
      mftTracksTree->GetEntry(i);

      if (mftTracks->size() > 0) {
        std::cout << "mftTracks.size(): " << mftTracks->size() << std::endl;
        mftTrackFound = true;
        break;
      }
    }

    if (mftTrackFound) {
      std::string copyCommand = std::format("alien_cp alien://{}/mftclusters.root file://./{}/mftclusters.root", folder, outFolder);
      std::cout << copyCommand << std::endl;
      std::system(copyCommand.c_str());

      if (!fs::exists(Form("%s/mftclusters.root", outFolder.c_str()))) {
        std::cout << Form("%s/mftclusters.root", outFolder.c_str()) << " not found" << std::endl;
        continue;
      }

      copyCommand = std::format("alien_cp alien://{}/o2clus_its.root file://./{}/o2clus_its.root", folder, outFolder);
      std::cout << copyCommand << std::endl;
      std::system(copyCommand.c_str());

      if (!fs::exists(Form("%s/o2clus_its.root", outFolder.c_str()))) {
        std::cout << Form("%s/o2clus_its.root", outFolder.c_str()) << " not found" << std::endl;
        continue;
      }

      copyCommand = std::format("alien_cp alien://{}/o2_tfidinfo.root file://./{}/o2_tfidinfo.root", folder, outFolder);
      std::cout << copyCommand << std::endl;
      std::system(copyCommand.c_str());

      if (!fs::exists(Form("%s/o2_tfidinfo.root", outFolder.c_str()))) {
        std::cout << Form("%s/o2_tfidinfo.root", outFolder.c_str()) << " not found" << std::endl;
        continue;
      }

      std::string mftClustersFileName = Form("%s/mftclusters.root", outFolder.c_str());
      std::string itsClustersFileName = Form("%s/o2clus_its.root", outFolder.c_str());

      // MFT tracks
      std::vector<o2::mft::TrackMFT> mftTracksFiltered;
      std::vector<o2::itsmft::ROFRecord>* mftTracksROF = nullptr;
      std::vector<o2::itsmft::ROFRecord> mftTracksROFFiltered;
      mftTracksTree->SetBranchAddress("MFTTracksROF", &mftTracksROF);
      std::vector<int>* mftTrackClusIdx = nullptr;
      std::vector<int> mftTrackClusIdxFiltered;
      mftTracksTree->SetBranchAddress("MFTTrackClusIdx", &mftTrackClusIdx);

      // MFT clusters
      TFile mftClustersFile(mftClustersFileName.data());
      TTree* mftClustersTree = (TTree*)mftClustersFile.Get("o2sim");
      std::vector<o2::itsmft::CompClusterExt>* mftClusters = nullptr;
      std::vector<o2::itsmft::CompClusterExt> mftClustersFiltered;
      mftClustersTree->SetBranchAddress("MFTClusterComp", &mftClusters);
      std::vector<o2::itsmft::ROFRecord>* mftClustersROF = nullptr;
      std::vector<o2::itsmft::ROFRecord> mftClustersROFFiltered;
      mftClustersTree->SetBranchAddress("MFTClustersROF", &mftClustersROF);
      std::vector<unsigned char>* mftClusterPatterns = nullptr;
      std::vector<unsigned char> mftClusterPatternsFiltered;
      mftClustersTree->SetBranchAddress("MFTClusterPatt", &mftClusterPatterns);

      // ITS clusters
      TFile itsClustersFile(itsClustersFileName.data());
      TTree* itsClustersTree = (TTree*)itsClustersFile.Get("o2sim");
      std::vector<o2::itsmft::CompClusterExt>* itsClusters = nullptr;
      std::vector<o2::itsmft::CompClusterExt> itsClustersFiltered;
      itsClustersTree->SetBranchAddress("ITSClusterComp", &itsClusters);
      std::vector<o2::itsmft::ROFRecord>* itsClustersROF = nullptr;
      std::vector<o2::itsmft::ROFRecord> itsClustersROFFiltered;
      itsClustersTree->SetBranchAddress("ITSClustersROF", &itsClustersROF);
      std::vector<unsigned char>* itsPatterns = nullptr;
      std::vector<unsigned char> itsPatternsFiltered;
      itsClustersTree->SetBranchAddress("ITSClusterPatt", &itsPatterns);

      // Filtered outputs
      TFile* mftTracksFileFiltered = nullptr;
      TTree* mftTracksTreeFiltered = nullptr;

      TFile* mftClustersFileFiltered = nullptr;
      TTree* mftClustersTreeFiltered = nullptr;

      TFile* itsClustersFileFiltered = nullptr;
      TTree* itsClustersTreeFiltered = nullptr;

      for (int i = 0; i < mftTracksTree->GetEntries(); i++) {
        mftTracksTree->GetEntry(i);
        mftClustersTree->GetEntry(i);
        itsClustersTree->GetEntry(i);

        for (const auto& mftTrack : (*mftTracks)) {

          if (mftTrack.getNumberOfPoints() < minClusters) {
            continue;
          }

          mftTracksFiltered = *mftTracks;
          mftTracksROFFiltered = *mftTracksROF;
          mftTrackClusIdxFiltered = *mftTrackClusIdx;

          mftClustersFiltered = *mftClusters;
          mftClustersROFFiltered = *mftClustersROF;
          mftClusterPatternsFiltered = *mftClusterPatterns;

          itsClustersFiltered = *itsClusters;
          itsClustersROFFiltered = *itsClustersROF;
          itsPatternsFiltered = *itsPatterns;

          if (!mftTracksFileFiltered) {
            mftTracksFileFiltered = new TFile(Form("%s/mfttracks-filtered.root", outFolder.c_str()), "RECREATE");
            mftClustersFileFiltered = new TFile(Form("%s/mftclusters-filtered.root", outFolder.c_str()), "RECREATE");
            itsClustersFileFiltered = new TFile(Form("%s/o2clus_its-filtered.root", outFolder.c_str()), "RECREATE");

            mftTracksFileFiltered->cd();
            mftTracksTreeFiltered = new TTree("o2sim", "MFT tracks");
            mftTracksTreeFiltered->Branch("MFTTrack", &mftTracksFiltered);
            mftTracksTreeFiltered->Branch("MFTTracksROF", &mftTracksROFFiltered);
            mftTracksTreeFiltered->Branch("MFTTrackClusIdx", &mftTrackClusIdxFiltered);

            mftClustersFileFiltered->cd();
            mftClustersTreeFiltered = new TTree("o2sim", "MFT clusters");
            mftClustersTreeFiltered->Branch("MFTClusterComp", &mftClustersFiltered);
            mftClustersTreeFiltered->Branch("MFTClustersROF", &mftClustersROFFiltered);
            mftClustersTreeFiltered->Branch("MFTClusterPatt", &mftClusterPatternsFiltered);

            itsClustersFileFiltered->cd();
            itsClustersTreeFiltered = new TTree("o2sim", "ITS clusters");
            itsClustersTreeFiltered->Branch("ITSClusterComp", &itsClustersFiltered);
            itsClustersTreeFiltered->Branch("ITSClustersROF", &itsClustersROFFiltered);
            itsClustersTreeFiltered->Branch("ITSClusterPatt", &itsPatternsFiltered);
          }

          if (mftTracksTreeFiltered) {
            mftTracksTreeFiltered->Fill();
          }
          if (mftClustersFileFiltered) {
            mftClustersTreeFiltered->Fill();
          }
          if (itsClustersFileFiltered) {
            itsClustersTreeFiltered->Fill();
          }
        }
      }

      if (mftTracksFileFiltered) {
        mftTracksFileFiltered->cd();
        mftTracksTreeFiltered->Write();
        mftTracksFileFiltered->Close();
        delete mftTracksFileFiltered;
      }
      if (mftClustersFileFiltered) {
        mftClustersFileFiltered->cd();
        mftClustersTreeFiltered->Write();
        mftClustersFileFiltered->Close();
        delete mftClustersFileFiltered;
      }
      if (itsClustersFileFiltered) {
        itsClustersFileFiltered->cd();
        itsClustersTreeFiltered->Write();
        itsClustersFileFiltered->Close();
        delete itsClustersFileFiltered;
      }
    }

    fs::remove(Form("%s/mfttracks.root", outFolder.c_str()));
    fs::remove(Form("%s/mftclusters.root", outFolder.c_str()));
    fs::remove(Form("%s/o2clus_its.root", outFolder.c_str()));

    std::ofstream { Form("%s/filtering.done", outFolder.c_str()) };

    //break;
  }

  std::cout << std::format("\n{}/{} folders processed", iFolder, nFolders) << std::endl;

  return;

  std::string mftTracksFileName = "mfttracks.root";
  std::string mftClustersFileName = "mftclusters.root";
  std::string itsClustersFileName = "o2clus_its.root";
  // MFT tracks

  TFile mftTracksFile(mftTracksFileName.data());
  TTree* mftTracksTree = (TTree*)mftTracksFile.Get("o2sim");
  std::vector<o2::mft::TrackMFT>* mftTracks = nullptr;
  std::vector<o2::mft::TrackMFT> mftTracksFiltered;
  mftTracksTree->SetBranchAddress("MFTTrack", &mftTracks);
  std::vector<o2::itsmft::ROFRecord>* mftTracksROF = nullptr;
  std::vector<o2::itsmft::ROFRecord> mftTracksROFFiltered;
  mftTracksTree->SetBranchAddress("MFTTracksROF", &mftTracksROF);
  std::vector<int>* mftTrackClusIdx = nullptr;
  std::vector<int> mftTrackClusIdxFiltered;
  mftTracksTree->SetBranchAddress("MFTTrackClusIdx", &mftTrackClusIdx);

  // MFT clusters
  TFile mftClustersFile(mftClustersFileName.data());
  TTree* mftClustersTree = (TTree*)mftClustersFile.Get("o2sim");
  std::vector<o2::itsmft::CompClusterExt>* mftClusters = nullptr;
  std::vector<o2::itsmft::CompClusterExt> mftClustersFiltered;
  mftClustersTree->SetBranchAddress("MFTClusterComp", &mftClusters);
  std::vector<o2::itsmft::ROFRecord>* mftClustersROF = nullptr;
  std::vector<o2::itsmft::ROFRecord> mftClustersROFFiltered;
  mftClustersTree->SetBranchAddress("MFTClustersROF", &mftClustersROF);
  std::vector<unsigned char>* mftClusterPatterns = nullptr;
  std::vector<unsigned char> mftClusterPatternsFiltered;
  mftClustersTree->SetBranchAddress("MFTClusterPatt", &mftClusterPatterns);

  // ITS clusters
  TFile itsClustersFile(itsClustersFileName.data());
  TTree* itsClustersTree = (TTree*)itsClustersFile.Get("o2sim");
  std::vector<o2::itsmft::CompClusterExt>* itsClusters = nullptr;
  std::vector<o2::itsmft::CompClusterExt> itsClustersFiltered;
  itsClustersTree->SetBranchAddress("ITSClusterComp", &itsClusters);
  std::vector<o2::itsmft::ROFRecord>* itsClustersROF = nullptr;
  std::vector<o2::itsmft::ROFRecord> itsClustersROFFiltered;
  itsClustersTree->SetBranchAddress("ITSClustersROF", &itsClustersROF);
  std::vector<unsigned char>* itsPatterns = nullptr;
  std::vector<unsigned char> itsPatternsFiltered;
  itsClustersTree->SetBranchAddress("ITSClusterPatt", &itsPatterns);

  // Filtered outputs
  TFile* mftTracksFileFiltered = nullptr;
  TTree* mftTracksTreeFiltered = nullptr;

  TFile* mftClustersFileFiltered = nullptr;
  TTree* mftClustersTreeFiltered = nullptr;

  TFile* itsClustersFileFiltered = nullptr;
  TTree* itsClustersTreeFiltered = nullptr;


  int nTracks = 0;
  for (int i = 0; i < mftTracksTree->GetEntries(); i++) {
    mftTracksTree->GetEntry(i);
    mftClustersTree->GetEntry(i);
    itsClustersTree->GetEntry(i);

    for (const auto& mftTrack : (*mftTracks)) {

      if (mftTrack.getNumberOfPoints() < minClusters) {
        continue;
      }

      mftTracksFiltered = *mftTracks;
      mftTracksROFFiltered = *mftTracksROF;
      mftTrackClusIdxFiltered = *mftTrackClusIdx;

      mftClustersFiltered = *mftClusters;
      mftClustersROFFiltered = *mftClustersROF;
      mftClusterPatternsFiltered = *mftClusterPatterns;

      itsClustersFiltered = *itsClusters;
      itsClustersROFFiltered = *itsClustersROF;
      itsPatternsFiltered = *itsPatterns;

      if (!mftTracksFileFiltered) {
        mftTracksFileFiltered = new TFile("mfttracks-filtered.root", "RECREATE");
        mftClustersFileFiltered = new TFile("mftclusters-filtered.root", "RECREATE");
        itsClustersFileFiltered = new TFile("o2clus_its-filtered.root", "RECREATE");

        mftTracksFileFiltered->cd();
        mftTracksTreeFiltered = new TTree("o2sim", "MFT tracks");
        mftTracksTreeFiltered->Branch("MFTTrack", &mftTracksFiltered);
        mftTracksTreeFiltered->Branch("MFTTracksROF", &mftTracksROFFiltered);
        mftTracksTreeFiltered->Branch("MFTTrackClusIdx", &mftTrackClusIdxFiltered);

        mftClustersFileFiltered->cd();
        mftClustersTreeFiltered = new TTree("o2sim", "MFT clusters");
        mftClustersTreeFiltered->Branch("MFTClusterComp", &mftClustersFiltered);
        mftClustersTreeFiltered->Branch("MFTClustersROF", &mftClustersROFFiltered);
        mftClustersTreeFiltered->Branch("MFTClusterPatt", &mftClusterPatternsFiltered);

        itsClustersFileFiltered->cd();
        itsClustersTreeFiltered = new TTree("o2sim", "ITS clusters");
        itsClustersTreeFiltered->Branch("ITSClusterComp", &itsClustersFiltered);
        itsClustersTreeFiltered->Branch("ITSClustersROF", &itsClustersROFFiltered);
        itsClustersTreeFiltered->Branch("ITSClusterPatt", &itsPatternsFiltered);
      }

      if (mftTracksTreeFiltered) {
        mftTracksTreeFiltered->Fill();
      }
      if (mftClustersFileFiltered) {
        mftClustersTreeFiltered->Fill();
      }
      if (itsClustersFileFiltered) {
        itsClustersTreeFiltered->Fill();
      }

      nTracks += 1;
    }
  }

  if (mftTracksFileFiltered) {
    mftTracksFileFiltered->cd();
    mftTracksTreeFiltered->Write();
  }
  if (mftClustersFileFiltered) {
    mftClustersFileFiltered->cd();
    mftClustersTreeFiltered->Write();
  }
  if (itsClustersFileFiltered) {
    itsClustersFileFiltered->cd();
    itsClustersTreeFiltered->Write();
  }

  std::cout << nTracks << std::endl;
}
