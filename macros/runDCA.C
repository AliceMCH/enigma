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
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

#include "MFTAlignment/TracksToRecords.h"

#endif

namespace fs = std::filesystem;

double vertexShift = 0.25;


//_________________________________
// alienv setenv O2Physics/latest -c root -l
// .L ~/cernbox/alice/enigma/macros/runTracksToRecords.C++
// runTracksToRecords()
//_________________________________
void runDCA(std::string alienFolder)
{

  ROOT::EnableImplicitMT(0);
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  // field off data
  double Bz = 0;

  // Histograms

  TFile outputFile("DCA.root", "RECREATE");

  TH1F* hVertexZ = new TH1F("vertex_z", "Vertex z", 1000, -10.0, 10.0);
  TH2F* hDCAxy = new TH2F("muon-global-alignment/DCA/MFT/DCA_y_vs_x", "DCA y vs. x", 400, -0.5, 0.5, 400, -0.5, 0.5);

  TDirectory* subDir = outputFile.mkdir("muon-global-alignment");
  subDir = subDir->mkdir("DCA");
  TDirectory* outputDir = subDir->mkdir("MFT");

  float mftLadderWidth = 1.7;
  Int_t bins[5] = {400, 20, 30 * 4, 24 * 4, 6};
  Double_t xmin[5] = {-0.5, -10.0, -mftLadderWidth * 15.f / 2.f, -12.0, 5};
  Double_t xmax[5] = { 0.5,  10.0,  mftLadderWidth * 15.f / 2.f,  12.0, 11};
  THnSparseF* hDCAx = new THnSparseF("DCA_x", "DCA(x) vs. vz, tx, ty, nclus", 5, bins, xmin, xmax);
  THnSparseF* hDCAy = new THnSparseF("DCA_y", "DCA(y) vs. vz, tx, ty, nclus", 5, bins, xmin, xmax);


  std::string tempFolder = alienFolder;
  std::string ctfFolder = fs::path(tempFolder).filename();

  tempFolder = fs::path(tempFolder).parent_path();
  std::string tenMinutesFolder = fs::path(tempFolder).filename();

  tempFolder = fs::path(tempFolder).parent_path();
  tempFolder = fs::path(tempFolder).parent_path();
  std::string runFolder = fs::path(tempFolder).filename();
  std::string localFolder = std::string("input/") + runFolder + "/" + tenMinutesFolder + "/" + ctfFolder;

  fs::create_directories(localFolder);


  std::string vertexFileName = "o2_primary_vertex.root";
  bool vertexFileIsTemp = false;
  if (!std::filesystem::exists(Form("%s/%s", localFolder.c_str(), vertexFileName.c_str()))) {
    // copy from alien
    std::string copyCommand = std::format("alien_cp alien://{}/o2_primary_vertex.root file://./{}/o2_primary_vertex.root", alienFolder, localFolder);
    std::cout << copyCommand << std::endl;
    std::system(copyCommand.c_str());
    vertexFileIsTemp = true;
  }

  std::string mfttracksFileName = "mfttracks.root";
  bool mfttracksFileIsTemp = false;
  if (!std::filesystem::exists(Form("%s/%s", localFolder.c_str(), mfttracksFileName.c_str()))) {
    // copy from alien
    std::string copyCommand = std::format("alien_cp alien://{}/mfttracks.root file://./{}/mfttracks.root", alienFolder, localFolder);
    std::cout << copyCommand << std::endl;
    std::system(copyCommand.c_str());
    mfttracksFileIsTemp = true;
  }

  // vertex tree
  TFile vertexFile(Form("%s/%s", localFolder.c_str(), vertexFileName.c_str()));
  TTree* vertexTree = (TTree*)vertexFile.Get("o2sim");

  std::vector<o2::dataformats::PrimaryVertex>* vertices = nullptr;
  vertexTree->SetBranchAddress("PrimaryVertex", &vertices);

  // tracks tree
  TFile mfttracksFile(Form("%s/%s", localFolder.c_str(), mfttracksFileName.c_str()));
  TTree* mfttracksTree = (TTree*)mfttracksFile.Get("o2sim");

  std::vector<o2::mft::TrackMFT>* mftTracks = nullptr;
  mfttracksTree->SetBranchAddress("MFTTrack", &mftTracks);
  std::vector<o2::itsmft::ROFRecord>* mftTrackROFs = nullptr;
  mfttracksTree->SetBranchAddress("MFTTracksROF", &mftTrackROFs);

  if (vertexTree->GetEntries() != mfttracksTree->GetEntries()) {
    std::cout << "Vrtex and track trees have diffent number of entries" << std::endl;
    return;
  }


  for (int iEntry = 0; mfttracksTree->LoadTree(iEntry) >= 0; ++iEntry) {
     // Load the data for the given tree entry
    vertexTree->GetEntry(iEntry);
    mfttracksTree->GetEntry(iEntry);

    //______________________________________________________
    for (const auto& oneRof : *mftTrackROFs) { // track ROF loop
      const auto& rofStart = oneRof.getBCData();
      auto rofEnd = rofStart + 198;

      int vertexId = -1;
      int vertexCount = -1;
      int nVerticesInRof = 0;
      //______________________________________________________
      for (const auto& oneVertex : *vertices) { // vertex loop
        vertexCount += 1;
        if (oneVertex.getIRMin() < rofStart) {
          continue;
        }
        if (oneVertex.getIRMax() > rofEnd) {
          continue;
        }
        nVerticesInRof += 1;
        vertexId = vertexCount;

        if (false) {
          std::cout << std::format("[TOTO] correlated vertex found with IR=({},{} -> {},{})  MFT IR=({},{} -> {},{})",
              oneVertex.getIRMin().orbit, oneVertex.getIRMin().bc,
              oneVertex.getIRMax().orbit, oneVertex.getIRMax().bc,
              rofStart.orbit, rofStart.bc,
              rofEnd.orbit, rofEnd.bc) << std::endl;
        }
      }
      if (false && nVerticesInRof > 0) {
        std::cout << std::format("[TOTO] number of vertices in MFT ROF: {}", nVerticesInRof) << std::endl;
        std::cout << std::format("[TOTO] selected vertex ID: {}  ", vertexId) << std::endl;
      }

      // only select ROFs with a single ITS vertex
      if (nVerticesInRof != 1) {
        continue;
      }

      const auto& theVertex = (*vertices)[vertexId];

      // discard vertices in the +/- 5cm region
      auto pvz = theVertex.getZ();
      if (false) {
      std::cout << std::format("[TOTO] vertex z={:0.3f} IR=({},{} -> {},{})  MFT IR=({},{} -> {},{})",
          pvz,
          theVertex.getIRMin().orbit, theVertex.getIRMin().bc,
          theVertex.getIRMax().orbit, theVertex.getIRMax().bc,
          rofStart.orbit, rofStart.bc,
          rofEnd.orbit, rofEnd.bc) << std::endl;
      }

      hVertexZ->Fill(pvz);

      for (int iTrack = 0; iTrack != oneRof.getNEntries(); iTrack++) {
        const auto& mftTrack = (*mftTracks)[iTrack + oneRof.getFirstEntry()];
        auto mftTrackAtDCA = mftTrack;
        mftTrackAtDCA.propagateToZ(pvz - vertexShift, Bz);

        double dcax = mftTrackAtDCA.getX() - theVertex.getX();
        double dcay = mftTrackAtDCA.getY() - theVertex.getY();
        hDCAxy->Fill(dcax, dcay);
        hDCAx->Fill(dcax, pvz, mftTrack.getX(), mftTrack.getY(), mftTrack.getNumberOfPoints());
        hDCAy->Fill(dcay, pvz, mftTrack.getX(), mftTrack.getY(), mftTrack.getNumberOfPoints());
      }
    }
    //break;
  }

  // the end

  hVertexZ->Write();
  outputDir->cd();
  hDCAxy->Write();
  hDCAx->Write();
  hDCAy->Write();
  outputFile.Write();

  if (vertexFileIsTemp) {
    fs::remove(localFolder + "/" + vertexFileName);
  }
  if (mfttracksFileIsTemp) {
    fs::remove(localFolder + "/" + mfttracksFileName);
  }

  std::chrono::steady_clock::time_point stop_time = std::chrono::steady_clock::now();

  std::cout << "----------------------------------- " << endl;
  std::cout << "Total Execution time: \t\t"
            << std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time).count()
            << " seconds" << endl;
}
