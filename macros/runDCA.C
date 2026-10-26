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

double vertexShift = 0.225;
//double vertexShift = 0.25;


std::pair<double, double> getTrackSlopes(const o2::mft::TrackMFT& mftTrack)
{
  double phi = mftTrack.getPhi();
  double tanLambda = mftTrack.getTanl();
  return {std::cos(phi) / tanLambda, std::sin(phi) / tanLambda};
}

struct Point3D
{
  double x{0};
  double y{0};
  double z{0};
};

struct Line3D
{
  Point3D p1;
  Point3D p2;
  double sx{0};
  double sy{0};
};

void PCA(const Line3D& l1, const Line3D& l2, Point3D& pca, double& dca)
{
  // from https://math.stackexchange.com/a/4347294/122921
  // and https://palitri.com/vault/stuff/maths/Rays%20closest%20point.pdf

  // track 1 first point
  double x1 = l1.p1.x;
  double y1 = l1.p1.y;
  double z1 = l1.p1.z;
  // track 1 direction
  double x21 = l1.p2.x - l1.p1.x;
  double y21 = l1.p2.y - l1.p1.y;
  double z21 = l1.p2.z - l1.p1.z;
  //double x21 = l1.sx;
  //double y21 = l1.sy;
  //double z21 = 1;
  //std::cout << std::format("  p1: ({} {} {})  d1: ({} {} {})", x1, y1, z1, x21, y21, z21) << std::endl;

  //std::cout << std::format("  sx={} sy={}", l1.sx, l1.sy) << std::endl;

  // track 2 first point
  double x3 = l2.p1.x;
  double y3 = l2.p1.y;
  double z3 = l2.p1.z;
  // track 2 direction
  double x43 = l2.p2.x - l2.p1.x;
  double y43 = l2.p2.y - l2.p1.y;
  double z43 = l2.p2.z - l2.p1.z;
  //double x43 = l2.sx;
  //double y43 = l2.sy;
  //double z43 = 1;
  //std::cout << std::format("  p2: ({} {} {})  d2: ({} {} {})", x3, y3, z3, x43, y43, z43) << std::endl;

  //std::cout << std::format("  sx={} sy={}", l2.sx, l2.sy) << std::endl;

  //std::cout << std::format("closest points on Line(({:0.3f},{:0.3f},{:0.3f}),({:0.3f},{:0.3f},{:0.3f})) and Line(({:0.3f},{:0.3f},{:0.3f}),({:0.3f},{:0.3f},{:0.3f}))",
  //    l1.p1.x, l1.p1.y, l1.p1.z, l1.p2.x, l1.p2.y, l1.p2.z,
  //    l2.p1.x, l2.p1.y, l2.p1.z, l2.p2.x, l2.p2.y, l2.p2.z) << std::endl;
  //std::cout << std::format("close_points(Line(({},{},{}),({},{},{})), Line(({},{},{}),({},{},{})))",
  //    l1.p1.x, l1.p1.y, l1.p1.z, l1.p2.x, l1.p2.y, l1.p2.z,
  //    l2.p1.x, l2.p1.y, l2.p1.z, l2.p2.x, l2.p2.y, l2.p2.z) << std::endl;

  // distance between first points of track1 and track2
  double x31 = x3 - x1;
  double y31 = y3 - y1;
  double z31 = z3 - z1;
  //std::cout << std::format("  d3: ({} {} {})", x31, y31, z31) << std::endl;

  double R1 = x21 * x21 + y21 * y21 + z21 * z21;
  double R2 = x43 * x43 + y43 * y43 + z43 * z43;
  double d4321 = x21 * x43 + y21 * y43 + z21 * z43;
  double d3121 = x21 * x31 + y21 * y31 + z21 * z31;
  double d4331 = x31 * x43 + y31 * y43 + z31 * z43;
  double den = d4321 * d4321 - R1 * R2;
  //std::cout << std::format("  R1: {}", R1) << std::endl;
  //std::cout << std::format("  R2: {}", R2) << std::endl;
  //std::cout << std::format("  d4321: {}", d4321) << std::endl;
  //std::cout << std::format("  d3121: {}", d3121) << std::endl;
  //std::cout << std::format("  d4331: {}", d4331) << std::endl;
  //std::cout << std::format("  den: {}", den) << std::endl;

  // R1*s - d4321*t - d3121 = 0
  // d4321*s - R2*t - d4331 = 0
  double s = (d4321 * d4331 - R2 * d3121) / den;
  double t = (R1 * d4331 - d4321 * d3121) / den;

  // ps = l1.p1 + s*l1.direction
  double sx = x1 + s * x21;
  double sy = y1 + s * y21;
  double sz = z1 + s * z21;

  // pt = l2.p1 + t*l2.direction
  double tx = x3 + t * x43;
  double ty = y3 + t * y43;
  double tz = z3 + t * z43;

  pca.x = (sx + tx) / 2.0;
  pca.y = (sy + ty) / 2.0;
  pca.z = (sz + tz) / 2.0;

  dca = std::sqrt((sx - tx) * (sx - tx) + (sy - ty) * (sy - ty) + (sz - tz) * (sz - tz));

  //std::cout << std::format("  s={:0.3f}\n  t={:0.3f}\n  ps=({:0.3f} {:0.3f} {:0.3f})\n  pt=({:0.3f} {:0.3f} {:0.3f})\n  D={:0.3f}", s, t, sx, sy, sz, tx, ty, tz, dca) << std::endl;
}

void PCA(const o2::mft::TrackMFT& mftTrack1, const o2::mft::TrackMFT& mftTrack2, Point3D& pca, double& dca)
{

  // track 1 first point
  double x1 = mftTrack1.getX();
  double y1 = mftTrack1.getY();
  double z1 = mftTrack1.getZ();
  // track 1 direction
  auto slopes1 = getTrackSlopes(mftTrack1);
  double x21 = slopes1.first;
  double y21 = slopes1.second;
  double z21 = 1;

  Line3D l1;
  l1.p1.x = x1;
  l1.p1.y = y1;
  l1.p1.z = z1;
  l1.sx = x21;
  l1.sy = y21;
  l1.p2.x = x1 + x21 * 10;
  l1.p2.y = y1 + y21 * 10;
  l1.p2.z = z1 + 10;

  // track 2 first point
  double x3 = mftTrack2.getX();
  double y3 = mftTrack2.getY();
  double z3 = mftTrack2.getZ();
  // track 2 direction
  auto slopes2 = getTrackSlopes(mftTrack2);
  double x43 = slopes2.first;
  double y43 = slopes2.second;
  double z43 = 1;

  Line3D l2;
  l2.p1.x = x3;
  l2.p1.y = y3;
  l2.p1.z = z3;
  l2.sx = x43;
  l2.sy = y43;
  l2.p2.x = x3 + x43 * 10;
  l2.p2.y = y3 + y43 * 10;
  l2.p2.z = z3 + 10;

  PCA(l1, l2, pca, dca);
  return;

  // distance between first points of track1 and track2
  double x31 = x3 - x1;
  double y31 = y3 - y1;
  double z31 = z3 - z1;

  double R1 = x21 * x21 + y21 * y21 + z21 * z21;
  double R2 = x43 * x43 + y43 * y43 + z43 * z43;
  double d4321 = x21 * x43 + y21 * y43 + z21 * z43;
  double d3121 = x21 * x31 + y21 * y31 + z21 * z31;
  double d4331 = x31 * x43 + y31 * y43 + z31 * z43;
  double den = d4321 * d4321 - R1 * R2;

  // R1*s - d4321*t - d3121 = 0
  // d4321*s - R2*t - d4331 = 0
  double s = (d4321 * d4331 - R2 * d3121) / den;
  double t = (R1 * d4331 - d4321 * d3121) / den;

  // ps = l1.p1 + s*l1.direction
  double sx = x1 + s * x21;
  double sy = y1 + s * y21;
  double sz = z1 + s * z21;

  // pt = l2.p1 + t*l2.direction
  double tx = x3 + t * x43;
  double ty = y3 + t * y43;
  double tz = z3 + t * z43;

  pca.x = (sx + tx) / 2.0;
  pca.y = (sy + ty) / 2.0;
  pca.z = (sz + tz) / 2.0;

  dca = std::sqrt((sx - tx) * (sx - tx) + (sy - ty) * (sy - ty) + (sz - tz) * (sz - tz));
}


//_________________________________
// alienv setenv O2Physics/latest -c root -l
// .L ~/cernbox/alice/enigma/macros/runTracksToRecords.C++
// runTracksToRecords()
//_________________________________
void runDCA(std::string alienFolder, double _vertexShift = 0.225)
{

  ROOT::EnableImplicitMT(0);
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  vertexShift = _vertexShift;

  // field off data
  double Bz = 0;

  // Histograms

  TFile outputFile("DCA.root", "RECREATE");

  TH1F* hVertexZ = new TH1F("vertex_z", "Vertex z", 1000, -10.0, 10.0);
  TH1F* hVertexSigmaZ = new TH1F("VertexSigmaZ", "Vertex #sigma(z)", 100, 0, 0.05);
  TH2F* hTrackChi2 = new TH2F("TrackChi2", "MFT track #chi^{2} vs. vertex z", 20, -10, 10, 10000, 0, 100);

  TDirectory* subDir = outputFile.mkdir("muon-global-alignment");
  TDirectory* rootDir = subDir;
  subDir = subDir->mkdir("DCA");
  TDirectory* outputDir = subDir->mkdir("MFT");
  TDirectory* trackPairsDir = outputDir->mkdir("track-pairs");

  outputDir->cd();
  // DCA vs. track (X,Y) analysis
  float mftLadderWidth = 1.7;
  Int_t bins[5] = {400, 1, 30 * 4 / 2, 24 * 4 / 2, 6};
  Double_t xmin[5] = {-0.5, -10.0, -mftLadderWidth * 15.f / 2.f, -12.0, 5};
  Double_t xmax[5] = { 0.5,  10.0,  mftLadderWidth * 15.f / 2.f,  12.0, 11};
  THnSparseF* hDCAx = new THnSparseF("DCA_x", "DCA(x) vs. vz, tx, ty, nclus", 5, bins, xmin, xmax);
  THnSparseF* hDCAy = new THnSparseF("DCA_y", "DCA(y) vs. vz, tx, ty, nclus", 5, bins, xmin, xmax);
  TH2F* hDCAxy = new TH2F("DCA_y_vs_x", "DCA y vs. x", 400, -0.5, 0.5, 400, -0.5, 0.5);

  // Binning:
  // 0: DCA
  // 1: slope
  // 2: track x
  // 3: track y
  // 4: vertex z
  // 5: z shift
  // 6: bottom/top/other
  Int_t bins2[7] = {100, 50, 40, 40, 20, 21, 3};
  Double_t xmin2[7] = {-0.5, -0.25, -10, -10, -10, -5.25, 0};
  Double_t xmax2[7] = { 0.5,  0.25,  10,  10,  10,  5.25, 3};
  // DCA vs. vertex Z analysis
  auto* hDCAxVsVertexShift = new THnSparseF("DCA_x_vs_zshift_sparse", "DCA(x) vs. z shift;DCA_{x} (cm);x slope;track x (cm);track y (cm);vertex z (cm);z shift (mm);bottom/top/other", 7, bins2, xmin2, xmax2);
  auto* hDCAyVsVertexShift = new THnSparseF("DCA_y_vs_zshift_sparse", "DCA(y) vs. z shift;DCA_{y} (cm);y slope;track x (cm);track y (cm);vertex z (cm);z shift (mm);bottom/top/other", 7, bins2, xmin2, xmax2);
  /*TH3F* hDCAxVsVertexShift = new TH3F("DCA_x_vs_slopex_vs_zshift",
                                      "DCA(x) vs. x slope vs. z shift;z shift (mm);x slope;DCA_{x} (cm)",
                                      21, -5.25, 5.25,
                                      50, -0.25, 0.25,
                                      400, -0.5, 0.5);
  TH3F* hDCAxVsVertexShiftTop = new TH3F("DCA_x_vs_slopex_vs_zshift_top",
                                         "DCA(x) vs. x slope vs. z shift (top);z shift (mm);x slope;DCA_{x} (cm)",
                                         21, -5.25, 5.25,
                                         50, -0.25, 0.25,
                                         400, -0.5, 0.5);
  TH3F* hDCAxVsVertexShiftBottom = new TH3F("DCA_x_vs_slopex_vs_zshift_bottom",
                                            "DCA(x) vs. x slope vs. z shift (bottom);z shift (mm);x slope;DCA_{x} (cm)",
                                            21, -5.25, 5.25,
                                            50, -0.25, 0.25,
                                            400, -0.5, 0.5);

  TH3F* hDCAyVsVertexShift = new TH3F("DCA_y_vs_slopey_vs_zshift",
                                      "DCA(y) vs. y slope vs. z shift;z shift (mm);y slope;DCA_{y} (cm)",
                                      21, -5.25, 5.25,
                                      50, -0.25, 0.25,
                                      400, -0.5, 0.5);
  TH3F* hDCAyVsVertexShiftTop = new TH3F("DCA_y_vs_slopex_vs_zshift_top",
                                         "DCA(y) vs. y slope vs. z shift (top);z shift (mm);y slope;DCA_{y} (cm)",
                                         21, -5.25, 5.25,
                                         50, -0.25, 0.25,
                                         400, -0.5, 0.5);
  TH3F* hDCAyVsVertexShiftBottom = new TH3F("DCA_y_vs_slopex_vs_zshift_bottom",
                                            "DCA(y) vs. y slope vs. z shift (bottom);z shift (mm);y slope;DCA_{y} (cm)",
                                            21, -5.25, 5.25,
                                            50, -0.25, 0.25,
                                            400, -0.5, 0.5);*/

  trackPairsDir->cd();
  TH3F* hPairsDCAx = new TH3F("DCA_x", "DCA(x);vertex z (cm);angle (rad);DCA_{x} (cm)",
      20, -10.0, 10.0, 50, 0, 0.5, 400, -0.5, 0.5);
  TH3F* hPairsDCAy = new TH3F("DCA_y", "DCA(y);vertex z (cm);angle (rad);DCA_{y} (cm)",
      20, -10.0, 10.0, 50, 0, 0.5, 400, -0.5, 0.5);
  TH3F* hPairsDCAz = new TH3F("DCA_z", "DCA(z);vertex z (cm);angle (rad);DCA_{z} (cm)",
      20, -10.0, 10.0, 50, 0, 0.5, 400, -0.5, 0.5);
  TH3F* hTDCA = new TH3F("TDCA", "Track DCA vs. vz;vertex z (cm);angle (rad);DCA_{z} (cm)",
      20, -10.0, 10.0, 50, 0, 0.5, 400, 0, 0.5);



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
  std::cout << "Loading " << Form("%s/%s", localFolder.c_str(), vertexFileName.c_str()) << std::endl;
  TFile vertexFile(Form("%s/%s", localFolder.c_str(), vertexFileName.c_str()));
  TTree* vertexTree = (TTree*)vertexFile.Get("o2sim");
  std::cout << std::format("Number of entries: {}", vertexTree->GetEntries()) << std::endl;

  std::vector<o2::dataformats::PrimaryVertex>* vertices = nullptr;
  vertexTree->SetBranchAddress("PrimaryVertex", &vertices);

  // tracks tree
  std::cout << "Loading " << Form("%s/%s", localFolder.c_str(), mfttracksFileName.c_str()) << std::endl;
  TFile mfttracksFile(Form("%s/%s", localFolder.c_str(), mfttracksFileName.c_str()));
  TTree* mfttracksTree = (TTree*)mfttracksFile.Get("o2sim");
  std::cout << std::format("Number of entries: {}", mfttracksTree->GetEntries()) << std::endl;

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

      // discard vertices outside the +/- 10cm region
      auto pvz = theVertex.getZ();
      if (std::fabs(pvz) > 10) {
        continue;
      }

      pvz -= vertexShift;

      if (false) {
        std::cout << std::format("[TOTO] vertex z={:0.3f} IR=({},{} -> {},{})  MFT IR=({},{} -> {},{})",
            pvz,
            theVertex.getIRMin().orbit, theVertex.getIRMin().bc,
            theVertex.getIRMax().orbit, theVertex.getIRMax().bc,
            rofStart.orbit, rofStart.bc,
            rofEnd.orbit, rofEnd.bc) << std::endl;
      }

      hVertexZ->Fill(pvz);
      hVertexSigmaZ->Fill(theVertex.getSigmaZ());

      // require a vertex resolution in the Z coordinate better than 200 um
      if (theVertex.getSigmaZ() > 0.02) {
        continue;
      }

      for (int iTrack = 0; iTrack != oneRof.getNEntries(); iTrack++) {
        const auto& mftTrack = (*mftTracks)[iTrack + oneRof.getFirstEntry()];

        // select tracks on the upper-right part of MFT
        //if (mftTrack.getY() < 8 || mftTrack.getX() < 0) {
        //  continue;
        //}

        // skip CA tracks for the moment
        if (mftTrack.isCA()) {
          continue;
        }

        auto mftNclusters = mftTrack.getNumberOfPoints();

        auto mftTrackAtDCA = mftTrack;
        mftTrackAtDCA.propagateToZ(pvz, Bz);

        double tx = mftTrack.getX();
        double ty = mftTrack.getY();
        double dcax = mftTrackAtDCA.getX() - theVertex.getX();
        double dcay = mftTrackAtDCA.getY() - theVertex.getY();
        hDCAxy->Fill(dcax, dcay);
        hDCAx->Fill(dcax, pvz, tx, ty, mftNclusters);
        hDCAy->Fill(dcay, pvz, tx, ty, mftNclusters);
        //std::cout << std::format("MFT track at vertex: z={:0.3f}  dx={:0.3f}  dy={:0.3f}", pvz, dcax, dcay) << std::endl;

        //if (mftTrack.getChi2OverNDF() < 999.f && std::fabs(theVertex.getZ()) < 1.f && mftNclusters >= 6) {
        if (mftTrack.getChi2OverNDF() < 999.f && mftNclusters >= 6) {
          float zshift[21] = {// in millimeters
                              -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0,
                              0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
          for (int zi = 0; zi < 21; zi++) {
            auto mftTrackAtDCAshifted = mftTrack;
            mftTrackAtDCAshifted.propagateToZ(theVertex.getZ() - (zshift[zi] / 10.f), Bz);
            double dcaxShifted = mftTrackAtDCAshifted.getX() - theVertex.getX();
            double dcayShifted = mftTrackAtDCAshifted.getY() - theVertex.getY();

            std::pair<double, double> slopes = getTrackSlopes(mftTrackAtDCAshifted);
            auto slopex = slopes.first;
            auto slopey = slopes.second;
            //hDCAxVsVertexShift->Fill(dcaxShifted, slopex, mftTrack.getX(), mftTrack.getY(), theVertex.getZ(), zshift[zi]);
            //hDCAyVsVertexShift->Fill(dcayShifted, slopey, mftTrack.getX(), mftTrack.getY(), theVertex.getZ(), zshift[zi]);

            // top/bottom index for horizontal coordinate
            int tbx = 2;
            // top/bottom index for vertical coordinate
            int tby = 2;

            // remove a 1 cm cross around the horizontal and vertical axes
            if (std::fabs(tx) < 0.5 || std::fabs(ty) < 0.5) {
              hDCAxVsVertexShift->Fill(dcaxShifted, slopex, mftTrack.getX(), mftTrack.getY(), theVertex.getZ(), zshift[zi], tbx);
              hDCAyVsVertexShift->Fill(dcayShifted, slopey, mftTrack.getX(), mftTrack.getY(), theVertex.getZ(), zshift[zi], tby);
              continue;
            }

            if (std::fabs(ty) < 2.5) {
              if (ty > 0) {
                // DCA(x) top
                //hDCAxVsVertexShift->Fill(dcaxShifted, slopex, mftTrack.getX(), mftTrack.getY(), theVertex.getZ(), zshift[zi], 1);
                tbx = 1;
              } else {
                // DCA(x) bottom
                //hDCAxVsVertexShift->Fill(dcaxShifted, slopex, mftTrack.getX(), mftTrack.getY(), theVertex.getZ(), zshift[zi], 0);
                tbx = 0;
              }
            }
            hDCAxVsVertexShift->Fill(dcaxShifted, slopex, mftTrack.getX(), mftTrack.getY(), theVertex.getZ(), zshift[zi], tbx);

            if (std::fabs(tx) < 2.5) {
              if (ty > 0) {
                // DCA(y) top
                //hDCAyVsVertexShift->Fill(dcayShifted, slopey, mftTrack.getX(), mftTrack.getY(), theVertex.getZ(), zshift[zi], 1);
                tby = 1;
              } else {
                // DCA(y) bottom
                //hDCAyVsVertexShift->Fill(dcayShifted, slopey, mftTrack.getX(), mftTrack.getY(), theVertex.getZ(), zshift[zi], 0);
                tby = 0;
              }
            }
            hDCAyVsVertexShift->Fill(dcayShifted, slopey, mftTrack.getX(), mftTrack.getY(), theVertex.getZ(), zshift[zi], tby);
          }
        }

        if (mftNclusters >= 6) {
          hTrackChi2->Fill(theVertex.getZ(), mftTrack.getChi2OverNDF());

          auto slopes = getTrackSlopes(mftTrack);
          std::array<double, 3> direction{slopes.first, slopes.second, 1};
          double norm = std::sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
          direction[0] /= norm;
          direction[1] /= norm;
          direction[2] /= norm;

          for (int jTrack = iTrack + 1; jTrack != oneRof.getNEntries(); jTrack++) {
            const auto& mftTrack2 = (*mftTracks)[jTrack + oneRof.getFirstEntry()];

            // skip CA tracks for the moment
            if (mftTrack2.isCA()) {
              continue;
            }

            // skip tracks with less than 6 clusters
            if (mftTrack2.getNumberOfPoints() < 6) {
              continue;
            }

            auto slopes2 = getTrackSlopes(mftTrack2);
            std::array<double, 3> direction2{slopes2.first, slopes2.second, 1};
            double norm2 = std::sqrt(direction2[0] * direction2[0] + direction2[1] * direction2[1] + direction2[2] * direction2[2]);
            direction2[0] /= norm2;
            direction2[1] /= norm2;
            direction2[2] /= norm2;

            double angle = std::acos(direction[0] * direction2[0] + direction[1] * direction2[1] + direction[2] * direction2[2]);

            Point3D pca;
            double tdca;
            PCA(mftTrack, mftTrack2, pca, tdca);

            Line3D l1;
            l1.p1.x = 0;
            l1.p1.y = 0;
            l1.p1.z = 0;
            l1.p2.x = 1;
            l1.p2.y = 1;
            l1.p2.z = 1;
            l1.sx = 1;
            l1.sy = 1;

            Line3D l2;
            l2.p1.x = 0;
            l2.p1.y = 1;
            l2.p1.z = 0;
            //l2.p2.x = 1./3.;
            //l2.p2.y = 1. + 1./3.;
            //l2.p2.z = 3./3.;
            l2.p2.x = 1;
            l2.p2.y = 2;
            l2.p2.z = 3;
            l2.sx = 1. / 3.;
            l2.sy = 1. / 3.;

            //PCA(l1, l2, pca, pdca);

            double dcax = pca.x - theVertex.getX();
            double dcay = pca.y - theVertex.getY();
            double dcaz = pvz - pca.z;
            hPairsDCAx->Fill(pvz, angle, dcax);
            hPairsDCAy->Fill(pvz, angle, dcay);
            hPairsDCAz->Fill(pvz, angle, dcaz);
            hTDCA->Fill(pvz, angle, tdca);

            //std::cout << std::format("  PVz: {:0.3f}  DCAz: {:0.3f}", pvz, pca.z) << std::endl;
            //std::cout << std::format("  Pair DCA: {:0.3f}\n  DCAz:     {:0.3f}", tdca, dcaz) << std::endl;
            //std::cout << std::format("  angle={:0.3f}", angle) << std::endl;

            //return;
          }
        }
      }
    }
    //break;
  }

  // the end

  rootDir->cd();
  hVertexZ->Write();
  hVertexSigmaZ->Write();
  hTrackChi2->Write();
  outputDir->cd();
  hDCAxy->Write();
  hDCAx->Write();
  hDCAy->Write();
  hDCAxVsVertexShift->Write();
  hDCAyVsVertexShift->Write();
  trackPairsDir->cd();
  hPairsDCAx->Write();
  hPairsDCAy->Write();
  hPairsDCAz->Write();
  hTDCA->Write();
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
