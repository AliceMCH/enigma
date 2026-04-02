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

#include "MFTTracking/ROframe.h"
#include "MFTTracking/IOUtils.h"
//#include "MFTTracking/Tracker.h"
#include "MFTTracking/TrackCA.h"

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"


#endif

namespace fs = std::filesystem;

using namespace o2::mft;


struct AlignConfigHelper {
  int minPoints = 6;                ///< mininum number of clusters in a track used for alignment
  int chi2CutNStdDev = 3;           ///< Number of standard deviations for chi2 cut
  double residualCutInitial = 100.; ///< Cut on residual on first iteration
  double residualCut = 100.;        ///< Cut on residual for other iterations
  double allowedVarDeltaX = 0.5;    ///< allowed max delta in x-translation (cm)
  double allowedVarDeltaY = 0.5;    ///< allowed max delta in y-translation (cm)
  double allowedVarDeltaZ = 0.5;    ///< allowed max delta in z-translation (cm)
  double allowedVarDeltaRz = 0.01;  ///< allowed max delta in rotation around z-axis (rad)
  double chi2CutFactor = 256.;      ///< used to reject outliers i.e. bad tracks with sum(chi2) > Chi2DoFLim(fNStdDev, nDoF) * fChi2CutFactor
};

enum class AlignmentStatus {
  idealgeo = 0,       // ideal geometry, no alignment corrections
  prealign,           // aligned geometry with pre-aligned corrections
  pass1,              // aligned geometry with pass 1 MillePede corrections computed from and applied to pre-aligned geometry
  newpass1,           // aligned geometry with pass 1 MillePede corrections expressed w.r.t. and applied to ideal geometry
  pass2,              // aligned geometry with pass 2 MillePede corrections computed from and applied to ideal geometry
  dcacorrpass2,       // pass2 aligned geometry updated with x,y translation of each half MFT to correct the DCA shifts computed from and applied to ideal geometry
  dcacorrpass2rotated // above + global rotation Rx, Ry of each half MFT to correct non-zero trend of DCA_x,y vs z_vertex
};

std::pair<double, double> GetTrackSlopes(const o2::mft::TrackMFT& fwdtrack)
{
  // Parameter conversion
  double alpha1, alpha3, x2, x3;

  x2 = fwdtrack.getPhi();
  x3 = fwdtrack.getTanl();

  auto sinx2 = TMath::Sin(x2);
  auto cosx2 = TMath::Cos(x2);

  alpha1 = cosx2 / x3;
  alpha3 = sinx2 / x3;

  return {alpha1, alpha3};
}


//_________________________________
bool addFilesToChains(TChain* itsclusterChain,
                      TChain* mftclusterChain,
                      TChain* mfttrackChain,
                      int fileStart, int fileStop)
{
  bool success = false;

  if (!itsclusterChain || !mftclusterChain || !mfttrackChain) {
    // LOG(error) << "at least one TChain (mftclusterChain, mfttrackChain) is a null pointer";
    success = false;
    // LOG(error) << "addFilesToChains() - aborted !";
    return success;
  }

  success = false;
  std::cout << "[addFilesToChains]" << " start " << fileStart << " , stop " << fileStop << std::endl;
  // LOG(info) << "ChainBuildMethod::continuous"
  //           << " , start " << fileStart
  //           << " , stop " << fileStop;

  int countFiles = 0;
  std::string generalPath = "ordered";
  for (int ii = fileStart; ii <= fileStop; ii++) {
    std::stringstream ss;

    ss << generalPath << "/"
        << std::setw(5) << std::setfill('0') << ii;

    std::string filePath = ss.str();
    if (std::filesystem::exists(Form("%s/o2clus_its-filtered.root", filePath.c_str())) &&
        std::filesystem::exists(Form("%s/mftclusters-filtered.root", filePath.c_str())) &&
        std::filesystem::exists(Form("%s/mfttracks-filtered.root", filePath.c_str()))) {
      itsclusterChain->Add(Form("%s/o2clus_its-filtered.root", filePath.c_str()));
      mftclusterChain->Add(Form("%s/mftclusters-filtered.root", filePath.c_str()));
      mfttrackChain->Add(Form("%s/mfttracks-filtered.root", filePath.c_str()));
      std::cout << "Add " << Form("%s/o2clus_its-filtered.root", filePath.c_str()) << std::endl;
      std::cout << "Add " << Form("%s/mftclusters-filtered.root", filePath.c_str()) << std::endl;
      std::cout << "Add " << Form("%s/mfttracks-filtered.root", filePath.c_str()) << std::endl;
      // LOG(info) << "Add " << Form("%s/mftclusters.root", filePath.c_str());
      // LOG(info) << "Add " << Form("%s/mfttracks.root", filePath.c_str());
      countFiles++;
    } else {
      std::cout << "mftclusters-filtered.root or mfttracks-filtered.root not found at : " << filePath << std::endl;
      // LOG(error) << "mftclusters.root or mfttracks.root not found at : " << filePath;
    }
    std::cout << "number of files per chain = " << countFiles << std::endl;
  }
  if (countFiles) {
    success = true;
  }

  // LOG(info) << "addFilesToChains() - successful for run " << runN;
  return success;
}


bool checkTreeEntries(std::string mftTracksFileName,
                      std::string mftClustersFileName,
                      std::string itsClustersFileName)
{

  // MFT clusters
  TFile mftTracksFile(mftTracksFileName.data());
  TTree* mftTracksTree = (TTree*)mftTracksFile.Get("o2sim");

  // MFT clusters
  TFile mftClustersFile(mftClustersFileName.data());
  TTree* mftClustersTree = (TTree*)mftClustersFile.Get("o2sim");

  // ITS clusters
  TFile itsClustersFile(itsClustersFileName.data());
  TTree* itsClustersTree = (TTree*)itsClustersFile.Get("o2sim");

  std::cout << "Number of entries in ROOT trees:" << std::endl
      << "  MFT tracks:   " << mftTracksTree->GetEntries() << std::endl
      << "  MFT clusters: " << mftClustersTree->GetEntries() << std::endl
      << "  ITS clusters: " << itsClustersTree->GetEntries() << std::endl;

  if (mftTracksTree->GetEntries() != mftClustersTree->GetEntries() ||
      mftTracksTree->GetEntries() != itsClustersTree->GetEntries()) {
    std::cout << "Mismatch in number of entries in ROOT trees" << std::endl;
    return false;
  }

  return true;
}


void fillMftTrackEntries(std::string treeFileName,
                         std::vector<std::string>& chunks)
{
  TFile treeFile(treeFileName.data());
  TTree* tracksTree = (TTree*)treeFile.Get("o2sim");
  std::vector<o2::mft::TrackMFT>* mftTracks = nullptr;
  tracksTree->SetBranchAddress("MFTTrack", &mftTracks);

  for (int i = 0; i < tracksTree->GetEntries(); i++) {
    chunks.push_back(treeFileName);
  }
}


//_________________________________
bool addFilesToChainsUnordered(TChain* itsclusterChain,
                      TChain* mftclusterChain,
                      TChain* mfttrackChain,
                      std::vector<std::string>& chunks)
{
  bool success = false;

  if (!itsclusterChain || !mftclusterChain || !mfttrackChain) {
    // LOG(error) << "at least one TChain (mftclusterChain, mfttrackChain) is a null pointer";
    success = false;
    // LOG(error) << "addFilesToChains() - aborted !";
    return success;
  }

  success = false;
  // LOG(info) << "ChainBuildMethod::continuous"
  //           << " , start " << fileStart
  //           << " , stop " << fileStop;

  int countFiles = 0;
  std::string generalPath = "input";
  for (const auto & runEntry : fs::directory_iterator(generalPath)) {
    if (runEntry.path() != "input/569138") {
    //  continue;
    }
    for (const auto & chunkEntry : fs::directory_iterator(runEntry.path())) {
      for (const auto & tfEntry : fs::directory_iterator(chunkEntry.path())) {
        std::string filePath = tfEntry.path().c_str();
        //std::string mfttracksFileName = "mfttracks-filtered.root";
        std::string mfttracksFileName = "mfttracks-realign.root";
        std::cout << "Loading trees from " << filePath << std::endl;
        if (std::filesystem::exists(Form("%s/o2clus_its-filtered.root", filePath.c_str())) &&
            std::filesystem::exists(Form("%s/mftclusters-filtered.root", filePath.c_str())) &&
            std::filesystem::exists(Form("%s/%s", filePath.c_str(), mfttracksFileName.c_str()))) {

          fillMftTrackEntries(Form("%s/%s", filePath.c_str(), mfttracksFileName.c_str()), chunks);
          if (!checkTreeEntries(Form("%s/%s", filePath.c_str(), mfttracksFileName.c_str()),
                                Form("%s/mftclusters-filtered.root", filePath.c_str()),
                                Form("%s/o2clus_its-filtered.root", filePath.c_str()))) {
            continue;
          }

          itsclusterChain->Add(Form("%s/o2clus_its-filtered.root", filePath.c_str()));
          mftclusterChain->Add(Form("%s/mftclusters-filtered.root", filePath.c_str()));
          mfttrackChain->Add(Form("%s/%s", filePath.c_str(), mfttracksFileName.c_str()));
          std::cout << "Add " << Form("%s/o2clus_its-filtered.root", filePath.c_str()) << std::endl;
          std::cout << "Add " << Form("%s/mftclusters-filtered.root", filePath.c_str()) << std::endl;
          std::cout << "Add " << Form("%s/%s", filePath.c_str(), mfttracksFileName.c_str()) << std::endl;
          // LOG(info) << "Add " << Form("%s/mftclusters.root", filePath.c_str());
          // LOG(info) << "Add " << Form("%s/mfttracks.root", filePath.c_str());
          countFiles++;
        } else {
          if (!std::filesystem::exists(Form("%s/o2clus_its-filtered.root", filePath.c_str()))) {
            std::cout << "o2clus_its-filtered.root not found at : " << filePath << std::endl;
          }
          if (!std::filesystem::exists(Form("%s/mftclusters-filtered.root", filePath.c_str()))) {
            std::cout << "mftclusters-filtered.root not found at : " << filePath << std::endl;
          }
          if (!std::filesystem::exists(Form("%s/%s", filePath.c_str(), mfttracksFileName.c_str()))) {
            std::cout << mfttracksFileName << " not found at : " << filePath << std::endl;
          }
          // LOG(error) << "mftclusters.root or mfttracks.root not found at : " << filePath;
        }
        std::cout << "number of files per chain = " << countFiles << std::endl;
      }
    }
  }

  if (countFiles) {
    success = true;
  }

  // LOG(info) << "addFilesToChains() - successful for run " << runN;
  return success;
}

/*bool fitTracks(std::vector<o2::itsmft::CompClusterExt>& compClusters,
               std::vector<o2::itsmft::ROFRecord>& rofs,
               std::vector<unsigned char> patterns)
{
  //ROFFilter filter = [](const o2::itsmft::ROFRecord& r) { return true; };

  //auto tracker = std::make_unique<o2::mft::Tracker<TrackLTFL>>(false);
  auto pattIt = patterns.begin();

  auto iROF = 0;
  for (const auto& rof : rofs) {
    o2::mft::ROframe<o2::mft::TrackLTFL> roFrameData;
    //int nclUsed = ioutils::loadROFrameData(rof, roFrameData, compClusters, pattIt, mDict, nullptr, tracker.get(), filter);
    //std::cout << "ROframeId: " << iROF << ", clusters loaded : " << nclUsed << std::endl;
    iROF++;
  }
  return true;
}*/

template<class T>
void extrapolateTrack(const T& mftTrack, double slopex, double slopey, double z, double& x, double& y)
{
  int step = 4;
  double zOffset = 0;
  double xStretch = 1;

  if (step == 1) {
    zOffset = 0.18;
  }

  if (step == 2) {
    zOffset = 0.18;
    xStretch = 1.0018;
  }

  if (step == 3) {
    //zOffset = 0.222;
    zOffset = 0.25;
    xStretch = 1.0018;
  }

  if (step == 4) {
    //zOffset = 0.222;
    zOffset = 0.25;
    xStretch = 1.0;
  }

  double zMFT = mftTrack.getZ() + zOffset;
  x = (mftTrack.getX() * xStretch + 0.0) + slopex * xStretch * (z - zMFT);
  y = (mftTrack.getY() - 0.0) + slopey * (z - zMFT);
}

//_________________________________
// alienv setenv O2Physics/latest -c root -l
// .L ~/cernbox/alice/enigma/macros/runTracksToRecords.C++
// runTracksToRecords()
//_________________________________
void runItsMftAlign(int folderIdMin = 0, int folderIdMax = 50000,
                    const AlignmentStatus alignStatus = AlignmentStatus::dcacorrpass2,
                    //const AlignmentStatus alignStatus = AlignmentStatus::dcacorrpass2rotated,
                    const int minClusters = 5,
                    const bool doControl = true,
                    const int nEntriesAutoSave = 10000)
{

  ROOT::EnableImplicitMT(0);
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  gStyle->SetOptFit(111111);

  // geometry

  ///< load geometry from file
  ///< When applyMisalignedment == false --> read from unaligned file
  ///< When preferAlignedFile == true and applyMisalignment == true : Prefer reading from existing aligned file

  bool applyMisalignment = true;
  const bool preferAlignedFile = true;
  bool useMilleAlignment = true;

  switch (alignStatus) {
    case AlignmentStatus::idealgeo:
      applyMisalignment = false;
      useMilleAlignment = false;
      break;
    case AlignmentStatus::prealign:
      applyMisalignment = true;
      useMilleAlignment = false;
      break;
    case AlignmentStatus::pass1:
      applyMisalignment = true;
      useMilleAlignment = true;
      break;
    case AlignmentStatus::newpass1:
      applyMisalignment = true;
      useMilleAlignment = true;
      break;
    case AlignmentStatus::pass2:
      applyMisalignment = true;
      useMilleAlignment = true;
      break;
    case AlignmentStatus::dcacorrpass2:
      applyMisalignment = true;
      useMilleAlignment = true;
      break;
    default:
      applyMisalignment = true;
      useMilleAlignment = true;
      break;
  }

  // LOG(info) << "applyMisalignment = " << applyMisalignment;
  // LOG(info) << "preferAlignedFile = " << preferAlignedFile;
  o2::base::GeometryManager::loadGeometry("", applyMisalignment, preferAlignedFile);
  o2::its::GeometryTGeo* geom = o2::its::GeometryTGeo::Instance();
  geom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                                                 o2::math_utils::TransformType::T2GRot,
                                                 o2::math_utils::TransformType::L2G));

  // dictionary

  std::string itsDictFileName = "o2_its_dictionary.root";
  o2::itsmft::TopologyDictionary* itsDict = nullptr;
  try {
    itsDict = o2::itsmft::TopologyDictionary::loadFrom(itsDictFileName);
  } catch (std::exception e) {
    std::cout << "Error " << e.what() << std::endl;
    // LOG(error) << "runTracksToRecords() - aborted !";
    return;
  }
  itsDict->readFromFile(itsDictFileName);

  std::string mftDictFileName = "o2_mft_dictionary.root";
  o2::itsmft::TopologyDictionary* mftDict = nullptr;
  try {
    mftDict = o2::itsmft::TopologyDictionary::loadFrom(mftDictFileName);
  } catch (std::exception e) {
    std::cout << "Error " << e.what() << std::endl;
    // LOG(error) << "runTracksToRecords() - aborted !";
    return;
  }
  mftDict->readFromFile(mftDictFileName);

  // cluster and track chains

  TChain* itsclusterChain = new TChain("o2sim");
  TChain* mftclusterChain = new TChain("o2sim");
  TChain* mfttrackChain = new TChain("o2sim");

  // feed cluster and track chains
  std::vector<std::string> chunks;
  bool success = addFilesToChainsUnordered(itsclusterChain,
                                  mftclusterChain,
                                  mfttrackChain,
                                  chunks);

  if (!success) {
    std::cout << "runTracksToRecords() - aborted !" << std::endl;
    // LOG(error) << "runTracksToRecords() - aborted !";
    return;
  }

  std::cout << "Number of entries in ROOT trees:" << std::endl
      << "  ITS clusters: " << itsclusterChain->GetEntries() << std::endl
      << "  MFT clusters: " << mftclusterChain->GetEntries() << std::endl
      << "  MFT tracks:   " << mfttrackChain->GetEntries() << std::endl;

  // ITS clusters
  std::vector<o2::itsmft::CompClusterExt>* itsClusArr = nullptr;
  itsclusterChain->SetBranchAddress("ITSClusterComp", &itsClusArr);
  std::vector<unsigned char> itsPatternsVec, *itsPatternsVecP = &itsPatternsVec;
  itsclusterChain->SetBranchAddress("ITSClusterPatt", &itsPatternsVecP);

  //std::vector<unsigned char>* itsPatternsPtr = nullptr;
  //auto itsPattBranch = itsclusterChain->GetBranch("ITSClusterPatt");
  //if (itsPattBranch) {
  //  itsPattBranch->SetAddress(&itsPatternsPtr);
  //}

  // ITS ROFrecords
  std::vector<o2::itsmft::ROFRecord> itsRofRecVec, *itsRofRecVecP = &itsRofRecVec;
  itsclusterChain->SetBranchAddress("ITSClustersROF", &itsRofRecVecP);

  // MFT tracks
  std::vector<o2::mft::TrackMFT>* mftTracks = nullptr;
  mfttrackChain->SetBranchAddress("MFTTrack", &mftTracks);
  std::vector<o2::itsmft::ROFRecord>* mftTracksROF = nullptr;
  mfttrackChain->SetBranchAddress("MFTTracksROF", &mftTracksROF);
  std::vector<int>* mftTrackClusIdx = nullptr;
  mfttrackChain->SetBranchAddress("MFTTrackClusIdx", &mftTrackClusIdx);

  // MFT clusters
  std::vector<o2::itsmft::CompClusterExt>* mftClusterComp = nullptr;
  mftclusterChain->SetBranchAddress("MFTClusterComp", &mftClusterComp);
  std::vector<o2::itsmft::ROFRecord>* mftClustersROF = nullptr;
  mftclusterChain->SetBranchAddress("MFTClustersROF", &mftClustersROF);
  std::vector<unsigned char>* mftClusterPatt = nullptr;
  mftclusterChain->SetBranchAddress("MFTClusterPatt", &mftClusterPatt);



  // Histograms
  TFile outRootFile("itsMftAlign.root", "RECREATE");
  TH1F hTrackNClus("mftTrackNClus", "MFT track # of clusters", 10, 0.5, 10.5);
  TH2F hTrackXY("mftTrackXY", "MFT track Y vs. X", 150, -15, 15, 150, -15, 15);
  TH1F hPhiMFT("phiMFT", "#phi pf MFT track;#phi (rad)", 100, -TMath::Pi(), TMath::Pi());
  TH2F hDeltaPhi("deltaPhi", "#Delta#phi vs. #phi;#phi(degrees);#Delta#phi (rad)", 100, -180, 180, 100, -0.005, 0.005);
  TH2F hDeltaPhiVsDeltar("deltaPhiVsDeltar", "#Delta#phi vs. #Deltar;#Deltar (cm);#Delta#phi (rad)", 50, 0, 0.5, 100, -0.005, 0.005);
  TH2F hDeltaPhiVsClusx("deltaPhiVsClusx", "#Delta#phi vs. cluster x;cluster x (cm);#Delta#phi (rad)", 20, -40, 40, 100, -0.005, 0.005);
  TH2F hDeltaPhiVsClusy("deltaPhiVsClusy", "#Delta#phi vs. cluster y;cluster y (cm);#Delta#phi (rad)", 20, -40, 40, 100, -0.005, 0.005);
  TH2F hDeltaPhiVsMfty("deltaPhiVsMfty", "#Delta#phi vs. MFT y;MFT y (cm);#Delta#phi (rad)", 100, -10, 10, 100, -0.01, 0.01);
  TH2F hDeltaxVsClusterx("deltaxVsClusterx", "#Deltax vs. cluster x;cluster x (cm);#Deltax (cm)", 20, -40, 40, 100, -0.2, 0.2);
  TH2F hDeltaxVsClustery("deltaxVsClustery", "#Deltax vs. cluster y;cluster y (cm);#Deltax (cm)", 20, -40, 40, 100, -0.2, 0.2);
  TH2F hDeltayVsClusterx("deltayVsClusterx", "#Deltay vs. cluster x;cluster x (cm);#Deltay (cm)", 20, -40, 40, 100, -0.5, 0.5);
  TH2F hDeltayVsClustery("deltayVsClustery", "#Deltay vs. cluster y;cluster y (cm);#Deltay (cm)", 20, -40, 40, 100, -0.5, 0.5);
  //THnSparseF hDeltaX("deltaX", "#Deltax;#Deltax (cm);", 100, -TMath::Pi(), TMath::Pi(), 100, -0.01, 0.01);

  TH2F hDeltaxVsSlopex("deltaxVsSlopex", "#Deltax vs. slope x;slope x;#Deltax (cm)", 10, -0.5, 0.5, 100, -0.25, 0.25);
  TH2F hDeltayVsSlopey("deltayVsSlopex", "#Deltay vs. slope y;slope y;#Deltay (cm)", 10, -0.5, 0.5, 100, -0.25, 0.25);

  TH2F hDeltaxVsz("deltaxVsz", "#Deltax vs. z;z (cm);#Deltax (cm)", 10, -50, 50, 100, -0.25, 0.25);
  TH2F hDeltayVsz("deltayVsz", "#Deltay vs. z;z (cm);#Deltay (cm)", 10, -50, 50, 100, -0.25, 0.25);

  int ntf = mfttrackChain->GetEntries();
  //if (ntf > 100000) {
  //  ntf = 100000;
  //}
  for (int i = 0; i < ntf; i++) {
    itsclusterChain->GetEntry(i);
    mfttrackChain->GetEntry(i);

    //______________________________________________________
    for (int iMftRof = 0; iMftRof < mftTracksROF->size(); iMftRof++) {

      const auto& mftRofRec = (*mftTracksROF)[iMftRof];
      //std::cout << "mftRofRec.getNEntries(): " << mftRofRec.getNEntries() << std::endl;

      for (int itrk = 0; itrk < mftRofRec.getNEntries(); itrk++) {
        int trkEntry = mftRofRec.getFirstEntry() + itrk; // entry of icl-th cluster of this ROF in the vector of clusters
        const auto& mftTrack = (*mftTracks)[trkEntry];
        hTrackNClus.Fill(mftTrack.getNumberOfPoints());
        if (mftTrack.getNumberOfPoints() < minClusters) {
          continue;
        }
        auto slopes = GetTrackSlopes(mftTrack);

        std::cout << std::format("Entry #{} MFT track nClus={} ROF: orbit={} bc={}", i, mftTrack.getNumberOfPoints(), mftRofRec.getBCData().orbit,  mftRofRec.getBCData().bc) << std::endl;
        std::cout << std::format("Entry #{} -> {}", i, chunks[i]) << std::endl;

        double xTrackMFT = mftTrack.getX();
        double yTrackMFT = mftTrack.getY();
        double zTrackMFT = mftTrack.getZ();
        double rTrackMFT = std::sqrt(xTrackMFT * xTrackMFT + yTrackMFT * yTrackMFT);
        //double phiTrackMFT = mftTrack.getPhi(); //std::atan2(yTrackMFT, xTrackMFT);

        double xTrackMFTback;
        double yTrackMFTback;
        extrapolateTrack(mftTrack, slopes.first, slopes.second, -77.5, xTrackMFTback, yTrackMFTback);

        // MFT tracks full contained in upper half
        //if (yTrackMFT < 0 || slopes.second > 0) {
        // MFT tracks full contained in upper half
        //if (yTrackMFT < 0) { continue; }
        //if (yTrackMFT > 0 && yTrackMFTback < 0) { continue; }

        // MFT tracks full contained in lower half
        //if (yTrackMFT > 0) { continue; }
        //if (yTrackMFT < 0 && yTrackMFTback > 0) { continue; }

        hTrackXY.Fill(xTrackMFT, yTrackMFT);
        hPhiMFT.Fill(mftTrack.getPhi());

        int nROFRec = (int)itsRofRecVec.size();
        std::cout << "itsRofRecVec.size():   " << itsRofRecVec.size() << std::endl;
        std::cout << "itsPatternsVec.size(): " << itsPatternsVec.size() << std::endl;
        auto pattIt = itsPatternsVec.cbegin();
        for (int irof = 0; irof < nROFRec; irof++) {
          const auto& rofRec = itsRofRecVec[irof];
          if (rofRec.getBCData().orbit != mftRofRec.getBCData().orbit || rofRec.getBCData().bc != mftRofRec.getBCData().bc) {
            continue;
          }
          std::cout << std::format("  ITS cluster ROF: orbit={} bc={} size={}", rofRec.getBCData().orbit, rofRec.getBCData().bc, rofRec.getNEntries()) << std::endl;
          for (int icl = 0; icl < rofRec.getNEntries(); icl++) {
            int clEntry = rofRec.getFirstEntry() + icl; // entry of icl-th cluster of this ROF in the vector of clusters
            const auto& cluster = (*itsClusArr)[clEntry];

            auto pattID = cluster.getPatternID();
            o2::math_utils::Point3D<float> locC;
            if (pattID == o2::itsmft::CompCluster::InvalidPatternID || itsDict->isGroup(pattID)) {
              //std::cout << "pattID: " << pattID << std::endl;
              o2::itsmft::ClusterPattern patt(pattIt);
              locC = itsDict->getClusterCoordinates(cluster, patt, false);
            } else {
              locC = itsDict->getClusterCoordinates(cluster);
            }
            auto chipID = cluster.getSensorID();
            // Transformation to the local --> global
            auto gloC = geom->getMatrixL2G(chipID) * locC;

            double xStretch = 1.0 / 1.0018;
            double xClus = gloC.X();
            double yClus = gloC.Y();
            double zClus = gloC.Z();
            double rClus = std::sqrt(xClus * xClus + yClus * yClus);

            // outer ITS layer
            if (rClus < 38) {
              continue;
            }
            xClus *= xStretch;

            //if (yClus > 0) { continue; }

            double xTrack; // = mftTrack.getX() + slopes.first * (gloC.Z() - mftTrack.getZ());
            double yTrack; // = mftTrack.getY() + slopes.second * (gloC.Z() - mftTrack.getZ());
            extrapolateTrack(mftTrack, slopes.first, slopes.second, zClus, xTrack, yTrack);
            double zTrack = gloC.Z();
            double rTrack = std::sqrt(xTrack * xTrack + yTrack * yTrack);
            double phiTrack = std::atan2(yTrack, xTrack);

            double phiClus = std::atan2(yClus, xClus);
            double dphi = phiTrack - phiClus;

            double dx = xTrack - xClus;
            double dy = yTrack - yClus;
            double dr = std::sqrt(dx * dx + dy * dy);

            // cut on distance
            //if (dr > 0.5) {
            //  continue;
            //}


            hDeltaPhi.Fill(phiClus * 180 / TMath::Pi(), dphi);

            hDeltaPhiVsMfty.Fill(yTrackMFT, dphi);

            hDeltaPhiVsClusx.Fill(xClus, dphi);
            hDeltaPhiVsClusy.Fill(yClus, dphi);

            hDeltaPhiVsDeltar.Fill(dr, dphi);

            hDeltaxVsClusterx.Fill(xClus, dx);
            hDeltaxVsClustery.Fill(yClus, dx);
            hDeltayVsClusterx.Fill(xClus, dy);
            hDeltayVsClustery.Fill(yClus, dy);

            hDeltaxVsSlopex.Fill(slopes.first, dx);
            hDeltayVsSlopey.Fill(slopes.second, dy);

            if (gloC.X() > 0) {
              hDeltaxVsz.Fill(gloC.Z(), dx);
              hDeltayVsz.Fill(gloC.Z(), dy);
            }

            //std::cout << std::format("    Cluster (loc): x={:+0.3f} y={:+0.3f} z={:+0.3f} chipID={}", locC.X(), locC.Y(), locC.Z(), chipID) << std::endl;
            //std::cout << std::format("    Track at MFT:  x={:+0.3f} y={:+0.3f} z={:+0.3f} r={:+0.3f} phi={:+0.3f} sx={:+0.3f} sy={:+0.3f}",
            //    xTrackMFT, yTrackMFT, zTrackMFT, rTrackMFT, phiTrackMFT, slopes.first, slopes.second) << std::endl;
            //std::cout << std::format("    Cluster (glo): x={:+0.3f} y={:+0.3f} z={:+0.3f} r={:+0.3f} phi={:+0.3f}", gloC.X(), gloC.Y(), gloC.Z(), rClus, phiClus) << std::endl;
            //std::cout << std::format("    Track at clus: x={:+0.3f} y={:+0.3f} z={:+0.3f} r={:+0.3f} phi={:+0.3f} sx={:+0.3f} sy={:+0.3f}",
            //    xTrack, yTrack, zTrack, rTrack, phiTrack, slopes.first, slopes.second) << std::endl;
            //std::cout << std::format("    Delta:         x={:+0.3f} y={:+0.3f} r={:+0.3f} phi={:+0.3f}", dx, dy, dr, dphi) << std::endl;
            //std:cout << "---------------------------------" << std::endl;
          }
        }
      }
    }

    //std::cout << "Entry " << i << " done." << std::endl;
  }

  // the end

  std::chrono::steady_clock::time_point stop_time = std::chrono::steady_clock::now();

  std::cout << "----------------------------------- " << endl;
  std::cout << "Total Execution time: \t\t"
            << std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time).count()
            << " seconds" << endl;

  outRootFile.Write();

  TCanvas c("c", "c", 1200, 800);

  hTrackXY.Draw("colz");
  c.SaveAs("itsMftAlign.pdf(");

  hTrackNClus.Draw();
  c.SaveAs("itsMftAlign.pdf");

  hPhiMFT.Draw();
  c.SaveAs("itsMftAlign.pdf");

  hDeltaPhi.Draw("colz");
  c.SaveAs("itsMftAlign.pdf");

  hDeltaPhiVsClusx.Draw("colz");
  c.SaveAs("itsMftAlign.pdf");

  TH1* hDeltaProj = (TH1*)hDeltaPhiVsClusx.ProjectionY();
  hDeltaProj->Draw("hist");
  hDeltaProj->Fit("gaus");
  c.SaveAs("itsMftAlign.pdf");

  hDeltaxVsClusterx.Draw("colz");
  c.SaveAs("itsMftAlign.pdf");

  hDeltaxVsSlopex.Draw("colz");
  c.SaveAs("itsMftAlign.pdf");

  TH1* hDeltaxProj = (TH1*)hDeltaxVsSlopex.ProjectionY();
  hDeltaxProj->Draw("hist");
  hDeltaxProj->Fit("gaus");
  c.SaveAs("itsMftAlign.pdf");

  hDeltayVsSlopey.Draw("colz");
  c.SaveAs("itsMftAlign.pdf");

  TH1* hDeltayProj = (TH1*)hDeltayVsSlopey.ProjectionY();
  hDeltayProj->Draw("hist");
  hDeltayProj->Fit("gaus");
  c.SaveAs("itsMftAlign.pdf");

  c.Clear();
  c.SaveAs("itsMftAlign.pdf");


  hDeltaPhiVsClusy.Draw("colz");
  c.SaveAs("itsMftAlign.pdf");

  hDeltaPhiVsMfty.Draw("colz");
  c.SaveAs("itsMftAlign.pdf");

  hDeltaPhiVsDeltar.Draw("colz");
  c.SaveAs("itsMftAlign.pdf");

  hDeltaxVsClusterx.Draw("colz");
  c.SaveAs("itsMftAlign.pdf");

  hDeltaxVsClustery.Draw("colz");
  c.SaveAs("itsMftAlign.pdf");

  hDeltayVsClusterx.Draw("colz");
  c.SaveAs("itsMftAlign.pdf");

  hDeltayVsClustery.Draw("colz");
  c.SaveAs("itsMftAlign.pdf");

  hDeltaxVsSlopex.Draw("colz");
  c.SaveAs("itsMftAlign.pdf");

  hDeltayVsSlopey.Draw("colz");
  c.SaveAs("itsMftAlign.pdf");

  hDeltaxVsz.Draw("colz");
  c.SaveAs("itsMftAlign.pdf");

  hDeltayVsz.Draw("colz");
  c.SaveAs("itsMftAlign.pdf");

  c.Clear();
  c.SaveAs("itsMftAlign.pdf)");
}
