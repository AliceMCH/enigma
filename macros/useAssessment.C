#include <array>
#include <iomanip>
#include <string>
#include <sstream>
#include <filesystem>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <TH1F.h>
#include <TH2F.h>

// B Off data reconstructed with various geometries
// LHC23e /alice/data/2023/LHC23e/535046/raw/0000
enum class AlignGeometry : unsigned int {
  kDCAcorrPass2 = 0,    // pass 2 + global x, y translations for each half MFT (L. Michelleti)
  kDCAcorrPass2Rotated, // above + global x, y rotation of each half MFT (A.Ferrero)
  kNumberOfAlignGeo
};

std::array<std::unique_ptr<TH1F>, AlignGeometry::kNumberOfAlignGeo> mTrackNumberOfClusters = {nullptr};
std::array<std::unique_ptr<TH1F>, AlignGeometry::kNumberOfAlignGeo> mTrackChi2 = {nullptr};
std::array<std::unique_ptr<TH1F>, AlignGeometry::kNumberOfAlignGeo> mTrackCharge = {nullptr};
std::array<std::unique_ptr<TH1F>, AlignGeometry::kNumberOfAlignGeo> mTrackInvQPt = {nullptr};
std::array<std::unique_ptr<TH2F>, AlignGeometry::kNumberOfAlignGeo> mTrackXYNCls = {nullptr};

constexpr unsigned int runN = 535046;
constexpr unsigned int fileStart[(unsigned int)AlignGeometry::kNumberOfAlignGeo] = {1, 1};
constexpr unsigned int fileStop[(unsigned int)AlignGeometry::kNumberOfAlignGeo] = {20, 44};
constexpr int minNClusters = 6;
std::string generalPath = "";

// ___________________________________
bool createHistos(const unsigned int ig = 0)
{
  bool success = false;
  if (ig >= (unsigned int)AlignGeometry::kNumberOfAlignGeo) {
    success = false;
    LOG(error) << "createHistos() - Wrong geometry choice, ig = " << ig << ", but max is " << AlignGeometry::kNumberOfAlignGeo;
    LOG(error) << "createHistos( ig = " << ig << ") - Aborted !";
    return success;
  }
  mTrackNumberOfClusters[ig] = std::make_unique<TH1F>(
    "mMFTTrackNumberOfClusters",
    "Number Of Clusters Per Track; # clusters; # entries",
    10, 0.5, 10.5);
  mTrackNumberOfClusters[ig]->Sumw2();
  mTrackChi2[ig] = std::make_unique<TH1F>(
    "mMFTTrackChi2",
    "Track #chi^{2}/NDF; #chi^{2}/NDF; # entries",
    210, -0.5, 20.5);
  mTrackChi2[ig]->Sumw2();
  mTrackCharge[ig] = std::make_unique<TH1F>(
    "mMFTTrackCharge",
    "Track Charge; q; # entries",
    3, -1.5, 1.5);
  mTrackCharge[ig]->Sumw2();
  mTrackInvQPt[ig] = std::make_unique<TH1F>(
    "mMFTTrackInvQPt",
    "Track q/p_{T}; q/p_{T} [1/GeV]; # entries",
    200, -10, 10);
  mTrackInvQPt[ig]->Sumw2();
  mTrackXYNCls[ig] = std::make_unique<TH2F>(
    Form("mMFTTrackXY_%d_MinClusters", minNClusters),
    Form("Track Position (NCls >= %d); x; y", minNClusters),
    320, -16, 16, 320, -16, 16);
  mTrackXYNCls[ig]->Sumw2();
  mTrackXYNCls[ig]->SetOption("COLZ");
  LOG(info) << "createHistos() - Done";
  success = true;
  return success;
}

// ___________________________________
bool getGeneralPath(const unsigned int ig = 0)
{
  bool success = false;

  // common path for MFTAssessment.root files

  std::string basePath = "/Volumes/ugreen4tb/data/mft/2023/LHC23e";
  std::string alignStatusPath = "";
  switch (ig) { // file path depends on the used aligned geometry
    case (unsigned int)AlignGeometry::kDCAcorrPass2:
      alignStatusPath = "apass1_updated_muon_alignment";
      break;
    case (unsigned int)AlignGeometry::kDCAcorrPass2Rotated:
      alignStatusPath = "apass1_muon_alignment3";
      break;
    default:
      success = false;
      LOG(error) << "getGeneralPath() - Wrong geometry choice, ig = " << ig << ", but max is " << AlignGeometry::kNumberOfAlignGeo;
      LOG(error) << "getGeneralPath( ig = " << ig << ") - Aborted !";
      return success;
  }
  std::stringstream generalPathSs;
  generalPathSs << basePath << "/" << runN << "/" << alignStatusPath;
  generalPath = generalPathSs.str();
  std::filesystem::directory_entry entry{generalPath};
  if (entry.exists()) {
    LOG(info) << "getGeneralPath() - Found requested path : " << generalPath;
    success = true;
  } else {
    LOG(error) << "getGeneralPath() - Path does not exist : " << generalPath;
    LOG(error) << "getGeneralPath( ig = " << ig << ") - Aborted !";
    success = false;
  }
  return success;
}

// ___________________________________
bool addHistos(const unsigned int ig = 0)
{

  // common path for MFTAssessment.root files

  bool success = false;
  if (!getGeneralPath(ig)) {
    return success;
  }

  // loop on MFTAssessment.root files

  unsigned int countFiles = 0;
  for (unsigned int ff = fileStart[ig]; ff < fileStop[ig] + 1; ff++) {
    std::stringstream ss;
    ss << generalPath << "/" << ff << "/MFTAssessment.root";
    std::string filePath = ss.str();
    if (!std::filesystem::exists(filePath.c_str())) {
      LOG(error) << "addHistos() - File not found : " << filePath;
      continue;
    }
    std::unique_ptr<TFile> mftFile(TFile::Open(filePath.c_str()));
    if (!mftFile || mftFile->IsZombie()) {
      LOG(error) << "addHistos() - Error opening file " << filePath;
      continue;
    }
    LOG(info) << "addHistos() - Add histos from " << filePath;

    TH1F* mMFTTrackNumberOfClusters = mftFile->Get("mMFTTrackNumberOfClusters");
    mTrackNumberOfClusters->Add(mMFTTrackNumberOfClustersOfClusters);

    TH1F* mMFTTrackChi2 = mftFile->Get("mMFTTrackChi2");
    mTrackChi2->Add(mMFTTrackChi2);

    TH1F* mMFTTrackCharge = mftFile->Get("mMFTTrackCharge");
    mTrackCharge->Add(mMFTTrackCharge);

    TH1F* mMFTTrackInvQPt = mftFile->Get("mMFTTrackInvQPt");
    mTrackInvQPt->Add(mMFTTrackInvQPt);

    TH2F* mMFTTrackXYMinCls = mftFile->Get(Form("mMFTTrackXY_%d_MinClusters", minNClusters));
    mTrackXYNCls->Add(mMFTTrackXYMinCls);

    countFiles++;
  }
  if (countFiles) {
    success = true;
    LOG(info) << "addHistos() - successful for ig = " << ig;
  }
  return success;
}

// ___________________________________
void useAssessment()
{
  bool success[(unsigned int)AlignGeometry::kNumberOfAlignGeo] = {false, false};
  for (unsigned int ig = 0; ig < (unsigned int)AlignGeometry::kNumberOfAlignGeo; ig++) {
    success[ig] = createHistos(ig);
    success[ig] *= addHistos(ig);
    if (!success[ig]) {
      LOG(error) << "useAssessment() - Not successful for ig = " << ig;
      continue;
    }
    if (!getGeneralPath(ig)) {
      continue;
    }
    std::stringstream ss;
    ss << generalPath << "/merged-" << ig << "-MFTAssessment.root";
    std::string filePath = ss.str();
    auto outFile = std::unique_ptr<TFile>(filePath.c_str(), "RECREATE");
    outFile->WriteObject(mMFTTrackNumberOfClusters[ig], "mMFTTrackNumberOfClusters");
    outFile->WriteObject(mMFTTrackChi2[ig], "mMFTTrackChi2");
    outFile->WriteObject(mMFTTrackCharge[ig], "mMFTTrackCharge");
    outFile->WriteObject(mMFTTrackInvQPt[ig], "mMFTTrackInvQPt");
    outFile->WriteObject(mTrackXYNCls, Form("mMFTTrackXY_%d_MinClusters", minNClusters));
    LOG(info) << "useAssessment() - Merged histograms written to output file " << filePath;
  }
}
