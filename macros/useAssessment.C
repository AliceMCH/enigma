#include <array>
#include <cstddef>
#include <iomanip>
#include <string>
#include <sstream>
#include <filesystem>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>

#include "Framework/Logger.h"

// B Off data reconstructed with various geometries
// LHC23e /alice/data/2023/LHC23e/535046/raw/0000
enum class AlignGeometry {
  kDCAcorrPass2 = 0,    // pass 2 + global x, y translations for each half MFT (L. Michelleti)
  kDCAcorrPass2Rotated, // above + global x, y rotation of each half MFT (A.Ferrero)
  kNumberOfAlignGeo
};

std::array<std::string, (unsigned int)AlignGeometry::kNumberOfAlignGeo> geoName = {
  "gDCAcorrPass2", "gDCAcorrPass2Rotated"};

std::array<std::unique_ptr<TH1F>, (unsigned int)AlignGeometry::kNumberOfAlignGeo> mTrackNumberOfClusters = {nullptr};
std::array<std::unique_ptr<TH1F>, (unsigned int)AlignGeometry::kNumberOfAlignGeo> mLTFTrackNumberOfClusters = {nullptr};

std::array<std::unique_ptr<TH1F>, (unsigned int)AlignGeometry::kNumberOfAlignGeo> mTrackChi2 = {nullptr};
std::array<std::unique_ptr<TH1F>, (unsigned int)AlignGeometry::kNumberOfAlignGeo> mTrackCharge = {nullptr};
std::array<std::unique_ptr<TH1F>, (unsigned int)AlignGeometry::kNumberOfAlignGeo> mTrackInvQPt = {nullptr};
std::array<std::unique_ptr<TH1F>, (unsigned int)AlignGeometry::kNumberOfAlignGeo> mTrackPhi = {nullptr};
std::array<std::unique_ptr<TH1F>, (unsigned int)AlignGeometry::kNumberOfAlignGeo> mTrackEta = {nullptr};

std::array<std::unique_ptr<TH2F>, (unsigned int)AlignGeometry::kNumberOfAlignGeo> mTrackXYNCls = {nullptr};
std::array<std::unique_ptr<TH2F>, (unsigned int)AlignGeometry::kNumberOfAlignGeo> mTrackEtaPhi = {nullptr};

constexpr unsigned int runN = 535046;
constexpr unsigned int fileStart[(unsigned int)AlignGeometry::kNumberOfAlignGeo] = {1, 1};
constexpr unsigned int fileStop[(unsigned int)AlignGeometry::kNumberOfAlignGeo] = {20, 44};
constexpr int minNClusters = 6;
std::string generalPath = "";

// ___________________________________
bool createHistos(const unsigned int ig = 0, const bool isZeroB = true)
{
  bool success = false;
  if (ig >= (unsigned int)AlignGeometry::kNumberOfAlignGeo) {
    success = false;
    LOG(error) << "createHistos() - Wrong geometry choice, ig = " << ig
               << ", but max is " << (unsigned int)AlignGeometry::kNumberOfAlignGeo;
    LOG(error) << "createHistos( ig = " << ig << ") - Aborted !";
    return success;
  }

  mTrackNumberOfClusters[ig] = std::make_unique<TH1F>(
    Form("mMFTTrackNumberOfClusters_%s", geoName[ig].c_str()),
    "Number Of Clusters Per Track; # clusters; # entries",
    10, 0.5, 10.5);
  mTrackNumberOfClusters[ig]->Sumw2();

  mLTFTrackNumberOfClusters[ig] = std::make_unique<TH1F>(
    Form("mMFTLTFTrackNumberOfClusters_%s", geoName[ig].c_str()),
    "Number Of Clusters Per LTF Track; # clusters; # entries",
    10, 0.5, 10.5);
  mTrackNumberOfClusters[ig]->Sumw2();

  mTrackChi2[ig] = std::make_unique<TH1F>(
    Form("mTrackChi2_%s", geoName[ig].c_str()),
    "Track #chi^{2}/NDF; #chi^{2}/NDF; # entries",
    210, -0.5, 20.5);
  mTrackChi2[ig]->Sumw2();

  if (!isZeroB) {
    mTrackCharge[ig] = std::make_unique<TH1F>(
      Form("mTrackCharge_%s", geoName[ig].c_str()),
      "Track Charge; q; # entries",
      3, -1.5, 1.5);
    mTrackCharge[ig]->Sumw2();

    mTrackInvQPt[ig] = std::make_unique<TH1F>(
      Form("mTrackInvQPt_%s", geoName[ig].c_str()),
      "Track q/p_{T}; q/p_{T} [1/GeV]; # entries",
      200, -10, 10);
    mTrackInvQPt[ig]->Sumw2();
  }

  mTrackPhi[ig] = std::make_unique<TH1F>(
    Form("mTrackPhi_%s", geoName[ig].c_str()),
    "Track #phi; #phi; # entries",
    100, -3.2, 3.2);
  mTrackPhi[ig]->Sumw2();

  mTrackEta[ig] = std::make_unique<TH1F>(
    Form("mTrackEta_%s", geoName[ig].c_str()),
    "Track #eta; #eta; # entries",
    50, -4, -2);
  mTrackEta[ig]->Sumw2();

  mTrackXYNCls[ig] = std::make_unique<TH2F>(
    Form("mMFTTrackXY_%d_MinClusters_%s", minNClusters, geoName[ig].c_str()),
    Form("Track Position (NCls >= %d); x; y", minNClusters),
    320, -16, 16, 320, -16, 16);
  mTrackXYNCls[ig]->Sumw2();
  mTrackXYNCls[ig]->SetOption("COLZ");

  mTrackEtaPhi[ig] = std::make_unique<TH2F>(
    Form("mMFTTrackEtaPhi_%d_MinClusters_%s", minNClusters, geoName[ig].c_str()),
    Form("Track #eta , #phi (NCls >= %d); #eta; #phi", minNClusters),
    50, -4, -2, 100, -3.2, 3.2);
  mTrackEtaPhi[ig]->Sumw2();
  mTrackEtaPhi[ig]->SetOption("COLZ");

  LOG(info) << "createHistos() - Done";
  success = true;
  return success;
}

// ___________________________________
bool getGeneralPath(const unsigned int ig = 0, const bool isZeroB = true)
{
  bool success = false;

  // common path for MFTAssessment.root files

  std::string basePath = "";
  std::string alignStatusPath = "";
  if (isZeroB) {
    basePath = "/Volumes/ugreen4tb/data/mft/2023/LHC23e";
  }
  switch (ig) { // file path depends on the used aligned geometry
    case (unsigned int)AlignGeometry::kDCAcorrPass2:
      alignStatusPath = "apass1_updated_muon_alignment";
      break;
    case (unsigned int)AlignGeometry::kDCAcorrPass2Rotated:
      alignStatusPath = "apass1_muon_alignment3";
      break;
    default:
      success = false;
      LOG(error) << "getGeneralPath() - Wrong geometry choice, ig = " << ig
                 << ", but max is " << (unsigned int)AlignGeometry::kNumberOfAlignGeo;
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
bool addHistos(const unsigned int ig = 0, const bool isZeroB = true)
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

    std::unique_ptr<TH1F> mMFTTrackNumberOfClusters(mftFile->Get<TH1F>("mMFTTrackNumberOfClusters"));
    mTrackNumberOfClusters[ig]->Add(mMFTTrackNumberOfClusters.get());

    std::unique_ptr<TH1F> mMFTTrackChi2(mftFile->Get<TH1F>("mMFTTrackChi2"));
    mTrackChi2[ig]->Add(mMFTTrackChi2.get());

    if (!isZeroB) {
      std::unique_ptr<TH1F> mMFTTrackCharge(mftFile->Get<TH1F>("mMFTTrackCharge"));
      mTrackCharge[ig]->Add(mMFTTrackCharge.get());

      std::unique_ptr<TH1F> mMFTTrackInvQPt(mftFile->Get<TH1F>("mMFTTrackInvQPt"));
      mTrackInvQPt[ig]->Add(mMFTTrackInvQPt.get());
    }

    std::unique_ptr<TH1F> mMFTTrackPhi(mftFile->Get<TH1F>("mMFTTrackPhi"));
    mTrackPhi[ig]->Add(mMFTTrackPhi.get());

    std::unique_ptr<TH1F> mMFTTrackEta(mftFile->Get<TH1F>("mMFTTrackEta"));
    mTrackEta[ig]->Add(mMFTTrackEta.get());

    std::unique_ptr<TH2F> mMFTTrackXYMinCls(mftFile->Get<TH2F>(
      Form("mMFTTrackXY_%d_MinClusters", minNClusters)));
    mTrackXYNCls[ig]->Add(mMFTTrackXYMinCls.get());

    std::unique_ptr<TH2F> mMFTTrackEtaPhiMinCls(mftFile->Get<TH2F>(
      Form("mMFTTrackEtaPhi_%d_MinClusters", minNClusters)));
    mTrackEtaPhi[ig]->Add(mMFTTrackEtaPhiMinCls.get());

    countFiles++;
  }
  if (countFiles) {
    success = true;
    LOG(info) << "addHistos() - successful for ig = " << ig;
  }
  return success;
}

// ___________________________________
void useAssessment(const bool isZeroB = true)
{
  bool success[(unsigned int)AlignGeometry::kNumberOfAlignGeo] = {false, false};
  for (unsigned int ig = 0; ig < (unsigned int)AlignGeometry::kNumberOfAlignGeo; ig++) {
    LOG(info) << "useAssessment() - ig = " << ig << " , " << geoName[ig];
    LOG(info) << "useAssessment() - isZeroB = " << isZeroB;
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
    std::unique_ptr<TFile> outFile(TFile::Open(filePath.c_str(), "RECREATE"));

    outFile->WriteObject(mTrackNumberOfClusters[ig].get(), "mMFTTrackNumberOfClusters");
    outFile->WriteObject(mTrackChi2[ig].get(), "mMFTTrackChi2");
    if (!isZeroB) {
      outFile->WriteObject(mTrackCharge[ig].get(), "mMFTTrackCharge");
      outFile->WriteObject(mTrackInvQPt[ig].get(), "mMFTTrackInvQPt");
    }
    outFile->WriteObject(mTrackPhi[ig].get(), "mMFTTrackPhi");
    outFile->WriteObject(mTrackEta[ig].get(), "mMFTTrackEta");
    outFile->WriteObject(mTrackXYNCls[ig].get(), Form("mMFTTrackXY_%d_MinClusters", minNClusters));
    outFile->WriteObject(mTrackEtaPhi[ig].get(), Form("mMFTTrackEtaPhi_%d_MinClusters", minNClusters));

    LOG(info) << "useAssessment() - Merged histograms written to output file " << filePath;
    outFile->Close();
  }
}
