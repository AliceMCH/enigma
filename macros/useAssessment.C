#include <array>
#include <cstddef>
#include <iomanip>
#include <random>
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

struct configHelper {
  unsigned int runN = 535046;
  unsigned int ig = 0;  // index of the current Alignment geometry
  int minNClusters = 6; // use tracks with at least 6 clusters
  bool isZeroB = true;  // is it B off data or not ?
  unsigned int fileStart = 1;
  unsigned int fileStop = 20;
  std::string generalPath = "";     // full path to subdirectories with MFTAssessment files
  std::string basePath = "";        // directory with all align passes for that runN
  std::string alignStatusPath = ""; // subdirectory name for the current Alignment geometry
};

configHelper mConfig;

std::unique_ptr<TH1F> mTrackNumberOfClusters = {nullptr};
std::unique_ptr<TH1F> mLTFTrackNumberOfClusters = {nullptr};

std::unique_ptr<TH1F> mTrackChi2 = {nullptr};
std::unique_ptr<TH1F> mTrackCharge = {nullptr};
std::unique_ptr<TH1F> mTrackInvQPt = {nullptr};
std::unique_ptr<TH1F> mTrackPhi = {nullptr};
std::unique_ptr<TH1F> mTrackEta = {nullptr};

std::unique_ptr<TH2F> mTrackXYNCls = {nullptr};
std::unique_ptr<TH2F> mTrackEtaPhi = {nullptr};

// ___________________________________
bool checkGeo(const unsigned int ig = 0)
{
  bool success = false;
  if (ig >= (unsigned int)AlignGeometry::kNumberOfAlignGeo) {
    success = false;
    LOG(error) << "checkGeo() - Wrong geometry choice, ig = " << ig
               << ", but max is " << (unsigned int)AlignGeometry::kNumberOfAlignGeo;
    return success;
  }
  success = true;
  return success;
};

// ___________________________________
bool getGeneralPath()
{
  bool success = false;
  if (!checkGeo(mConfig.ig)) {
    LOG(error) << "getGeneralPath( ig = " << mConfig.ig << ") - Aborted !";
    return success;
  }

  // common path for MFTAssessment.root files

  std::stringstream generalPathSs;
  generalPathSs << mConfig.basePath << "/" << mConfig.runN << "/" << mConfig.alignStatusPath;
  mConfig.generalPath = generalPathSs.str();
  std::filesystem::directory_entry entry{mConfig.generalPath};
  if (entry.exists()) {
    LOG(info) << "getGeneralPath() - Found requested path : " << mConfig.generalPath;
    success = true;
  } else {
    LOG(error) << "getGeneralPath() - Path does not exist : " << mConfig.generalPath;
    LOG(error) << "getGeneralPath( ig = " << mConfig.ig << ") - Aborted !";
    success = false;
  }
  return success;
}

// ___________________________________
bool fillConfig(const unsigned int ig = 0,
                const unsigned int runN = 535046,
                const int minNClusters = 6)
{
  bool success = false;
  mConfig.ig = ig;
  if (!checkGeo(mConfig.ig)) {
    LOG(error) << "fillConfig( ig = " << mConfig.ig << ") - Aborted !";
    return success;
  }
  mConfig.runN = runN;
  mConfig.minNClusters = minNClusters;
  if (mConfig.runN == 535046) {
    mConfig.isZeroB = true;
    mConfig.basePath = "/Volumes/ugreen4tb/data/mft/2023/LHC23e";
    switch (mConfig.ig) {
      case (unsigned int)AlignGeometry::kDCAcorrPass2:
        mConfig.fileStart = 1;
        mConfig.fileStop = 20;
        mConfig.alignStatusPath = "apass1_updated_muon_alignment";
        break;
      default: // i.e. AlignGeometry::kDCAcorrPass2Rotated:
        mConfig.fileStart = 1;
        mConfig.fileStop = 44;
        mConfig.alignStatusPath = "apass1_muon_alignment3";
    };
  }
  success = true;
  success *= getGeneralPath();
  return success;
}

// ___________________________________
bool createHistos()
{
  bool success = false;
  if (!checkGeo(mConfig.ig)) {
    LOG(error) << "createHistos( ig = " << mConfig.ig << ") - Aborted !";
    return success;
  }

  mTrackNumberOfClusters = std::make_unique<TH1F>(
    Form("mMFTTrackNumberOfClusters_%s", geoName[mConfig.ig].c_str()),
    "Number Of Clusters Per Track; # clusters; # entries",
    10, 0.5, 10.5);
  mTrackNumberOfClusters->Sumw2();

  mLTFTrackNumberOfClusters = std::make_unique<TH1F>(
    Form("mMFTLTFTrackNumberOfClusters_%s", geoName[mConfig.ig].c_str()),
    "Number Of Clusters Per LTF Track; # clusters; # entries",
    10, 0.5, 10.5);
  mLTFTrackNumberOfClusters->Sumw2();

  mTrackChi2 = std::make_unique<TH1F>(
    Form("mTrackChi2_%s", geoName[mConfig.ig].c_str()),
    "Track #chi^{2}/NDF; #chi^{2}/NDF; # entries",
    210, -0.5, 20.5);
  mTrackChi2->Sumw2();

  if (!mConfig.isZeroB) {
    mTrackCharge = std::make_unique<TH1F>(
      Form("mTrackCharge_%s", geoName[mConfig.ig].c_str()),
      "Track Charge; q; # entries",
      3, -1.5, 1.5);
    mTrackCharge->Sumw2();

    mTrackInvQPt = std::make_unique<TH1F>(
      Form("mTrackInvQPt_%s", geoName[mConfig.ig].c_str()),
      "Track q/p_{T}; q/p_{T} [1/GeV]; # entries",
      200, -10, 10);
    mTrackInvQPt->Sumw2();
  }

  mTrackPhi = std::make_unique<TH1F>(
    Form("mTrackPhi_%s", geoName[mConfig.ig].c_str()),
    "Track #phi; #phi; # entries",
    100, -3.2, 3.2);
  mTrackPhi->Sumw2();

  mTrackEta = std::make_unique<TH1F>(
    Form("mTrackEta_%s", geoName[mConfig.ig].c_str()),
    "Track #eta; #eta; # entries",
    50, -4, -2);
  mTrackEta->Sumw2();

  mTrackXYNCls = std::make_unique<TH2F>(
    Form("mMFTTrackXY_%d_MinClusters_%s", mConfig.minNClusters, geoName[mConfig.ig].c_str()),
    Form("Track Position (NCls >= %d); x; y", mConfig.minNClusters),
    320, -16, 16, 320, -16, 16);
  mTrackXYNCls->Sumw2();
  mTrackXYNCls->SetOption("COLZ");

  mTrackEtaPhi = std::make_unique<TH2F>(
    Form("mMFTTrackEtaPhi_%d_MinClusters_%s", mConfig.minNClusters, geoName[mConfig.ig].c_str()),
    Form("Track #eta , #phi (NCls >= %d); #eta; #phi", mConfig.minNClusters),
    50, -4, -2, 100, -3.2, 3.2);
  mTrackEtaPhi->Sumw2();
  mTrackEtaPhi->SetOption("COLZ");

  LOG(info) << "createHistos() - Done";
  success = true;
  return success;
}

// ___________________________________
bool addHistos()
{
  // common path for MFTAssessment.root files

  bool success = false;
  if (!getGeneralPath()) {
    return success;
  }

  // loop on MFTAssessment.root files

  unsigned int countFiles = 0;
  for (unsigned int ff = mConfig.fileStart; ff < mConfig.fileStop + 1; ff++) {
    std::stringstream ss;
    ss << mConfig.generalPath << "/" << ff << "/MFTAssessment.root";
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
    mTrackNumberOfClusters->Add(mMFTTrackNumberOfClusters.get());

    std::unique_ptr<TH1F> mMFTLTFTrackNumberOfClusters(mftFile->Get<TH1F>("mMFTLTFTrackNumberOfClusters"));
    mLTFTrackNumberOfClusters->Add(mMFTLTFTrackNumberOfClusters.get());

    std::unique_ptr<TH1F> mMFTTrackChi2(mftFile->Get<TH1F>("mMFTTrackChi2"));
    mTrackChi2->Add(mMFTTrackChi2.get());

    if (!mConfig.isZeroB) {
      std::unique_ptr<TH1F> mMFTTrackCharge(mftFile->Get<TH1F>("mMFTTrackCharge"));
      mTrackCharge->Add(mMFTTrackCharge.get());

      std::unique_ptr<TH1F> mMFTTrackInvQPt(mftFile->Get<TH1F>("mMFTTrackInvQPt"));
      mTrackInvQPt->Add(mMFTTrackInvQPt.get());
    }

    std::unique_ptr<TH1F> mMFTTrackPhi(mftFile->Get<TH1F>("mMFTTrackPhi"));
    mTrackPhi->Add(mMFTTrackPhi.get());

    std::unique_ptr<TH1F> mMFTTrackEta(mftFile->Get<TH1F>("mMFTTrackEta"));
    mTrackEta->Add(mMFTTrackEta.get());

    std::unique_ptr<TH2F> mMFTTrackXYMinCls(mftFile->Get<TH2F>(
      Form("mMFTTrackXY_%d_MinClusters", mConfig.minNClusters)));
    mTrackXYNCls->Add(mMFTTrackXYMinCls.get());

    std::unique_ptr<TH2F> mMFTTrackEtaPhiMinCls(mftFile->Get<TH2F>(
      Form("mMFTTrackEtaPhi_%d_MinClusters", mConfig.minNClusters)));
    mTrackEtaPhi->Add(mMFTTrackEtaPhiMinCls.get());

    countFiles++;
  }
  if (countFiles) {
    success = true;
    LOG(info) << "addHistos() - successful for ig = " << mConfig.ig;
  }
  return success;
}

// ___________________________________
void useAssessment(const unsigned int ig = 0,
                   const unsigned int runN = 535046,
                   const int minNClusters = 6)
{
  if (!fillConfig(ig, runN, minNClusters)) {
    return;
  }

  bool success = false;
  LOG(info) << "useAssessment() - ig = " << mConfig.ig << " , " << geoName[mConfig.ig];
  LOG(info) << "useAssessment() - isZeroB = " << mConfig.isZeroB;
  success = createHistos();
  success *= addHistos();
  if (!success) {
    LOG(error) << "useAssessment() - Not successful for ig = " << mConfig.ig;
    LOG(error) << "useAssessment() - Aborted !";
    return;
  }
  if (!getGeneralPath()) {
    LOG(error) << "useAssessment() - Aborted !";
    return;
  }
  std::stringstream ss;
  ss << mConfig.generalPath << "/merged-" << mConfig.ig << "-MFTAssessment.root";
  std::string filePath = ss.str();
  std::unique_ptr<TFile> outFile(TFile::Open(filePath.c_str(), "RECREATE"));

  outFile->WriteObject(mTrackNumberOfClusters.get(), "mMFTTrackNumberOfClusters");
  outFile->WriteObject(mLTFTrackNumberOfClusters.get(), "mMFTLTFTrackNumberOfClusters");
  outFile->WriteObject(mTrackChi2.get(), "mMFTTrackChi2");
  if (!mConfig.isZeroB) {
    outFile->WriteObject(mTrackCharge.get(), "mMFTTrackCharge");
    outFile->WriteObject(mTrackInvQPt.get(), "mMFTTrackInvQPt");
  }
  outFile->WriteObject(mTrackPhi.get(), "mMFTTrackPhi");
  outFile->WriteObject(mTrackEta.get(), "mMFTTrackEta");
  outFile->WriteObject(mTrackXYNCls.get(),
                       Form("mMFTTrackXY_%d_MinClusters", mConfig.minNClusters));
  outFile->WriteObject(mTrackEtaPhi.get(),
                       Form("mMFTTrackEtaPhi_%d_MinClusters", mConfig.minNClusters));

  LOG(info) << "useAssessment() - Merged histograms written to output file " << filePath;
  outFile->Close();
}
