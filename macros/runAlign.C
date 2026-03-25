#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

#include <Rtypes.h>
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>

#include "CommonUtils/NameConf.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

#include "MFTAlignment/AlignSensorHelper.h"
#include "MFTAlignment/RecordsToAlignParams.h"

#endif

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

// alienv setenv O2Physics/latest -c root -l
// .L ~/cernbox/alice/enigma/macros/runAlign.C++
// runAlign()

void runAlign(const double chi2CutFactor = 256, // 65536 // 256
              const bool preferAlignedFile = true,
              const bool useMilleAlignment = true, // true,
              const bool doControl = true,
              const int nEntriesAutoSave = 5000)
{
  ROOT::EnableImplicitMT(0);
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  // geometry

  const bool applyMisalignment = false;
  o2::base::GeometryManager::loadGeometry("", applyMisalignment, preferAlignedFile);
  o2::mft::GeometryTGeo* geom = o2::mft::GeometryTGeo::Instance();
  geom->fillMatrixCache(
    o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                             o2::math_utils::TransformType::L2G));

  // instantiate and configure the aligner

  AlignConfigHelper alignConfigParam;
  alignConfigParam.chi2CutFactor = chi2CutFactor;

  o2::mft::RecordsToAlignParams aligner;

  aligner.setChi2CutNStdDev(alignConfigParam.chi2CutNStdDev);
  aligner.setResidualCutInitial(alignConfigParam.residualCutInitial);
  aligner.setResidualCut(alignConfigParam.residualCut);
  aligner.setAllowedVariationDeltaX(alignConfigParam.allowedVarDeltaX);
  aligner.setAllowedVariationDeltaY(alignConfigParam.allowedVarDeltaY);
  aligner.setAllowedVariationDeltaZ(alignConfigParam.allowedVarDeltaZ);
  aligner.setAllowedVariationDeltaRz(alignConfigParam.allowedVarDeltaRz);
  aligner.setChi2CutFactor(alignConfigParam.chi2CutFactor);

  aligner.setWithControl(doControl);
  aligner.setNEntriesAutoSave(nEntriesAutoSave);

  // TODO: fix det. elements here

  // init Millipede

  aligner.init();

  // get chain of records

  std::string generalPath = ".";
  std::string alignParamFileName = "mft_alignment";
  std::string recordFileName = "mft_mille_records";
  TChain* recordChain = new TChain("o2sim");

  LOGF(info, "Load records from file %s/%s.root", generalPath.c_str(), recordFileName.c_str());
  recordChain->Add(Form("%s/%s.root", generalPath.c_str(), recordFileName.c_str()));

  // compute alignment parameters

  aligner.connectRecordReaderToChain(recordChain);
  aligner.globalFit();

  // save alignment parameters to file

  std::vector<o2::detectors::AlignParam> alignParams;
  aligner.getAlignParams(alignParams);
  LOGF(info, "Storing MFT alignment params in local file %s/%s.root",
       generalPath.c_str(), alignParamFileName.c_str());
  TFile outAlignParamFile(Form("%s/%s.root", generalPath.c_str(), alignParamFileName.c_str()),
                          "recreate", "", 505);
  outAlignParamFile.WriteObjectAny(&alignParams, "std::vector<o2::detectors::AlignParam>", o2::base::NameConf::CCDBOBJECT.data());
  outAlignParamFile.Close();

  // apply alignment

  bool isAlignApplied = o2::base::GeometryManager::applyAlignment(alignParams);
  if (isAlignApplied) {
    std::cout << "Successfully applied alignment parameters" << std::endl;

    // generate aligned geometry file (o2sim_geometry-aligned.root)
    auto alignedgeomfile = o2::base::NameConf::getAlignedGeomFileName();
    gGeoManager->Export(Form("%s/new-%s", generalPath.c_str(), alignedgeomfile.c_str()));
    std::cout << "New geometry file generated : "
              << Form("%s/new-%s", generalPath.c_str(), alignedgeomfile.c_str())
              << std::endl;
  }

  // save Pede outputs to ROOT and CSV files

  const int nDofPerSensor = aligner.getNDofPerSensor();
  std::vector<double> vPedeOutParams;
  std::vector<double> vPedeOutParamsErrors;
  std::vector<double> vPedeOutParamsPulls;
  aligner.getPedeOutParams(vPedeOutParams);
  aligner.getPedeOutParamsErrors(vPedeOutParamsErrors);
  aligner.getPedeOutParamsPulls(vPedeOutParamsPulls);
  TFile outPedefile(Form("%s/mft_pede_outputs.root", generalPath.c_str()), "recreate", "", 505);
  std::ofstream OutStream;
  OutStream.open(Form("%s/%s.csv", generalPath.c_str(), alignParamFileName.c_str()));
  OutStream << "half,disk,layer,zone,con,tr,chipid,dx,dy,dz,dRx,dRy,dRz,dxErr,dyErr,dzErr,dRxErr,dRyErr,dRzErr,dxPull,dyPull,dzPull,dRxPull,dRyPull,dRzPull"
            << endl;
  o2::mft::AlignSensorHelper chipHelper;
  double dx = 0., dy = 0., dz = 0.;
  double dxErr = 0., dyErr = 0., dzErr = 0.;
  double dxPull = 0., dyPull = 0., dzPull = 0.;
  double dRx = 0., dRy = 0., dRz = 0.;
  double dRxErr = 0., dRyErr = 0., dRzErr = 0.;
  double dRxPull = 0., dRyPull = 0., dRzPull = 0.;
  TTree* tree = new TTree("pede", "the Pede output tree");
  tree->Branch("dx", &dx, "dx/D");
  tree->Branch("dy", &dy, "dy/D");
  tree->Branch("dz", &dz, "dz/D");
  tree->Branch("dRz", &dRz, "dRz/D");
  tree->Branch("dxErr", &dxErr, "dxErr/D");
  tree->Branch("dyErr", &dyErr, "dyErr/D");
  tree->Branch("dzErr", &dzErr, "dzErr/D");
  tree->Branch("dRzErr", &dRzErr, "dRzErr/D");
  tree->Branch("dxPull", &dxPull, "dxPull/D");
  tree->Branch("dyPull", &dyPull, "dyPull/D");
  tree->Branch("dzPull", &dzPull, "dzPull/D");
  tree->Branch("dRzPull", &dRzPull, "dRzPull/D");
  o2::itsmft::ChipMappingMFT chipMappingMFT;
  int NChips = o2::itsmft::ChipMappingMFT::NChips;
  for (int iChip = 0; iChip < NChips; iChip++) {
    chipHelper.setSensorOnlyInfo(iChip);
    dx = vPedeOutParams[iChip * nDofPerSensor + 0];
    dy = vPedeOutParams[iChip * nDofPerSensor + 1];
    dz = vPedeOutParams[iChip * nDofPerSensor + 3];
    dRz = vPedeOutParams[iChip * nDofPerSensor + 2];
    dxErr = vPedeOutParamsErrors[iChip * nDofPerSensor + 0];
    dyErr = vPedeOutParamsErrors[iChip * nDofPerSensor + 1];
    dzErr = vPedeOutParamsErrors[iChip * nDofPerSensor + 3];
    dRzErr = vPedeOutParamsErrors[iChip * nDofPerSensor + 2];
    dxPull = vPedeOutParamsPulls[iChip * nDofPerSensor + 0];
    dyPull = vPedeOutParamsPulls[iChip * nDofPerSensor + 1];
    dzPull = vPedeOutParamsPulls[iChip * nDofPerSensor + 3];
    dRzPull = vPedeOutParamsPulls[iChip * nDofPerSensor + 2];
    tree->Fill();
    OutStream << chipHelper.half() << ","
              << chipHelper.disk() << ","
              << chipHelper.layer() << ","
              << chipHelper.zone() << ","
              << chipHelper.connector() << ","
              << chipHelper.transceiver() << ","
              << iChip << ","
              << dx << "," << dy << "," << dz << ","
              << dRx << "," << dRy << "," << dRz << ","
              << dxErr << "," << dyErr << "," << dzErr << ","
              << dRxErr << "," << dRyErr << "," << dRzErr << ","
              << dxPull << "," << dyPull << "," << dzPull << ","
              << dRxPull << "," << dRyPull << "," << dRzPull
              << endl;
  }
  tree->Write();
  outPedefile.Close();
  OutStream.close();

  // the end

  std::chrono::steady_clock::time_point stop_time = std::chrono::steady_clock::now();

  std::cout << "----------------------------------- " << endl;
  std::cout << "Total Execution time: \t\t"
            << std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time).count()
            << " seconds" << endl;
}
