#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

#include <Rtypes.h>
//#include <ROOT/RDataFrame.hxx>
//#include <TFile.h>
//#include <TChain.h>
//#include <TTree.h>

#include "CommonUtils/NameConf.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

#include "MFTAlignment/AlignSensorHelper.h"
//#include "MFTAlignment/RecordsToAlignParams.h"

#endif

// alienv setenv O2Physics/latest -c root -l
// .L ~/cernbox/alice/enigma/macros/applyZShiftOnDisks.C++
// runAlign()

void applyZShiftOnDisks()
{
  ROOT::EnableImplicitMT(0);
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  // geometry

  const bool applyMisalignment = false;
  const bool preferAlignedFile = false;
  std::cout << std::format("GeometryManager::loadGeometry(\"\", {}, {})", applyMisalignment, preferAlignedFile) << std::endl;
  o2::base::GeometryManager::loadGeometry("", applyMisalignment, preferAlignedFile);
  o2::mft::GeometryTGeo* geom = o2::mft::GeometryTGeo::Instance();
  geom->fillMatrixCache(
    o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                             o2::math_utils::TransformType::L2G));


  std::string generalPath = ".";
  std::string alignParamFileName = "mft_alignment_wrt_ideal";

  // create vectors to store the dummy and new align params
  std::vector<o2::detectors::AlignParam> alignParamDummy;
  std::vector<o2::detectors::AlignParam> alignParamShifted;

  // fill the vectors of align param using loaded align params
  // and the disk z shifts
  // Currently only disk #2 is shifted upstream by 1.5mm

  std::array<double, 5> diskZShift{0, 0, 0.15, 0, 0};

  o2::itsmft::ChipMappingMFT chipMappingMFT;
  int NChips = o2::itsmft::ChipMappingMFT::NChips;
  o2::mft::AlignSensorHelper chipHelper;
  for (int iChip = 0; iChip < NChips; iChip++) {

    chipHelper.setSensorOnlyInfo(iChip);
    double dx = 0;
    double dy = 0;
    double dz = 0;
    double dzShifted = diskZShift[chipHelper.disk()];
    double dRx = 0;
    double dRy = 0;
    double dRz = 0;

    const bool global = true;
    alignParamDummy.emplace_back(
      chipHelper.geoSymbolicName(),
      chipHelper.sensorUid(),
      dx, dy, dz,
      dRx, dRy, dRz,
      global);
    alignParamShifted.emplace_back(
      chipHelper.geoSymbolicName(),
      chipHelper.sensorUid(),
      dx, dy, dzShifted,
      dRx, dRy, dRz,
      global);
  }

  LOGF(info, "Storing dummy MFT alignment params in local file %s/%s.root",
       generalPath.c_str(), alignParamFileName.c_str());
  TFile outAlignParamFile(Form("%s/%s.root", generalPath.c_str(), alignParamFileName.c_str()),
                          "recreate", "", 505);
  outAlignParamFile.WriteObjectAny(&alignParamShifted, "std::vector<o2::detectors::AlignParam>", o2::base::NameConf::CCDBOBJECT.data());
  outAlignParamFile.Close();

  // apply alignment

  bool isAlignApplied = o2::base::GeometryManager::applyAlignment(alignParamShifted);
  if (isAlignApplied) {
    std::cout << "Successfully applied alignment parameters" << std::endl;

    // generate aligned geometry file (o2sim_geometry-aligned.root)
    auto alignedgeomfile = o2::base::NameConf::getGeomFileName();
    gGeoManager->Export(Form("%s/new-%s", generalPath.c_str(), alignedgeomfile.c_str()));
    std::cout << "New geometry file generated : "
              << Form("%s/new-%s", generalPath.c_str(), alignedgeomfile.c_str())
              << std::endl;
  }

  // the end

  std::chrono::steady_clock::time_point stop_time = std::chrono::steady_clock::now();

  std::cout << "----------------------------------- " << endl;
  std::cout << "Total Execution time: \t\t"
            << std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time).count()
            << " seconds" << endl;
}
