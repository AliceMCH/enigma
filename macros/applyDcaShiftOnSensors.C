#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <iomanip>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <TFile.h>
#include <TObject.h>
#include <TString.h>
#include <Rtypes.h>

#include "CommonUtils/NameConf.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

#include "MFTAlignment/AlignSensorHelper.h"

#endif

#include "helperGeom.h"

// root -l
// .L ~/cernbox/alice/enigma/macros/applyDcaShiftOnSensors.C++
// applyDcaShiftOnSensors()
constexpr double mean_dca(double a, double b)
{
  return ((a + b) / 2.); // compute mean DCA shift of a half MFT from 2 quadrants
}

void applyDcaShiftOnSensors(
  double dcaXTop, double dcaXBottom, double dcaYTop, double dcaYBottom,
  std::string generalPath = ".", //"/Users/andry/cernbox/alice/enigma/common-input/lhc2022h/", //"/Users/andry/cernbox/alice/mft/pilotbeam/505713/out-mille/pass1",
  std::string alignParamFileName = "pass2_mft_alignment")
{
  // load ideal geometry (o2sim_geometry.root)

  const bool applyMisalignment = false;
  const bool preferAlignedFile = false;
  o2::base::GeometryManager::loadGeometry("", applyMisalignment, preferAlignedFile);
  o2::mft::GeometryTGeo* geom = o2::mft::GeometryTGeo::Instance();

  // cache matrix elements in the geometry

  geom->fillMatrixCache(
    o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                             o2::math_utils::TransformType::L2G));

  // load alignement parameters from alignParamFileName

  std::vector<o2::detectors::AlignParam> alignParameters;
  try {
    alignParameters = loadAlignParam(
      Form("%s/%s", generalPath.c_str(), alignParamFileName.c_str()));
  } catch (std::exception e) {
    LOG(fatal) << "Abort, " << e.what();
    return;
  }

  // values of the DCA shifts in x and y
  // from Luca Micheletti
  // pp collisions LHC23zs_pass1_mch (phi angle fixed)
  // global muons, muon cut chi2(MCH-MFT) < 45
  // Bottom half MFT (-pi < phi < 0) and Top half MFT (0 < phi < pi)
  // phi (rad)    : [-pi, -pi/2[ , [-pi/2, 0[ , [0, pi/2[ , [pi/2, pi[
  // <DCA x> (cm) : -0.0908      , -0.1004    , 0.0291    , 0.0369
  // <DCA y> (cm) : 0.0689       , 0.0716     , 0.0757    , 0.0705

  const int NQuadrant = 4; // number of MFT quadrants i.e. left, right of Top, Bottom MFT
  const int NHalf = 2; // number of half MFT (Top, Bottom)

  /*const std::array<double, NQuadrant> dcaX{
    -0.0908, -0.1004, // left, right of Bottom MFT (phi < 0)
    0.0291, 0.0369    // left, right of Top MFT (phi > 0)
  };
  const std::array<double, NQuadrant> dcaY{
    0.0689, 0.0716, // left, right of Bottom MFT (phi < 0)
    0.0757, 0.0705  // left, right of Top MFT (phi > 0)
  };

  const std::array<double, NHalf> meanDcaX{
    mean_dca(dcaX[0], dcaX[1]), // Bottom MFT (phi < 0)
    mean_dca(dcaX[2], dcaX[3])  // Top MFT (phi > 0)
  };
  const std::array<double, NHalf> meanDcaY{
    mean_dca(dcaY[0], dcaY[1]), // Bottom MFT (phi < 0)
    mean_dca(dcaY[2], dcaY[3])  // Top MFT (phi > 0)
  };*/

  const std::array<double, NHalf> meanDcaX{
    dcaXBottom, // Bottom MFT (phi < 0)
    dcaXTop     // Top MFT (phi > 0)
  };
  const std::array<double, NHalf> meanDcaY{
    dcaYBottom, // Bottom MFT (phi < 0)
    dcaYTop     // Top MFT (phi > 0)
  };

  // create vectors to store the new align params

  std::vector<o2::detectors::AlignParam> alignParamShifted;

  // fill the vectors of align param using loaded align params
  // and the dca shift

  o2::itsmft::ChipMappingMFT chipMappingMFT;
  int NChips = o2::itsmft::ChipMappingMFT::NChips;

  o2::mft::AlignSensorHelper chipHelper;

  double dx = 0., dy = 0., dz = 0., dRx = 0., dRy = 0., dRz = 0.;

  for (int iChip = 0; iChip < NChips; iChip++) {

    chipHelper.setSensorOnlyInfo(iChip);
    dx = alignParameters[iChip].getX();
    dy = alignParameters[iChip].getY();
    dz = alignParameters[iChip].getZ();
    dRx = alignParameters[iChip].getPsi();
    dRy = alignParameters[iChip].getTheta();
    dRz = alignParameters[iChip].getPhi();

    const bool global = true;
    int hIndex = 0;

    if (chipHelper.half()) {
      // Top MFT h = 1
      hIndex = 1;
    } else {
      // Bottom MFT h = 0
      hIndex = 0;
    }

    alignParamShifted.emplace_back(
      chipHelper.geoSymbolicName(),
      chipHelper.sensorUid(),
      dx - meanDcaX[hIndex], dy - meanDcaY[hIndex], dz,
      dRx, dRy, dRz,
      global);
  }

  // save shifted alignment parameters to files

  LOGF(info, "Storing corrected align param vector in local file %s/%s_sensor_shifted.root",
       generalPath.c_str(), alignParamFileName.c_str());
  TFile outAlignParamShiftedFile(
    Form("%s/%s_sensor_shifted.root", generalPath.c_str(), alignParamFileName.c_str()),
    "recreate", "", 505);
  outAlignParamShiftedFile.WriteObjectAny(
    &alignParamShifted, "std::vector<o2::detectors::AlignParam>",
    o2::base::NameConf::CCDBOBJECT.data());
  outAlignParamShiftedFile.Close();

  // print

  const bool wTranslation = true;
  const bool wRotation = true;
  const bool wDeg = false;

  const bool printScreen = true;

  printAlignParam(
    Form("%s/%s_sensor_shifted", generalPath.c_str(), alignParamFileName.c_str()),
    alignParamShifted,
    printScreen, wTranslation, wRotation, wDeg);
}
