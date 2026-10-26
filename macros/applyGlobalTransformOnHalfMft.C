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
// .L ~/cernbox/alice/enigma/macros/applyGlobalTransformOnHalfMft.C++
// applyGlobalTransformOnHalfMft()
//
// root -q -b '/eos/home-a/arakotoz/alice/enigma/macros/applyGlobalTransformOnHalfMft.C+(true)' > log.txt 2>&1

// compute mean DCA shift of a half MFT from 2 quadrants
//_________________________________
constexpr double mean_dca(double a, double b)
{
  return ((a + b) / 2.); // compute mean DCA shift of a half MFT from 2 quadrants
}

//_________________________________
void applyGlobalTransformOnHalfMft(
  double dcaXTop, double dcaXBottom, double dcaYTop, double dcaYBottom, double dcaZTop, double dcaZBottom,
  double RhTop, double RhBottom, double RvTop, double RvBottom,
  double zScalingTop, double zScalingBottom,
  //const bool isLxplus = false,
  std::string generalPath = ".", //"/Users/andry/cernbox/alice/enigma/common-input/lhc2022h/", //"/Users/andry/cernbox/alice/mft/pilotbeam/505713/out-mille/pass1",
  std::string alignParamFileName = "pass2_mft_alignment")
{
  /*std::string generalPath = "/Users/andry/cernbox/alice/enigma/common-input/lhc2022h/";
  if (isLxplus) {
    generalPath = "/eos/home-a/arakotoz/alice/enigma/common-input/lhc2022h/";
  }*/

  // load ideal geometry (o2sim_geometry-aligned.root)

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
  const std::array<double, NHalf> meanDcaZ{
    dcaZBottom, // Bottom MFT (phi < 0)
    dcaZTop     // Top MFT (phi > 0)
  };

  // values of the global angular corrections
  // from Andrea Ferrero (email from 2025, June 5)
  // global rotation around x axis:
  // - bottom half MFT -0.0006130 rad
  // - top half MFT -0.0023705 rad
  // global rotation around y axis:
  // - bottom half MFT 0.0010285 rad
  // - top half MFT -0.00061585 rad

  /*const std::array<double, NHalf> globalRx{
    -0.0006130,
    -0.0023705};
  const std::array<double, NHalf> globalRy{
    0.0010285,
    -0.00061585};*/

  const std::array<double, NHalf> globalRx{
    -RvBottom, -RvTop};
  const std::array<double, NHalf> globalRy{
    RhBottom, RhTop};

  // create vectors to store the new align params

  std::vector<o2::detectors::AlignParam> updatedAlignParam;

  // fill the vectors of align param for each half MFT
  // with global translations and global rotations

  double dx = 0., dy = 0., dz = 0., dRx = 0., dRy = 0., dRz = 0.;
  const bool global = true;

  for (int hIndex = 0; hIndex < NHalf; hIndex++) {
    TString geoSymbolicName = geom->composeSymNameHalf(hIndex);
    Int_t uid = -1; // default value of alignable ID (set for sensors only), see AlignParam.h
    dx = -meanDcaX[hIndex];
    dy = -meanDcaY[hIndex];
    dz = -meanDcaZ[hIndex];
    dRx = -globalRx[hIndex];
    dRy = -globalRy[hIndex];
    updatedAlignParam.emplace_back(
      geoSymbolicName,
      uid,
      dx, dy, dz,
      dRx, dRy, dRz,
      global);
    LOG(info) << "applyGlobalTransformOnHalfMft() - half " << hIndex
              << " name " << geoSymbolicName
              << " uid " << uid
              << " dx " << dx
              << " dy " << dy
              << " dRx " << dRx
              << " dRy " << dRy;
  }

  // fill the vectors of align param using loaded align params
  // for all chips

  const std::array<double, NHalf> zScaling{
    zScalingBottom, // Bottom MFT (phi < 0)
    zScalingTop     // Top MFT (phi > 0)
  };

  o2::itsmft::ChipMappingMFT chipMappingMFT;
  int NChips = o2::itsmft::ChipMappingMFT::NChips;

  o2::mft::AlignSensorHelper chipHelper;

  dx = 0., dy = 0., dz = 0., dRx = 0., dRy = 0., dRz = 0.;
  int iChip = 0;
  const double z0 = -45.3;

  for (auto const& aParam : alignParameters) {

    //chipHelper.setSensorOnlyInfo(iChip);
    chipHelper.setSensor(iChip);
    dx = aParam.getX();
    dy = aParam.getY();
    dz = aParam.getZ();
    dRx = aParam.getPsi();
    dRy = aParam.getTheta();
    dRz = aParam.getPhi();

    int hIndex = 0;

    if (chipHelper.half()) {
      // Top MFT h = 1
      hIndex = 1;
    } else {
      // Bottom MFT h = 0
      hIndex = 0;
    }
    double zAligned = chipHelper.translateZ() + dz;
    // rescale the z position according to the scaling factor
    double zAlignedScaled = (1. + zScaling[hIndex]) * (zAligned - z0);
    double zScalingAdjustment = zAlignedScaled - zAligned + z0;
    std::cout << std::format("zAligned={}  z0={}  zScaling={}  zAlignedScaled={}  zScalingAdjustment={}",
        zAligned, z0, (1. + zScaling[hIndex]), zAlignedScaled, zScalingAdjustment) << std::endl;
    updatedAlignParam.emplace_back(
      chipHelper.geoSymbolicName(),
      chipHelper.sensorUid(),
      dx, dy, dz + zScalingAdjustment,
      dRx, dRy, dRz,
      global);
    iChip++;
  }

  // save shifted alignment parameters to files

  LOGF(info, "Storing corrected align param vector in local file %s/%s_half_mft_rotated_and_shifted.root",
       generalPath.c_str(), alignParamFileName.c_str());
  TFile outAlignParamShiftedFile(
    Form("%s/%s_half_mft_rotated_and_shifted.root", generalPath.c_str(), alignParamFileName.c_str()),
    "recreate", "", 505);
  outAlignParamShiftedFile.WriteObjectAny(
    &updatedAlignParam, "std::vector<o2::detectors::AlignParam>",
    o2::base::NameConf::CCDBOBJECT.data());
  outAlignParamShiftedFile.Close();

  // print

  const bool wTranslation = true;
  const bool wRotation = true;
  const bool wDeg = false;

  const bool printScreen = false;
  printAlignParam(
    Form("%s/%s_half_mft_rotated_and_shifted", generalPath.c_str(), alignParamFileName.c_str()),
    updatedAlignParam,
    printScreen, wTranslation, wRotation, wDeg);
}
