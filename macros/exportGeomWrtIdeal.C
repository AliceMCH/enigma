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
// .L ~/cernbox/alice/enigma/macros/exploreGeom.C++
// exploreGeom()

//_________________________________
void exportGeomWrtIdeal(std::string alignParamFileName = "mft_alignment",
                        std::string OutAlignParamFileName = "mft_alignment_wrt_ideal")
{

  // load prealigned geometry (o2sim_geometry-aligned.root)

  const bool applyMisalignment = false;
  const bool preferAlignedFile = true;
  o2::base::GeometryManager::loadGeometry("", applyMisalignment, preferAlignedFile);
  o2::mft::GeometryTGeo* geom = o2::mft::GeometryTGeo::Instance();

  // load pass1 alignement parameters w.r.t prealigned geometry
  std::vector<o2::detectors::AlignParam> alignParameters;
  try {
    alignParameters = loadAlignParam(alignParamFileName);
  } catch (std::exception e) {
    LOG(fatal) << "Abort, " << e.what();
    return;
  }
  printAlignParam(alignParamFileName, alignParameters);

  // compute pass 1 alignment parameters w.r.t. ideal geometry

  std::vector<o2::detectors::AlignParam> outAlignParams;

  o2::mft::AlignSensorHelper chipHelper;
  double dRx = 0., dRy = 0., dRz = 0.; // delta rotations
  double dx = 0., dy = 0., dz = 0.;    // delta translations
  bool global = true;                  // delta in global ref. system

  o2::itsmft::ChipMappingMFT chipMappingMFT;
  const int numberOfSensors = o2::itsmft::ChipMappingMFT::NChips;

  for (int chipId = 0; chipId < numberOfSensors; chipId++) {

    chipHelper.setSensorOnlyInfo(chipId);
    dx = alignParameters[chipId].getX();
    dy = alignParameters[chipId].getY();
    dz = alignParameters[chipId].getZ();
    dRx = alignParameters[chipId].getPsi();
    dRy = alignParameters[chipId].getTheta();
    dRz = alignParameters[chipId].getPhi();

    TGeoHMatrix alignParamMatrix = [](auto dx, auto dy, auto dz, auto dRx, auto dRy, auto dRz) {
      TGeoHMatrix tmp;
      double tra[3] = {dx, dy, dz};
      tmp.SetTranslation(tra);
      double rot[9] = {};
      double sinpsi = std::sin(dRx);
      double cospsi = std::cos(dRx);
      double sinthe = std::sin(dRy);
      double costhe = std::cos(dRy);
      double sinphi = std::sin(dRz);
      double cosphi = std::cos(dRz);
      rot[0] = costhe * cosphi;
      rot[1] = -costhe * sinphi;
      rot[2] = sinthe;
      rot[3] = sinpsi * sinthe * cosphi + cospsi * sinphi;
      rot[4] = -sinpsi * sinthe * sinphi + cospsi * cosphi;
      rot[5] = -costhe * sinpsi;
      rot[6] = -cospsi * sinthe * cosphi + sinpsi * sinphi;
      rot[7] = cospsi * sinthe * sinphi + sinpsi * cosphi;
      rot[8] = costhe * cospsi;
      tmp.SetRotation(rot);
      return tmp;
    }(dx, dy, dz, dRx, dRy, dRz);

    TGeoHMatrix idealGlobalTransform, idealGlobalTransformInv;
    o2::base::GeometryManager::getOriginalMatrix(
      o2::detectors::DetID::MFT, chipId, idealGlobalTransform);
    idealGlobalTransformInv = idealGlobalTransform.Inverse();

    geom->fillMatrixCache(
      o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                               o2::math_utils::TransformType::L2G));
    TGeoHMatrix prealignedGlobalTransform = geom->getMatrixL2G(chipId);

    TGeoHMatrix new_matrix;

    // copy pass 1 w.r.t. prealigned

    new_matrix.SetTranslation(alignParamMatrix.GetTranslation());
    new_matrix.SetRotation(alignParamMatrix.GetRotationMatrix());

    // compute pass 1 w.r.t. ideal

    new_matrix.Multiply(&prealignedGlobalTransform);
    new_matrix.Multiply(&idealGlobalTransformInv);

    Double_t* tra = new_matrix.GetTranslation();
    dx = tra[0];
    dy = tra[1];
    dz = tra[2];

    Double_t* rot = new_matrix.GetRotationMatrix();
    dRx = std::atan2(-rot[5], rot[8]);
    dRy = std::asin(rot[2]);
    dRz = std::atan2(-rot[1], rot[0]);

    // store pass 1 w.r.t. ideal

    outAlignParams.emplace_back(
      chipHelper.geoSymbolicName(),
      chipHelper.sensorUid(),
      dx, dy, dz,
      dRx, dRy, dRz,
      global);
  }

  LOGF(info, "Storing *combined* MFT alignment params in local file %s.root",
       OutAlignParamFileName.c_str());
  TFile outAlignParamFile(Form("%s.root", OutAlignParamFileName.c_str()),
                          "recreate", "", 505);
  outAlignParamFile.WriteObjectAny(&outAlignParams,
                                   "std::vector<o2::detectors::AlignParam>",
                                   "ccdb_object");
  outAlignParamFile.Close();

  printAlignParam(OutAlignParamFileName, outAlignParams);
}
