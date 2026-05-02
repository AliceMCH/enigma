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

//_________________________________
std::vector<o2::detectors::AlignParam> loadAlignParam(std::string alignParamFileName)
{
  if (alignParamFileName.empty()) {
    LOG(fatal) << "No input align params file name provided !";
    throw std::exception();
  }
  LOG(info) << "Loading alignment parameters from " << alignParamFileName << ".root";
  TFile algFile(Form("%s.root", alignParamFileName.c_str()));
  if (algFile.IsZombie()) {
    LOG(fatal) << "Bad align param file " << alignParamFileName << ".root";
    throw std::exception();
  }
  std::vector<o2::detectors::AlignParam>* alignment;

  // use the standard object name
  algFile.GetObject(o2::base::NameConf::CCDBOBJECT.data(), alignment);

  // test alternative object names
  if (alignment == nullptr) {
    algFile.GetObject("alignment", alignment); // pass 1 and pass 2 millepede
  }
  if (alignment == nullptr) {
    algFile.GetObject("ccdb_object", alignment); // mft_survey_disk
  }

  // check that a valid align param vector was found
  algFile.Close();
  if (!alignment) {
    LOG(fatal) << "Empty vector of align params !";
    throw std::exception();
  }
  std::vector<o2::detectors::AlignParam> alignParameters = *alignment;
  return alignParameters;
}

//_________________________________
void printAlignParam(std::string alignParamFileName,
                     std::vector<o2::detectors::AlignParam> alignParameters,
                     const bool printScreen = false,
                     const bool wTranslation = true,
                     const bool wRotation = true,
                     const bool wDeg = false)
{
  o2::itsmft::ChipMappingMFT chipMappingMFT;
  int NChips = o2::itsmft::ChipMappingMFT::NChips;

  int alignParamsSize = static_cast<int>(alignParameters.size());
  LOG(info) << "printAlignParam() - vector of AlignParam, size = " << alignParamsSize;
  if (alignParamsSize == NChips) {
    LOG(info) << "printAlignParam() - same size as NChips = " << NChips;
  }

  std::ofstream outStream;
  outStream.open(Form("%s.csv", alignParamFileName.c_str()));
  outStream << "half,disk,layer,zone,con,tr,chipid,dx,dy,dz,dRx,dRy,dRz" << endl;

  o2::mft::AlignSensorHelper chipHelper;

  double dx = 0., dy = 0., dz = 0., dRx = 0., dRy = 0., dRz = 0.;
  int iChip = 0;
  int uid = -1;

  for (auto const& aParam : alignParameters) {

    dx = aParam.getX();
    dy = aParam.getY();
    dz = aParam.getZ();
    dRx = aParam.getPsi();
    dRy = aParam.getTheta();
    dRz = aParam.getPhi();

    uid = aParam.getAlignableID();
    if (uid > -1) { // --------------------------------------------------
      // it is a sensor
      if (iChip >= NChips) {
        LOGF(error, "printAlignParam() - chip index %d >= %d", iChip, NChips);
        continue;
      }

      chipHelper.setSensorOnlyInfo(iChip);
      outStream << chipHelper.half() << ","
                << chipHelper.disk() << ","
                << chipHelper.layer() << ","
                << chipHelper.zone() << ","
                << chipHelper.connector() << ","
                << chipHelper.transceiver() << ","
                << iChip << ","
                << dx << "," << dy << "," << dz << ","
                << dRx << "," << dRy << "," << dRz
                << endl;

      if (printScreen) {
        std::streamsize ss = std::cout.precision();
        bool wSymName = true;
        std::stringstream name = chipHelper.getSensorFullName(wSymName);
        std::cout << name.str().c_str()
                  << " (" << aParam.getSymName() << ") ";
        if (wTranslation) {
          std::cout << std::scientific << std::setprecision(2)
                    << " (cm) dx " << dx
                    << " dy " << dy
                    << " dz " << dz;
        }
        if (wRotation) {
          constexpr double rad2deg = 180.0 / 3.14159265358979323846;
          if (wDeg) {
            dRx *= rad2deg;
            dRy *= rad2deg;
            dRz *= rad2deg;
            std::cout << " (deg)";
          }
          std::cout << std::scientific << std::setprecision(2)
                    << " dRx " << dRx << " dRy " << dRy << " dRz " << dRz;
        }
        std::cout << std::setprecision(ss) << std::endl;
      }
      iChip++;

    } else { // -----------------------------------------------------------
      // not a sensor
      if (printScreen) {
        std::streamsize ss = std::cout.precision();
        std::cout << "not a sensor uid " << uid
                  << " (" << aParam.getSymName() << ") ";
        if (wTranslation) {
          std::cout << std::scientific << std::setprecision(2)
                    << " (cm) dx " << dx
                    << " dy " << dy
                    << " dz " << dz;
        }
        if (wRotation) {
          constexpr double rad2deg = 180.0 / 3.14159265358979323846;
          if (wDeg) {
            dRx *= rad2deg;
            dRy *= rad2deg;
            dRz *= rad2deg;
            std::cout << " (deg)";
          }
          std::cout << std::scientific << std::setprecision(2)
                    << " dRx " << dRx << " dRy " << dRy << " dRz " << dRz;
        }
        std::cout << std::setprecision(ss) << std::endl;
      }
    }
  }

  outStream.close();
}

//_________________________________
void printSensorGlobalTransform(std::ofstream& outStream,
                                const int iChip,
                                o2::mft::AlignSensorHelper chipHelper,
                                const bool wSymName = true,
                                const bool wTranslation = true,
                                const bool wRotation = true,
                                const bool wDeg = true)
{
  std::streamsize ss = std::cout.precision();
  std::stringstream name = chipHelper.getSensorFullName(wSymName);
  std::cout << name.str().c_str();
  outStream << chipHelper.half() << ","
            << chipHelper.disk() << ","
            << chipHelper.layer() << ","
            << chipHelper.zone() << ","
            << chipHelper.connector() << ","
            << chipHelper.transceiver() << ","
            << iChip;
  if (wTranslation) {
    outStream << "," << chipHelper.translateX()
              << "," << chipHelper.translateY()
              << "," << chipHelper.translateZ();
    std::cout << std::scientific << std::setprecision(2)
              << " (cm) dx " << chipHelper.translateX()
              << " dy " << chipHelper.translateY()
              << " dz " << chipHelper.translateZ();
  }
  if (wRotation) {
    constexpr double rad2deg = 180.0 / 3.14159265358979323846;
    double rotX = chipHelper.angleRx();
    double rotY = chipHelper.angleRy();
    double rotZ = chipHelper.angleRz();
    if (wDeg) {
      rotX *= rad2deg;
      rotY *= rad2deg;
      rotZ *= rad2deg;
      std::cout << " (deg)";
    }
    outStream << "," << rotX
              << "," << rotY
              << "," << rotZ;
    std::cout << std::scientific << std::setprecision(2)
              << " Rx " << rotX << " Ry " << rotY << " Rz " << rotZ;
  }
  outStream << endl;
  std::cout << std::setprecision(ss) << std::endl;
}
