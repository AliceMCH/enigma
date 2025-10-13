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

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

#include "MFTAlignment/TracksToRecords.h"

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

enum class ChainBuildMethod {
  continuous = 0,
  gridOk,
  gridFailed
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

enum class OneFourthChainPartId {
  first = 0,
  second,
  third,
  fourth,
  fifth,
  sixth,
  seventh,
  eighth,
  nineth,
  tenth,
  last,
  all
};

//_________________________________
bool addFilesToChains(TChain* mftclusterChain,
                      TChain* mfttrackChain,
                      const int runN = 505713,
                      const AlignmentStatus alignStatus = AlignmentStatus::pass2,
                      const bool useNewCTFs = true,
                      const OneFourthChainPartId partId = OneFourthChainPartId::first)
{
  bool success = false;

  if (!mftclusterChain || !mfttrackChain) {
    LOG(error) << "at least one TChain (mftclusterChain, mfttrackChain) is a null pointer";
    success = false;
    LOG(error) << "addFilesToChains() - aborted !";
    return success;
  }

  LOG(info) << "addFilesToChains() - start ... ";
  switch (alignStatus) {
    case AlignmentStatus::idealgeo:
      LOG(info) << "run " << runN << " , useNewCTFs = " << useNewCTFs << " , AlignmentStatus::idealgeo";
      break;
    case AlignmentStatus::prealign:
      LOG(info) << "run " << runN << " , useNewCTFs = " << useNewCTFs << " , AlignmentStatus::prealign";
      break;
    case AlignmentStatus::pass1:
      LOG(info) << "run " << runN << " , useNewCTFs = " << useNewCTFs << " , AlignmentStatus::pass1";
      break;
    case AlignmentStatus::newpass1:
      LOG(info) << "run " << runN << " , useNewCTFs = " << useNewCTFs << " , AlignmentStatus::newpass1";
      break;
    case AlignmentStatus::pass2:
      LOG(info) << "run " << runN << " , useNewCTFs = " << useNewCTFs << " , AlignmentStatus::pass2";
      break;
    case AlignmentStatus::dcacorrpass2:
      LOG(info) << "run " << runN << " , useNewCTFs = " << useNewCTFs << " , AlignmentStatus::dcacorrpass2";
      break;
    case AlignmentStatus::dcacorrpass2rotated:
      LOG(info) << "run " << runN << " , useNewCTFs = " << useNewCTFs << " , AlignmentStatus::dcacorrpass2rotated";
      break;
    default:
      LOG(info) << "run " << runN << " , useNewCTFs = " << useNewCTFs << " , AlignmentStatus " << int(alignStatus);
      break;
  }

  ChainBuildMethod mychoice = ChainBuildMethod::continuous;
  std::string generalPath = "";
  std::string basePath = "";
  std::string alignStatusPath = "";
  int fileStart = 1;
  int fileStop = 8;
  bool foundPath = false;

  if (runN == 505713) { // 2021 pilot beam
    foundPath = false;
    if (useNewCTFs) { // new CTFs
      if (alignStatus == AlignmentStatus::newpass1) {
        basePath = "/Volumes/ugreen4tb/data/mft/2021";
        alignStatusPath = "reco-with-mille/new-ctf/new-pass1";
        mychoice = ChainBuildMethod::continuous;
        fileStart = 1;
        fileStop = 819;
        foundPath = true;
      } else if (alignStatus == AlignmentStatus::prealign) {
        basePath = "/Volumes/ugreen4tb/data/mft/pilotbeam";
        alignStatusPath = "prealigned/new-ctf";
        mychoice = ChainBuildMethod::continuous;
        fileStart = 1;
        fileStop = 819;
        foundPath = true;
      } else if (alignStatus == AlignmentStatus::idealgeo) {
        basePath = "/Volumes/ugreen4tb/data/mft/pilotbeam";
        alignStatusPath = "ideal-geo/new-ctf";
        mychoice = ChainBuildMethod::gridOk;
        fileStart = 1;
        fileStop = 763;
        foundPath = true;
      } else {
        LOG(error) << "run " << runN
                   << " , useNewCTFs = " << useNewCTFs
                   << " , alignStatus = " << int(alignStatus)
                   << " , condition not found";
        foundPath = false;
      }
    } else { // old CTFs
      if (alignStatus == AlignmentStatus::pass1) {
        basePath = "/Volumes/ugreen4tb/data/mft/pilotbeam";
        alignStatusPath = "reco-with-mille/old-ctf/pass1";
        mychoice = ChainBuildMethod::continuous;
        fileStart = 1;
        fileStop = 44;
        foundPath = true;
      } else if (alignStatus == AlignmentStatus::prealign) {
        basePath = "/Volumes/ugreen4tb/data/mft/pilotbeam";
        alignStatusPath = "prealigned/old-ctf";
        mychoice = ChainBuildMethod::continuous;
        fileStart = 1;
        fileStop = 4315;
        foundPath = true;
      } else if (alignStatus == AlignmentStatus::pass2) {
        basePath = "/Volumes/ugreen4tb/data/mft/pilotbeam";
        basePath = "/Volumes/samt5/data/mft/pilotbeam";
        alignStatusPath = "reco-with-mille/old-ctf/pass2";
        mychoice = ChainBuildMethod::continuous;
        fileStart = 1;
        fileStop = 173;
        foundPath = true;
      } else {
        LOG(error) << "run " << runN
                   << " , useNewCTFs = " << useNewCTFs
                   << " , alignStatus = " << int(alignStatus)
                   << " , condition not found";
        foundPath = false;
      }
    }
  }

  if (runN == 505720) { // 2021 pilot beam
    LOG(info) << "run " << runN;
    foundPath = false;
    if (useNewCTFs && alignStatus == AlignmentStatus::idealgeo) {
      basePath = "/Volumes/ugreen4tb/data/mft/pilotbeam";
      alignStatusPath = "ideal-geo/new-ctf";
      mychoice = ChainBuildMethod::gridFailed;
      fileStart = 1;
      fileStop = 1619;
      foundPath = true;
    } else {
      LOG(error) << "run " << runN
                 << " , useNewCTFs = " << useNewCTFs
                 << " , alignStatus = " << int(alignStatus)
                 << " , condition not found";
      foundPath = false;
    }
  }

  if (runN == 520495) { // LHC22h
    LOG(info) << "run " << runN;
    foundPath = false;
    if (alignStatus == AlignmentStatus::pass1) {
      basePath = "/Volumes/ugreen4tb/data/mft/2022/LHC22h";
      alignStatusPath = "reco-with-mille/pass1";
      mychoice = ChainBuildMethod::continuous;
      fileStart = 1;
      fileStop = 8;
      foundPath = true;
    } else if (alignStatus == AlignmentStatus::pass2) {
      basePath = "/Volumes/ugreen4tb/data/mft/2022/LHC22h";
      alignStatusPath = "reco-with-mille/pass2";
      mychoice = ChainBuildMethod::continuous;
      fileStart = 1;
      fileStop = 8;
      foundPath = true;
    } else if (alignStatus == AlignmentStatus::idealgeo) {
      basePath = "/Volumes/ugreen4tb/data/mft/2022/LHC22h";
      alignStatusPath = "ideal-geo";
      mychoice = ChainBuildMethod::continuous;
      fileStart = 1;
      fileStop = 8;
      foundPath = true;
    } else {
      LOG(error) << "run " << runN
                 << " , useNewCTFs = " << useNewCTFs
                 << " , alignStatus = " << int(alignStatus)
                 << " , condition not found";
      foundPath = false;
    }
  }

  if (runN == 535046) { // LHC23e
    LOG(info) << "run " << runN;
    foundPath = false;
    if (alignStatus == AlignmentStatus::pass2) {
      basePath = "/Volumes/ugreen4tb/data/mft/2023/LHC23e";
      alignStatusPath = "apass1_muon_alignment";
      mychoice = ChainBuildMethod::continuous;
      fileStart = 188;
      fileStop = 198;
      // fileStop = 242;
      foundPath = true;
    } else if (alignStatus == AlignmentStatus::dcacorrpass2) {
      basePath = "/Volumes/ugreen4tb/data/mft/2023/LHC23e";
      alignStatusPath = "apass1_updated_muon_alignment";
      mychoice = ChainBuildMethod::continuous;
      switch (partId) {
        case OneFourthChainPartId::first:
          fileStart = 1;
          fileStop = 5;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::first";
          break;
        case OneFourthChainPartId::second:
          fileStart = 6;
          fileStop = 10;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::second";
          break;
        case OneFourthChainPartId::third:
          fileStart = 11;
          fileStop = 16;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::third";
          break;
        case OneFourthChainPartId::last:
          fileStart = 17;
          fileStop = 22;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::last";
          break;
        default:
          fileStart = 1;
          fileStop = 22;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::all";
          break;
      }
      foundPath = true;
    } else if (alignStatus == AlignmentStatus::dcacorrpass2rotated) {
      basePath = "/Volumes/ugreen4tb/data/mft/2023/LHC23e";
      alignStatusPath = "apass1_muon_alignment3";
      mychoice = ChainBuildMethod::continuous;
      switch (partId) {
        case OneFourthChainPartId::first:
          fileStart = 1;
          fileStop = 4;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::first";
          break;
        case OneFourthChainPartId::second:
          fileStart = 5;
          fileStop = 8;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::second";
          break;
        case OneFourthChainPartId::third:
          fileStart = 9;
          fileStop = 12;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::third";
          break;
        case OneFourthChainPartId::fourth:
          fileStart = 13;
          fileStop = 16;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::fourth";
          break;
        case OneFourthChainPartId::fifth:
          fileStart = 17;
          fileStop = 20;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::fifth";
          break;
        case OneFourthChainPartId::sixth:
          fileStart = 21;
          fileStop = 24;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::sixth";
          break;
        case OneFourthChainPartId::seventh:
          fileStart = 25;
          fileStop = 28;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::seventh";
          break;
        case OneFourthChainPartId::eighth:
          fileStart = 29;
          fileStop = 32;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::eighth";
          break;
        case OneFourthChainPartId::nineth:
          fileStart = 33;
          fileStop = 36;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::nineth";
          break;
        case OneFourthChainPartId::tenth:
          fileStart = 37;
          fileStop = 40;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::tenth";
          break;
        case OneFourthChainPartId::last:
          fileStart = 41;
          fileStop = 44;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::last";
          break;
        default:
          fileStart = 1;
          fileStop = 44;
          LOG(info) << "run " << runN
                    << " , OneFourthChainPartId::all";
          break;
      }
      foundPath = true;
    } else {
      LOG(error) << "run " << runN
                 << " , useNewCTFs = " << useNewCTFs
                 << " , alignStatus = " << int(alignStatus)
                 << " , condition not found";
      foundPath = false;
    }
  }

  if (foundPath) {
    std::stringstream generalPathSs;
    generalPathSs << basePath << "/" << runN << "/" << alignStatusPath;
    generalPath = generalPathSs.str();
    std::filesystem::directory_entry entry{generalPath};
    if (entry.exists()) {
      foundPath = true;
      LOG(info) << "Found requested path : " << generalPath;
    } else {
      foundPath = false;
      LOG(error) << "Path does not exist : " << generalPath;
    }
  }

  if (!foundPath) {
    success = false;
    LOG(error) << "addFilesToChains() - aborted !";
    return success;
  }

  success = false;
  if (mychoice == ChainBuildMethod::continuous) {

    LOG(info) << "ChainBuildMethod::continuous"
              << " , start " << fileStart
              << " , stop " << fileStop;

    int countFiles = 0;
    for (int ii = fileStart; ii <= fileStop; ii++) {
      std::stringstream ss;

      if (runN == 535046) {
        ss << generalPath << "/" << ii;
      } else {
        if (ii < 100) {
          ss << generalPath << "/"
             << std::setw(3) << std::setfill('0') << ii;
        } else {
          ss << generalPath << "/" << ii;
        }
      }
      std::string filePath = ss.str();
      if (std::filesystem::exists(Form("%s/mftclusters.root", filePath.c_str())) &&
          std::filesystem::exists(Form("%s/mfttracks.root", filePath.c_str()))) {
        mftclusterChain->Add(Form("%s/mftclusters.root", filePath.c_str()));
        mfttrackChain->Add(Form("%s/mfttracks.root", filePath.c_str()));
        LOG(info) << "Add " << Form("%s/mftclusters.root", filePath.c_str());
        LOG(info) << "Add " << Form("%s/mfttracks.root", filePath.c_str());
        countFiles++;
      } else {
        LOG(error) << "mftclusters.root or mfttracks.root not found at : " << filePath;
      }
    }
    LOG(info) << "number of files per chain = " << countFiles;
    if (countFiles) {
      success = true;
    }
  } else if (mychoice == ChainBuildMethod::gridFailed) {

    LOG(info) << "ChainBuildMethod::gridFailed";
    if (runN == 505720 && useNewCTFs && alignStatus == AlignmentStatus::idealgeo) {
      static constexpr int nFiles = 13;
      static constexpr std::array<int, nFiles> gridSubJob{
        33, 538, 1577, 1578, 1579,
        1580, 1581, 1582, 1583, 1604,
        1554, 1563, 1585};

      int idFailed = 0;
      int countFiles = 0;
      for (int ii = fileStart; ii <= fileStop; ii++) {
        std::stringstream ss;
        if (ii < 100) {
          ss << generalPath << "/"
             << std::setw(3) << std::setfill('0') << ii;
        } else {
          ss << generalPath << "/" << ii;
        }
        std::string filePath = ss.str();
        if (ii == gridSubJob[idFailed]) {
          idFailed++;
          LOG(info) << "--- Skip " << Form("%s/mftclusters.root", filePath.c_str());
          LOG(info) << "--- Skip " << Form("%s/mfttracks.root", filePath.c_str());
          continue;
        }
        if (std::filesystem::exists(Form("%s/mftclusters.root", filePath.c_str())) &&
            std::filesystem::exists(Form("%s/mfttracks.root", filePath.c_str()))) {
          mftclusterChain->Add(Form("%s/mftclusters.root", filePath.c_str()));
          mfttrackChain->Add(Form("%s/mfttracks.root", filePath.c_str()));
          LOG(info) << "Add " << Form("%s/mftclusters.root", filePath.c_str());
          LOG(info) << "Add " << Form("%s/mfttracks.root", filePath.c_str());
          countFiles++;
        }
      }
      LOG(info) << "number of files per chain = " << countFiles;
      if (countFiles) {
        success = true;
      }
    } else {
      LOG(error) << "ChainBuildMethod::gridFailed could not be used";
    }
  } else {

    assert(("Should'nt be executed when mychoice != ChainBuildMethod::gridOk",
            mychoice == ChainBuildMethod::gridOk));
    LOG(info) << "ChainBuildMethod::gridOk";
    if (runN == 505713 && useNewCTFs && alignStatus == AlignmentStatus::idealgeo) {
      static constexpr int nFiles = 763;
      static constexpr std::array<int, nFiles> gridSubJob{
        // first array of 587 successful jobs done on the grid
        1, 2, 4, 7, 10, 12, 13, 14, 16, 18, 19, 20, 22, 24, 26,
        28, 30, 31, 34, 38, 39, 41, 42, 44, 45, 46, 47, 49, 50, 51,
        52, 55, 56, 57, 59, 60, 62, 64, 65, 66, 67, 68, 69, 70, 76,
        77, 78, 79, 80, 82, 83, 84, 86, 87, 88, 89, 90, 92, 93, 94,
        95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 109, 110,
        111, 112, 113, 114, 116, 119, 120, 121, 122, 126, 127, 128, 130, 131, 133,
        134, 135, 136, 137, 138, 139, 140, 141, 142, 144, 145, 146, 147, 148, 149,
        150, 151, 152, 153, 154, 155, 156, 157, 158, 160, 161, 162, 164, 165, 166,
        167, 169, 171, 172, 173, 174, 176, 179, 180, 182, 184, 185, 186, 188, 189,
        190, 191, 192, 193, 195, 196, 197, 198, 199, 201, 202, 203, 204, 206, 207,
        208, 209, 211, 212, 213, 214, 215, 216, 218, 219, 221, 222, 223, 224, 225,
        226, 227, 228, 230, 231, 232, 233, 235, 236, 238, 239, 241, 242, 243, 244,
        245, 246, 247, 248, 250, 251, 252, 253, 256, 258, 259, 260, 261, 263, 264,
        265, 268, 269, 270, 271, 272, 273, 276, 277, 278, 281, 282, 284, 285, 286,
        288, 289, 290, 292, 293, 294, 296, 297, 298, 299, 301, 302, 303, 304, 306,
        307, 308, 309, 310, 311, 312, 314, 315, 316, 317, 318, 319, 320, 322, 323,
        324, 325, 326, 327, 328, 329, 330, 332, 333, 335, 336, 337, 338, 339, 341,
        342, 343, 345, 347, 348, 349, 352, 353, 354, 355, 357, 359, 360, 362, 364,
        367, 368, 372, 373, 374, 375, 378, 379, 380, 381, 382, 383, 384, 385, 387,
        390, 392, 397, 399, 400, 403, 404, 405, 408, 409, 412, 414, 417, 418, 420,
        421, 422, 426, 427, 430, 432, 433, 438, 441, 442, 443, 444, 446, 447, 448,
        449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 464,
        465, 466, 467, 468, 469, 471, 472, 474, 475, 477, 478, 480, 482, 483, 485,
        486, 487, 488, 490, 491, 492, 493, 494, 496, 497, 498, 499, 500, 501, 502,
        503, 504, 505, 507, 509, 510, 512, 513, 514, 515, 516, 517, 518, 520, 523,
        524, 525, 526, 527, 528, 529, 531, 532, 533, 534, 535, 536, 537, 538, 539,
        540, 541, 543, 544, 545, 549, 550, 552, 553, 554, 555, 556, 557, 559, 561,
        562, 563, 570, 571, 572, 573, 574, 577, 578, 581, 582, 583, 586, 587, 588,
        590, 593, 594, 597, 598, 599, 600, 602, 603, 605, 608, 609, 610, 615, 616,
        617, 618, 621, 622, 623, 624, 625, 627, 628, 629, 631, 632, 633, 634, 636,
        637, 638, 640, 642, 643, 644, 647, 650, 651, 653, 654, 657, 660, 664, 666,
        669, 670, 671, 673, 674, 677, 678, 679, 680, 681, 683, 684, 685, 686, 687,
        688, 689, 691, 692, 693, 695, 696, 698, 699, 700, 701, 702, 703, 704, 705,
        706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 717, 718, 719, 720, 721,
        725, 727, 728, 729, 730, 731, 732, 734, 735, 736, 738, 739, 740, 742, 744,
        745, 746, 747, 748, 749, 750, 751, 753, 754, 755, 756, 757, 758, 760, 763,
        766, 768, 769, 770, 772, 773, 777, 778, 779, 780, 781, 782, 783, 784, 786,
        788, 789, 790, 791, 793, 794, 796, 797, 798, 799, 800, 801, 802, 803, 804,
        805, 807, 808, 810, 811, 812, 816, 817, 818, 819, 823, 824, 826, 829, 830,
        832, 833,
        // second array of 176 successful jobs done on the grid
        5, 8, 9, 11, 17, 21, 23, 25, 33, 35, 36, 37, 43, 48, 53,
        54, 58, 61, 63, 71, 72, 73, 74, 75, 81, 85, 115, 117, 123, 124,
        125, 143, 163, 168, 170, 175, 177, 178, 181, 183, 187, 194, 200, 220, 229,
        237, 240, 249, 254, 266, 267, 274, 275, 280, 283, 287, 291, 305, 313, 321,
        331, 334, 340, 344, 351, 356, 358, 361, 363, 365, 366, 369, 370, 371, 376,
        377, 386, 388, 389, 395, 396, 401, 402, 407, 411, 413, 415, 416, 419, 423,
        424, 425, 428, 431, 434, 435, 436, 437, 440, 445, 463, 476, 481, 489, 508,
        511, 521, 530, 542, 546, 547, 551, 558, 560, 564, 565, 566, 568, 569, 579,
        580, 589, 591, 592, 595, 601, 606, 607, 611, 612, 613, 614, 619, 620, 626,
        630, 635, 639, 641, 645, 646, 648, 649, 652, 656, 659, 661, 662, 663, 665,
        667, 675, 676, 682, 690, 694, 697, 722, 723, 724, 726, 743, 759, 761, 764,
        767, 775, 776, 787, 792, 795, 806, 809, 814, 815, 831};

      int countFiles = 0;
      for (const auto& ii : gridSubJob) {
        if (countFiles > fileStop) {
          break;
        }
        std::stringstream ss;
        if (ii < 100) {
          ss << generalPath << "/"
             << std::setw(3) << std::setfill('0') << ii;
        } else {
          ss << generalPath << "/" << ii;
        }
        std::string filePath = ss.str();
        if (std::filesystem::exists(Form("%s/mftclusters.root", filePath.c_str())) &&
            std::filesystem::exists(Form("%s/mfttracks.root", filePath.c_str()))) {
          mftclusterChain->Add(Form("%s/mftclusters.root", filePath.c_str()));
          mfttrackChain->Add(Form("%s/mfttracks.root", filePath.c_str()));
          LOG(info) << "Add " << Form("%s/mftclusters.root", filePath.c_str());
          LOG(info) << "Add " << Form("%s/mfttracks.root", filePath.c_str());
          countFiles++;
        }
      }
      LOG(info) << "number of files per chain = " << countFiles;
      if (countFiles) {
        success = true;
      }
    } else {
      LOG(error) << "ChainBuildMethod::gridOk could not be used";
    }
  }

  LOG(info) << "addFilesToChains() - successful for run " << runN;
  return success;
}

//_________________________________
// alienv setenv O2Physics/latest -c root -l
// .L ~/cernbox/alice/enigma/macros/runTracksToRecords.C++
// runTracksToRecords()
//_________________________________
void runTracksToRecords(const OneFourthChainPartId partId = OneFourthChainPartId::first,
                        const int runN = 535046,
                        const AlignmentStatus alignStatus = AlignmentStatus::dcacorrpass2rotated,
                        const int minPoints = 6,
                        const bool useNewCTFs = true,
                        const bool doControl = true,
                        const int nEntriesAutoSave = 10000)
{

  ROOT::EnableImplicitMT(0);
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

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

  LOG(info) << "applyMisalignment = " << applyMisalignment;
  LOG(info) << "preferAlignedFile = " << preferAlignedFile;
  o2::base::GeometryManager::loadGeometry("", applyMisalignment, preferAlignedFile);
  o2::mft::GeometryTGeo* geom = o2::mft::GeometryTGeo::Instance();
  geom->fillMatrixCache(
    o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                             o2::math_utils::TransformType::L2G));

  // dictionary

  std::string dictFileName = "MFTdictionary.bin";
  if (useNewCTFs) {
    dictFileName = "o2_mft_dictionary.root";
  }
  o2::itsmft::TopologyDictionary* dict = nullptr;
  try {
    dict = o2::itsmft::TopologyDictionary::loadFrom(dictFileName);
  } catch (std::exception e) {
    std::cout << "Error " << e.what() << std::endl;
    LOG(error) << "runTracksToRecords() - aborted !";
    return;
  }
  dict->readFromFile(dictFileName);

  // cluster and track chains

  TChain* mftclusterChain = new TChain("o2sim");
  TChain* mfttrackChain = new TChain("o2sim");

  // feed cluster and track chains

  bool success = addFilesToChains(mftclusterChain,
                                  mfttrackChain,
                                  runN,
                                  alignStatus,
                                  useNewCTFs,
                                  partId);

  if (!success) {
    LOG(error) << "runTracksToRecords() - aborted !";
    return;
  }

  // instantiate and configure the aligner

  AlignConfigHelper alignConfigParam;
  alignConfigParam.minPoints = minPoints;

  o2::mft::TracksToRecords aligner;

  aligner.setRunNumber(runN);
  aligner.setBz(0.);

  aligner.setClusterDictionary(dict);
  aligner.setMinNumberClusterCut(alignConfigParam.minPoints);

  aligner.setWithControl(doControl);
  aligner.setNEntriesAutoSave(nEntriesAutoSave);

  // TODO: fix det. elements here

  // init Millipede

  aligner.init();

  // compute Mille records

  aligner.startRecordWriter();
  aligner.processROFs(mfttrackChain, mftclusterChain);
  aligner.printProcessTrackSummary();
  aligner.endRecordWriter();

  // the end

  std::chrono::steady_clock::time_point stop_time = std::chrono::steady_clock::now();

  std::cout << "----------------------------------- " << endl;
  std::cout << "Total Execution time: \t\t"
            << std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time).count()
            << " seconds" << endl;
}
