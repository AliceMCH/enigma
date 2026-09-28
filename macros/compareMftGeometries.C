struct MftChipGeom
{
  int half{-1};
  int disk{-1};
  int layer{-1};
  int zone{-1};
  int con{-1};
  int tr{-1};
  int chipId{-1};
  float dx{0};
  float dy{0};
  float dz{0};
  float dRx{0};
  float dRy{0};
  float dRz{0};
};

void compareMftGeometries(std::string csvGeom1, std::string csvGeom2)
{
  std::array<std::ifstream, 2> csvStr;
  csvStr[0].open(csvGeom1);
  csvStr[1].open(csvGeom2);
  std::array<std::vector<MftChipGeom>, 2> mftGeom;
  
  TCanvas c("c", "c", 1200, 800);
  
  for (int i = 0; i < 2; i++) {
    std::string line;
    while (std::getline(csvStr[i], line)) {
      // replace commas with spaces
      std::replace( line.begin(), line.end(), ',', ' ');
      std::stringstream ss(line);
      
      std::cout << line << std::endl;
      
      mftGeom[i].emplace_back(MftChipGeom{});
      MftChipGeom& chip = mftGeom[i].back();

      // half,disk,layer,zone,con,tr,chipid,dx,dy,dz,dRx,dRy,dRz
      ss >> chip.half >> chip.disk >> chip.layer >> chip.zone
      >> chip.con >> chip.tr >> chip.chipId
      >> chip.dx >> chip.dy >> chip.dz
      >> chip.dRx >> chip.dRy >> chip.dRz;
    }
  }
  
  std::cout << std::format("mftGeom sizes: {} {}", mftGeom[0].size(), mftGeom[1].size()) << std::endl;
  
  if (mftGeom[0].size() != mftGeom[1].size()) {
    return;
  }
  
  std::array<TGraph, 6> gr;
  gr[0].SetTitle("dx difference");
  gr[1].SetTitle("dy difference");
  gr[2].SetTitle("dz difference");
  gr[3].SetTitle("dRx difference");
  gr[4].SetTitle("dRy difference");
  gr[5].SetTitle("dRz difference");
  for (size_t i = 0; i < mftGeom[0].size(); i++) {
    int chipId1 = mftGeom[0][i].chipId;
    int chipId2 = mftGeom[1][i].chipId;
    if (chipId1 != chipId2) {
      std::cout << std::format("Chip ID mistatch for element {}: {} / {}", i, chipId1, chipId2) << std::endl;
    }
    
    gr[0].AddPoint(chipId1, mftGeom[1][i].dx - mftGeom[0][i].dx);
    gr[1].AddPoint(chipId1, mftGeom[1][i].dy - mftGeom[0][i].dy);
    gr[2].AddPoint(chipId1, mftGeom[1][i].dz - mftGeom[0][i].dz);
    gr[3].AddPoint(chipId1, mftGeom[1][i].dRx - mftGeom[0][i].dRx);
    gr[4].AddPoint(chipId1, mftGeom[1][i].dRy - mftGeom[0][i].dRy);
    gr[5].AddPoint(chipId1, mftGeom[1][i].dRz - mftGeom[0][i].dRz);
  }
  
  gr[0].Draw("A*");
  c.SaveAs("mftGeomComp.pdf(");
  gr[0].SetMinimum(-0.1);
  gr[0].SetMaximum(0.1);
  c.SaveAs("mftGeomComp.pdf");
  gr[0].SetMinimum(-0.01);
  gr[0].SetMaximum(0.01);
  c.SaveAs("mftGeomComp.pdf");
  
  gr[1].Draw("A*");
  c.SaveAs("mftGeomComp.pdf");
  gr[1].SetMinimum(-0.1);
  gr[1].SetMaximum(0.1);
  c.SaveAs("mftGeomComp.pdf");
  gr[1].SetMinimum(-0.01);
  gr[1].SetMaximum(0.01);
  c.SaveAs("mftGeomComp.pdf");
  
  gr[2].Draw("A*");
  c.SaveAs("mftGeomComp.pdf");
  gr[2].SetMinimum(-0.1);
  gr[2].SetMaximum(0.1);
  c.SaveAs("mftGeomComp.pdf");
  gr[2].SetMinimum(-0.01);
  gr[2].SetMaximum(0.01);
  c.SaveAs("mftGeomComp.pdf");

  gr[3].Draw("A*");
  c.SaveAs("mftGeomComp.pdf");
  gr[3].SetMinimum(-1.e-5);
  gr[3].SetMaximum(1.e-5);
  c.SaveAs("mftGeomComp.pdf");
  
  gr[4].Draw("A*");
  c.SaveAs("mftGeomComp.pdf");
  gr[4].SetMinimum(-1.e-5);
  gr[4].SetMaximum(1.e-5);
  c.SaveAs("mftGeomComp.pdf");
  
  gr[5].Draw("A*");
  c.SaveAs("mftGeomComp.pdf");
  gr[5].SetMinimum(-1.e-5);
  gr[5].SetMaximum(1.e-5);
  c.SaveAs("mftGeomComp.pdf)");
}