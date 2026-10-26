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

  std::array<TMultiGraph, 6> mgr;
  mgr[0].SetTitle("dx difference");
  mgr[1].SetTitle("dy difference");
  mgr[2].SetTitle("dz difference");
  mgr[3].SetTitle("dRx difference");
  mgr[4].SetTitle("dRy difference");
  mgr[5].SetTitle("dRz difference");

  //std::array<int, 10> colors{kRed - 9, kCyan - 9, kRed - 6, kCyan - 6, kRed - 2, kCyan - 2, kRed + 3, kCyan +3, kRed + 4, kCyan + 4};
  std::array<int, 10> colors{kRed + 3, kCyan + 3, kRed + 2, kCyan + 2, kRed + 1, kCyan + 1, kRed + 0, kCyan + 0, kRed - 9, kCyan - 9};
  std::array<std::array<TGraph, 10>, 6> grLayer;
  for (int l = 0;l < 10; l++) {
    grLayer[0][l].SetTitle(std::format("dx difference - L{}", l).c_str());
    grLayer[0][l].SetMarkerColor(colors[l]);
    mgr[0].Add(&(grLayer[0][l]));

    grLayer[1][l].SetTitle(std::format("dy difference - L{}", l).c_str());
    grLayer[1][l].SetMarkerColor(colors[l]);
    mgr[1].Add(&(grLayer[1][l]));

    grLayer[2][l].SetTitle(std::format("dz difference - L{}", l).c_str());
    grLayer[2][l].SetMarkerColor(colors[l]);
    mgr[2].Add(&(grLayer[2][l]));

    grLayer[3][l].SetTitle(std::format("dRx difference - L{}", l).c_str());
    grLayer[3][l].SetMarkerColor(colors[l]);
    mgr[3].Add(&(grLayer[3][l]));

    grLayer[4][l].SetTitle(std::format("dRy difference - L{}", l).c_str());
    grLayer[4][l].SetMarkerColor(colors[l]);
    mgr[4].Add(&(grLayer[4][l]));

    grLayer[5][l].SetTitle(std::format("dRz difference - L{}", l).c_str());
    grLayer[5][l].SetMarkerColor(colors[l]);
    mgr[5].Add(&(grLayer[5][l]));
  }

  for (size_t i = 0; i < mftGeom[0].size(); i++) {
    int chipId1 = mftGeom[0][i].chipId;
    int chipId2 = mftGeom[1][i].chipId;
    if (chipId1 != chipId2) {
      std::cout << std::format("Chip ID mistatch for element {}: {} / {}", i, chipId1, chipId2) << std::endl;
    }
    int l = mftGeom[0][i].layer;
    
    if (chipId1 < 0 || l < 0) continue;

    //std::cout << std::format("Chip = {}  Layer = {}", chipId1, l) << std::endl;

    if (chipId1 == 930) {
      std::cout << std::format("Chip = {}  Layer = {}  dRz = {} / {}", chipId1, l, mftGeom[0][i].dRz, mftGeom[1][i].dRz) << std::endl;
    }

    gr[0].AddPoint(chipId1, mftGeom[1][i].dx - mftGeom[0][i].dx);
    grLayer[0][l].AddPoint(chipId1, mftGeom[1][i].dx - mftGeom[0][i].dx);
    gr[1].AddPoint(chipId1, mftGeom[1][i].dy - mftGeom[0][i].dy);
    grLayer[1][l].AddPoint(chipId1, mftGeom[1][i].dy - mftGeom[0][i].dy);
    gr[2].AddPoint(chipId1, mftGeom[1][i].dz - mftGeom[0][i].dz);
    grLayer[2][l].AddPoint(chipId1, mftGeom[1][i].dz - mftGeom[0][i].dz);
    gr[3].AddPoint(chipId1, mftGeom[1][i].dRx - mftGeom[0][i].dRx);
    grLayer[3][l].AddPoint(chipId1, mftGeom[1][i].dRx - mftGeom[0][i].dRx);
    gr[4].AddPoint(chipId1, mftGeom[1][i].dRy - mftGeom[0][i].dRy);
    grLayer[4][l].AddPoint(chipId1, mftGeom[1][i].dRy - mftGeom[0][i].dRy);
    gr[5].AddPoint(chipId1, mftGeom[1][i].dRz - mftGeom[0][i].dRz);
    grLayer[5][l].AddPoint(chipId1, mftGeom[1][i].dRz - mftGeom[0][i].dRz);
  }
  
  TPaveText* title = new TPaveText(0.1, 0.4, 0.9, 0.6, "NDC");
  title->Clear(); title->AddText("Delta x"); title->Draw();
  c.SaveAs("mftGeomComp.pdf(");

  TLine* l;
  int topChipIdMin = 468;

  mgr[0].Draw("A*");
  l = new TLine(topChipIdMin, mgr[0].GetYaxis()->GetXmin(), topChipIdMin, mgr[0].GetYaxis()->GetXmax()); l->SetLineStyle(kDashed); l->Draw();
  c.SaveAs("mftGeomComp.pdf");
  mgr[0].SetMinimum(-0.1);
  mgr[0].SetMaximum(0.1);
  l = new TLine(topChipIdMin, -0.1, topChipIdMin, 0.1); l->SetLineStyle(kDashed); l->Draw();
  c.SaveAs("mftGeomComp.pdf");
  mgr[0].SetMinimum(-0.01);
  mgr[0].SetMaximum(0.01);
  l = new TLine(topChipIdMin, -0.01, topChipIdMin,0.01); l->SetLineStyle(kDashed); l->Draw();
  c.SaveAs("mftGeomComp.pdf");
  
  c.Clear(); title->Clear(); title->AddText("Delta y"); title->Draw();
  c.SaveAs("mftGeomComp.pdf");

  mgr[1].Draw("A*");
  l = new TLine(topChipIdMin, mgr[0].GetYaxis()->GetXmin(), topChipIdMin, mgr[0].GetYaxis()->GetXmax()); l->SetLineStyle(kDashed); l->Draw();
  c.SaveAs("mftGeomComp.pdf");
  mgr[1].SetMinimum(-0.1);
  mgr[1].SetMaximum(0.1);
  l = new TLine(topChipIdMin, -0.1, topChipIdMin, 0.1); l->SetLineStyle(kDashed); l->Draw();
  c.SaveAs("mftGeomComp.pdf");
  mgr[1].SetMinimum(-0.01);
  mgr[1].SetMaximum(0.01);
  l = new TLine(topChipIdMin, -0.01, topChipIdMin,0.01); l->SetLineStyle(kDashed); l->Draw();
  c.SaveAs("mftGeomComp.pdf");
  
  c.Clear(); title->Clear(); title->AddText("Delta z"); title->Draw();
  c.SaveAs("mftGeomComp.pdf");

  mgr[2].Draw("A*");
  l = new TLine(topChipIdMin, mgr[0].GetYaxis()->GetXmin(), topChipIdMin, mgr[0].GetYaxis()->GetXmax()); l->SetLineStyle(kDashed); l->Draw();
  c.SaveAs("mftGeomComp.pdf");
  mgr[2].GetXaxis()->SetRangeUser(835, 845);
  c.SaveAs("mftGeomComp.pdf");
  mgr[2].GetXaxis()->SetRange(0, -1);
  mgr[2].SetMinimum(-0.1);
  mgr[2].SetMaximum(0.1);
  l = new TLine(topChipIdMin, -0.1, topChipIdMin, 0.1); l->SetLineStyle(kDashed); l->Draw();
  c.SaveAs("mftGeomComp.pdf");
  mgr[2].SetMinimum(-0.01);
  mgr[2].SetMaximum(0.01);
  l = new TLine(topChipIdMin, -0.01, topChipIdMin,0.01); l->SetLineStyle(kDashed); l->Draw();
  c.SaveAs("mftGeomComp.pdf");

  c.Clear(); title->Clear(); title->AddText("Delta Rx"); title->Draw();
  c.SaveAs("mftGeomComp.pdf");

  mgr[3].Draw("A*");
  l = new TLine(topChipIdMin, mgr[0].GetYaxis()->GetXmin(), topChipIdMin, mgr[0].GetYaxis()->GetXmax()); l->SetLineStyle(kDashed); l->Draw();
  c.SaveAs("mftGeomComp.pdf");
  mgr[3].SetMinimum(-1.e-5);
  mgr[3].SetMaximum(1.e-5);
  c.SaveAs("mftGeomComp.pdf");
  
  c.Clear(); title->Clear(); title->AddText("Delta Ry"); title->Draw();
  c.SaveAs("mftGeomComp.pdf");

  mgr[4].Draw("A*");
  l = new TLine(topChipIdMin, mgr[0].GetYaxis()->GetXmin(), topChipIdMin, mgr[0].GetYaxis()->GetXmax()); l->SetLineStyle(kDashed); l->Draw();
  c.SaveAs("mftGeomComp.pdf");
  mgr[4].SetMinimum(-1.e-5);
  mgr[4].SetMaximum(1.e-5);
  c.SaveAs("mftGeomComp.pdf");
  
  c.Clear(); title->Clear(); title->AddText("Delta Rz"); title->Draw();
  c.SaveAs("mftGeomComp.pdf");

  mgr[5].Draw("A*");
  l = new TLine(topChipIdMin, mgr[0].GetYaxis()->GetXmin(), topChipIdMin, mgr[0].GetYaxis()->GetXmax()); l->SetLineStyle(kDashed); l->Draw();
  c.SaveAs("mftGeomComp.pdf");
  mgr[5].SetMinimum(-1.);
  mgr[5].SetMaximum(1.);
  c.SaveAs("mftGeomComp.pdf");
  mgr[5].SetMinimum(-1.e-1);
  mgr[5].SetMaximum(1.e-1);
  c.SaveAs("mftGeomComp.pdf");
  mgr[5].SetMinimum(-0.02);
  mgr[5].SetMaximum(0.02);
  c.SaveAs("mftGeomComp.pdf");
  mgr[5].SetMinimum(-1.e-5);
  mgr[5].SetMaximum(1.e-5);
  c.SaveAs("mftGeomComp.pdf)");
}
