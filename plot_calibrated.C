void plot_calibrated(string infname, int index=0, string outfname=""){

  gStyle->SetOptStat(0);
  gStyle->SetLabelFont(132, "XYZ");
  gStyle->SetTitleFont(132, "XYZ");

  TChain* tree = new TChain("tree");
  tree->Add(infname.c_str());
  TTreeReader reader(tree);
  TTreeReaderValue<vector<int> > channel(reader, "channel");
  TTreeReaderValue<vector<double> > baseSigma(reader, "baseSigma");
  TTreeReaderValue<vector<double> > brmsSigma(reader, "brmsSigma");
  TTreeReaderValue<vector<double> > time(reader, "time");
  TTreeReaderValue<vector<double> > deltat(reader, "deltat");
  TTreeReaderValue<vector<double> > t0(reader, "t0");
  TTreeReaderValue<vector<double> > trise(reader, "trise");
  TTreeReaderValue<vector<double> > trapECal(reader, "trapECal");
  TTreeReaderValue<vector<double> > trapEFCal(reader, "trapEFCal");
  TTreeReaderValue<vector<double> > trapEFRCal(reader, "trapEFRCal");
  TTreeReaderValue<vector<double> > avse(reader, "avse");
  TTreeReaderValue<vector<double> > dcr(reader, "dcr");

  TH1D* he = new TH1D("he", "", 5000, 0.0, 5000.0);
  he->SetXTitle("Energy (keV)");
  he->SetYTitle("Entries");
  he->SetLineColor(2);
  TH1D* hef = (TH1D*) he->Clone("hef");
  hef->SetLineColor(1);
  TH1D* hefa = (TH1D*) he->Clone("hefa");
  hefa->SetLineColor(8);
  TH1D* hefd = (TH1D*) he->Clone("hefd");
  hefd->SetLineColor(4);
  TH1D* hefr = (TH1D*) he->Clone("hefr");
  hefr->SetLineColor(6);
  
  TH2D* havse = new TH2D("havse", "", 1000, 0.0, 5000.0, 400, -20.0, 20.0);
  havse->SetXTitle("Energy (keV)");
  havse->SetYTitle("AvsE");
  TH2D* hdcre = new TH2D("hdcre", "", 1000, 0.0, 5000.0, 400, -2.0, 2.0);
  hdcre->SetXTitle("Energy (keV)");
  hdcre->SetYTitle("DCR");
  
  while(reader.Next()){
    if(index > channel->size()) continue;
    if(abs(baseSigma->at(index)) > 3) continue;
    if(abs(brmsSigma->at(index)) > 3) continue;
    he->Fill(trapECal->at(index));
    hef->Fill(trapEFCal->at(index));
    havse->Fill(trapEFCal->at(index), avse->at(index));
    if(avse->at(index) > -1.0){
      hefa->Fill(trapEFCal->at(index));
      hdcre->Fill(trapEFCal->at(index), dcr->at(index));
      if(dcr->at(index) >= -1.0 && dcr->at(index) < 1.0){
	hefd->Fill(trapEFCal->at(index));
	hefr->Fill(trapEFRCal->at(index));
      }
    }
  }

  TCanvas* ce = new TCanvas("ce", "", 0, 0, 1200, 800);
  ce->cd(1)->SetLogy(true);
  ce->cd(1);
  he->Draw();
  hefr->Draw("same");
  hef->Draw("same");
  hefa->Draw("same");
  hefd->Draw("same");
  TCanvas* cavse = new TCanvas("cavse", "", 400, 0, 1200, 800);
  cavse->cd(1)->SetLogz(true);
  cavse->cd(1);
  havse->Draw("colz");
  cavse->Update();
  TCanvas* cdcre = new TCanvas("cdcre", "", 0, 400, 1200, 800);
  cdcre->cd(1)->SetLogz(true);
  cdcre->cd(1);
  hdcre->Draw("colz");

  if(outfname != ""){
    getchar();
    TFile* outfile = new TFile(outfname.c_str(), "recreate");
    outfile->cd();
    he->Write();
    hef->Write();
    hefr->Write();
    havse->Write();
    ce->Write();
    cavse->Write();
    outfile->Close();
  }

}
  
  
