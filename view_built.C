void view_built(string fname, int nwf=300, int schan=0, int nchan=1,
		bool wait=false, bool rndm=false, bool same=false){
  gStyle->SetOptStat(0);
  gStyle->SetLabelFont(132, "XYZ");
  gStyle->SetTitleFont(132, "XYZ");
  gRandom->SetSeed(0);
  int iwf = 0;
  TFile* file = TFile::Open(fname.c_str());
  TTree* tree = (TTree*) file->Get("MGTree");
  MGTEvent* event = new MGTEvent();
  tree->SetBranchAddress("event", &event);
  int dx = std::min(3, nchan);
  int dy = 1+(nchan/3);
  TCanvas* c = new TCanvas("c", "", 0, 0, 800*dx, 600*dy);
  c->Divide(dx, dy);
  vector<TH1D*> vh;
  for(int iev=0; iev<(int)tree->GetEntries(); iev++){
    int jev = iev;
    if(rndm)
      jev = (int) (gRandom->Rndm() * tree->GetEntries());
    tree->GetEvent(jev);
    vh.assign(event->GetNWaveforms(), NULL);
    if(vh.size() == 0) continue;
    double maxval = -1.0e9;
    double minval =  1.0e9;
    for(int ich=0; ich<(int)event->GetNWaveforms(); ich++){
      vh[ich] = event->GetWaveform(ich)->GimmeUniqueHist();
      vh[ich]->SetTitle("");
      vh[ich]->SetLineColor(4);
      maxval = std::max(vh[ich]->GetMaximum(), maxval);
      minval = std::min(vh[ich]->GetMinimum(), minval);
      iwf ++;
    }
    vh[0]->GetYaxis()->SetRangeUser(minval-0.05*(maxval-minval),
    				    maxval+0.05*(maxval-minval));
    for(int ich=0; ich<(int)event->GetNWaveforms(); ich++){
      if(same){
	c->cd(1);
	vh[ich]->SetLineColor(ich+1);
	if(ich == 0) vh[ich]->Draw();
	else vh[ich]->Draw("same");
      }
      else{
	c->cd(ich+1);
	vh[ich]->Draw();
      }
    }
    c->Update();
    if(wait) getchar();
    for(int i=0; i<(int)vh.size(); i++) delete vh[i];
    if(iwf >= nwf - 1)
      break;
  }
  file->Close();
}
