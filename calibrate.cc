#include "utils.hh"
#include <MGTVector.hh>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <json/json.h>
#include <json/value.h>
#include <json/reader.h>
#include <json/writer.h>
#include <utility>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <assert.h>

using namespace std;

TDirectory* tdir;

// convenience function for printing calibration parameters
void print_value(string name, int prec, double val, double uncert){
  cout << "  " << name << ": " << setprecision(prec) << val
       << "+/-" << setprecision(prec) << uncert << endl;
}

// just a placeholder for a better algorithm
void Th_scale(TH1D* h, double& offset, double& scale){
  double max0 = 0.0, max1 = 0.0;
  double maxbin0 = 0, maxbin1 = 0;
  for(int bin=1; bin<=(int)h->GetNbinsX(); bin++){
    double val = h->GetBinContent(bin);
    double adc = h->GetBinCenter(bin);
    if(val > max0){
      max0 = val;
      maxbin0 = adc;
    }
    if(val > max1 && adc > 1200){
      max1 = val;
      maxbin1 = adc;
    }
  }
  scale = (2615-238) / (maxbin1-maxbin0);
  offset = 2615 - scale*maxbin1;
}

// should be identical as above, but use the Co peaks
void Co_scale(TH1D* h, double& offset, double& scale){
  h->GetNbinsX();
  offset = 0.0;
  scale = 0.0;
}

int main(int argc, char* argv[]){

  // root setup
  gStyle->SetOptStat(0);
  gStyle->SetLabelFont(132, "XYZ");
  gStyle->SetTitleFont(132, "XYZ");
  tdir = gROOT->CurrentDirectory();

  // Tl peaks to use for the rise time correction
  const vector<double> tl_peaks({338.0, 510.5, 583.0, 726.5,
	794.0, 860.0, 911.0, 968.0, 2615.0});

  //option handling
  map<int, int> chan_map;
  vector<string> ifname;
  string ofname = "";
  string source = "";
  bool json_config = false;
  string jfile = "calibrate.json";
  Json::Value jvalue;
  bool use_fixedt = true;
  static struct option opts[]{
    {"help",             no_argument, NULL, 'h'},
    {"channel",    required_argument, NULL, 'c'},
    {"infile",     required_argument, NULL, 'i'},
    {"outfile",    required_argument, NULL, 'o'},
    {"source",     required_argument, NULL, 's'},
    {"jsonconfig", required_argument, NULL, 'j'},
    {"jsonfile",   required_argument, NULL, 'J'},
    {"fixedtime",        no_argument, NULL, 'f'}
  };
  int opt = getopt_long(argc, argv, "hc:i:o:s:j:J:f", opts, NULL);
  while(opt != -1){
    switch(opt){
    case 'h':
      cout << "options:"                               << endl;
      cout << "  -c channel to be analyzed"            << endl;
      cout << "  -i input filename"                    << endl;
      cout << "  -o output filename"                   << endl;
      cout << "  -s source type (Th or Co)"            << endl;
      cout << "  -j name of input json config file"    << endl;
      cout << "  -J name of output json config file"   << endl;
      cout << "  -f do not use the fixed time pickoff" << endl;
    case 'c': chan_map[atoi(optarg)] = (int) chan_map.size()-1; break;
    case 'i': ifname.push_back(string(optarg));                 break;
    case 'o': ofname = string(optarg);                          break;
    case 's': source = string(optarg);                          break;
    case 'j':{
      ifstream jfile(optarg);
      Json::CharReaderBuilder reader;
      string errors;
      if(!parseFromStream(reader, jfile, &jvalue, &errors)){
	cout << errors << endl;
	return 2;
      }
      json_config = true;
      break;
    }
    case 'J': jfile = string(optarg); break;
    case 'f': use_fixedt = false;     break;
    default: return 1;
    }
    opt = getopt_long(argc, argv, "hc:i:o:s:j:J:f", opts, NULL);
  }
  assert(chan_map.size()>0 && ofname != "" && ifname.size() > 0);
  assert(json_config || source == "Th" || source == "Co");
  
  // read in input files, setup input tree
  TChain* intree = new TChain("tree");
  vector<TH1D*> hdeltat(chan_map.size(), NULL);
  vector<TH1D*> henergy(chan_map.size(), NULL);
  vector<TH1D*> henergyf(chan_map.size(), NULL);
  vector<TH1D*> hbase(chan_map.size(), NULL);
  vector<TH1D*> hbrms(chan_map.size(), NULL);
  vector<TH1D*> hdecay(chan_map.size(), NULL);
  vector<TH2D*> hbase_energy(chan_map.size(), NULL);
  vector<TH2D*> hbrms_energy(chan_map.size(), NULL);
  vector<TH2D*> hdecay_energy(chan_map.size(), NULL);
  vector<TH2D*> hdecay_deltat(chan_map.size(), NULL);
  vector<TH2D*> hamp_energy(chan_map.size(), NULL);
  vector<TH2D*> hamps_energy(chan_map.size(), NULL);
  vector<TH2D*> haoe_energy(chan_map.size(), NULL);
  vector<TH2D*> hdcr_energy(chan_map.size(), NULL);
  vector<TH2D*> hdcrs_energy(chan_map.size(), NULL);
  vector<TH2D*> hrise_energy(chan_map.size(), NULL);
  vector<TH2D*> hrises_energy(chan_map.size(), NULL);
  cout << "reading input histograms" << endl;
  for(auto const& n : ifname){
    intree->Add(n.c_str());
    TFile* f = TFile::Open(n.c_str());
    tdir->cd();
    for(auto const& p : chan_map){
      string s = "_" + to_string(p.first);
      int i = p.second;
      hdeltat[i]      = Get1DHistOrSum(f, "hdeltat"      +s, hdeltat[i]);
      henergy[i]      = Get1DHistOrSum(f, "henergy"      +s, henergy[i]);
      henergyf[i]     = Get1DHistOrSum(f, "henergyf"     +s, henergyf[i]);
      hbase_energy[i] = Get2DHistOrSum(f, "hbase_energy" +s, hbase_energy[i]);
      hbrms_energy[i] = Get2DHistOrSum(f, "hbrms_energy" +s, hbrms_energy[i]);
      hdecay_energy[i]= Get2DHistOrSum(f, "hdecay_energy"+s, hdecay_energy[i]);
      hdecay_deltat[i]= Get2DHistOrSum(f, "hdecay_deltat"+s, hdecay_deltat[i]);
      hamp_energy[i]  = Get2DHistOrSum(f, "hamp_energy"  +s, hamp_energy[i]);
      hdcr_energy[i]  = Get2DHistOrSum(f, "hdcr_energy"  +s, hdcr_energy[i]);
      hrise_energy[i] = Get2DHistOrSum(f, "hrise_energy" +s, hrise_energy[i]);
    }
    f->Close();
  }
  intree->SetBranchStatus("*", false);
  intree->SetBranchStatus("channel", true);
  intree->SetBranchStatus("imax", true);
  intree->SetBranchStatus("trappick", true);
  intree->SetBranchStatus("trapmax", true);
  intree->SetBranchStatus("dcrslope", true);
  TTreeReader reader(intree);
  TTreeReaderValue<vector<int> > channel(reader, "channel");
  TTreeReaderValue<vector<double> > trappick(reader, "trappick");
  TTreeReaderValue<vector<double> > trapmax(reader, "trapmax");
  TTreeReaderValue<vector<double> > imax(reader, "imax");
  TTreeReaderValue<vector<double> > dcrslope(reader, "dcrslope");
  TTreeReaderValue<vector<double> > t1(reader, "t1");
  TTreeReaderValue<vector<double> > t99(reader, "t99");
  
  // values and histograms to output
  vector<double> base_mean(chan_map.size(), 0.0);
  vector<double> brms_mean(chan_map.size(), 0.0);
  vector<double> base_rms(chan_map.size(), 0.0);
  vector<double> brms_rms(chan_map.size(), 0.0);
  vector<double> pz_mean(chan_map.size(), 0.0);
  vector<double> pz_rms(chan_map.size(), 0.0);
  vector<double> escale(chan_map.size(), 0.0);
  vector<double> efscale(chan_map.size(), 0.0);
  vector<double> eoffset(chan_map.size(), 0.0);
  vector<double> efoffset(chan_map.size(), 0.0);
  vector<TH1D*> hecal(chan_map.size(), NULL);
  vector<vector<double> > avse_param(chan_map.size(), vector<double>(4, 0.0));
  vector<vector<double> > avse_uncert(chan_map.size(), vector<double>(4, 0.0));
  vector<TGraphErrors*> gavse(chan_map.size(), NULL);
  vector<double> trise_val(chan_map.size(), 0.0);
  vector<double> trise_uncert(chan_map.size(), 0.0);
  vector<double> trise_slope(chan_map.size(), 0.0);
  vector<double> trise_suncert(chan_map.size(), 0.0);
  vector<TProfile*> hrises_energy_proj(chan_map.size(), NULL);
  vector<TGraphErrors*> hrise_slope(chan_map.size(), NULL);
  vector<double> dcre_slope(chan_map.size(), 0.0);
  vector<double> dcre_uncert(chan_map.size(), 0.0);
  vector<double> dcr_cut_lo(chan_map.size(), 0.0);
  vector<double> dcr_cut_hi(chan_map.size(), 0.0);

  // output tree setup
  vector<int> chan;
  vector<double> baseSigma, brmsSigma, timestamp, dt, T0, trise;
  vector<double> trapECal, trapEFCal, trapEFRCal, avse, aoe, dcr;
  TFile* outfile = new TFile(ofname.c_str(), "recreate");
  tdir->cd();
  TTree* outtree = new TTree("tree", "tree");
  outtree->Branch("channel", &chan);
  outtree->Branch("baseSigma", &baseSigma);
  outtree->Branch("brmsSigma", &brmsSigma);
  outtree->Branch("time", &timestamp);
  outtree->Branch("deltat", &dt);
  outtree->Branch("t0", &T0);
  outtree->Branch("trise", &trise);
  outtree->Branch("trapECal", &trapECal);
  outtree->Branch("trapEFCal", &trapEFCal);
  outtree->Branch("trapEFRCal", &trapEFRCal);
  outtree->Branch("avse", &avse);
  outtree->Branch("aoe", &aoe);
  outtree->Branch("dcr", &dcr);
  outtree->SetDirectory(outfile);

  // grab the calibration parameters if reading from an input configuration
  // otherwise compute them from the input data
  if(json_config){
    SetJson(jvalue, "base_mean",  base_mean);
    SetJson(jvalue, "base_rms",   base_rms);
    SetJson(jvalue, "pz_mean",    pz_mean);
    SetJson(jvalue, "escale",     escale);
    SetJson(jvalue, "efscale",    efscale);
    SetJson(jvalue, "eoffset",    eoffset);
    SetJson(jvalue, "efoffset",   efoffset);
    vector<double> p0, p1, p2, j;
    SetJson(jvalue, "avse_p0",    p0);
    SetJson(jvalue, "avse_p1",    p1);
    SetJson(jvalue, "avse_p2",    p1);
    SetJson(jvalue, "avse_j",     j);
    SetJson(jvalue, "rise_t",     trise_val);
    SetJson(jvalue, "rise_m",     trise_slope);
    SetJson(jvalue, "dcre_slope", dcre_slope);
    SetJson(jvalue, "dcr_cut_lo", dcr_cut_lo);
    SetJson(jvalue, "dcr_cut_hi", dcr_cut_hi);
  }
  else{
    cout << "getting initial energy scale calibration" << endl;
    for(auto const& p : chan_map){
      int i = p.second;
      // baseline, baseline rms, and pz constant
      int ebin0 = hbase_energy[i]->GetXaxis()->FindBin(200);
      int ebin1 = hbase_energy[i]->GetXaxis()->FindBin(5000);
      hbase[i] =  hbase_energy[i]->ProjectionY(("hbase_"+to_string(i)).c_str(),
					       ebin0, ebin1);
      hbrms[i] =  hbrms_energy[i]->ProjectionY(("hbrms_"+to_string(i)).c_str(),
					       ebin0, ebin1);
      base_mean[i] = hbase[i]->GetMean();
      brms_mean[i] = hbrms[i]->GetMean();
      base_rms[i]  = hbase[i]->GetRMS();
      brms_rms[i]  = hbrms[i]->GetRMS();
      int tbin0 = hdecay_deltat[i]->GetXaxis()->FindBin(50);
      int tbin1 = hdecay_deltat[i]->GetXaxis()->FindBin(250);
      hdecay[i] = hdecay_deltat[i]->ProjectionY(("hdecay_"+to_string(i)).c_str(),
						tbin0, tbin1);
      pz_mean[i] = hdecay[i]->GetMean();
      pz_rms[i]  = hdecay[i]->GetRMS();
      // initial energy scale correction
      if(source == "Th"){
	Th_scale(henergy[i], eoffset[i], escale[i]);
	Th_scale(henergyf[i], efoffset[i], efscale[i]);
      }
      else if(source == "Co"){
	Co_scale(henergy[i], eoffset[i], escale[i]);
	Co_scale(henergyf[i], efoffset[i], efscale[i]);
      }
      int nebins = henergy[i]->GetNbinsX();
      hecal[i] = new TH1D(("hecal_"+to_string(i)).c_str(), "", nebins,eoffset[i],
			  henergy[i]->GetXaxis()->GetBinUpEdge(nebins)*
			  escale[i] + eoffset[i]);
      for(int bin=1; bin<=nebins; bin++)
	hecal[i]->SetBinContent(bin, henergy[i]->GetBinContent(bin));
      hecal[i]->SetXTitle("Energy (keV)");
      hecal[i]->SetYTitle("Entries");
      // dcr slope correction
      TProfile* dcr_prof = hdcr_energy[i]->ProfileX(("hdcr_prof_" +
						     to_string(i)).c_str());
      TF1* fdcr = new TF1("fdcr", "[0]+[1]*x",0.0,dcr_prof->GetXaxis()->GetXmax());
      fdcr->FixParameter(0, 0.0);
      dcr_prof->Fit(fdcr, "QFR+");
      dcre_slope[i]   = fdcr->GetParameter(1);
      dcre_uncert[i] = fdcr->GetParError(1);
    }

    // calibrate AvsE and DCR
    cout << "calibrating avse and dcr" << endl;
    for(int ich=0; ich<(int)hamps_energy.size(); ich++){
      // 2d histograms to populate with the calibrated energy
      hamps_energy[ich] = new TH2D(("hamps_energy_"+to_string(ich)).c_str(),"",
				   hecal[ich]->GetXaxis()->GetNbins(),
				   hecal[ich]->GetXaxis()->GetXmin(),
				   hecal[ich]->GetXaxis()->GetXmax(),
				   hamp_energy[ich]->GetYaxis()->GetNbins(),
				   hamp_energy[ich]->GetYaxis()->GetXmin(),
				   hamp_energy[ich]->GetYaxis()->GetXmax());
      hamps_energy[ich]->SetXTitle("Energy (keV)");
      hamps_energy[ich]->SetYTitle("A * E/E_{unc}");
      hdcrs_energy[ich] = new TH2D(("hdcrs_energy_"+to_string(ich)).c_str(),"",
				   hecal[ich]->GetXaxis()->GetNbins(),
				   hecal[ich]->GetXaxis()->GetXmin(),
				   hecal[ich]->GetXaxis()->GetXmax(),
				   hdcr_energy[ich]->GetYaxis()->GetNbins(),
				   hdcr_energy[ich]->GetYaxis()->GetXmin(),
				   hdcr_energy[ich]->GetYaxis()->GetXmax());
      hdcrs_energy[ich]->SetXTitle("Energy (keV)");
      hdcrs_energy[ich]->SetYTitle("DCR Slope");
      haoe_energy[ich] = new TH2D(("haove_energy_"+to_string(ich)).c_str(), "",
				  hecal[ich]->GetXaxis()->GetNbins(),
				  hecal[ich]->GetXaxis()->GetXmin(),
				  hecal[ich]->GetXaxis()->GetXmax(),
				  hecal[ich]->GetXaxis()->GetNbins(),
				  hecal[ich]->GetXaxis()->GetXmin(),
				  hecal[ich]->GetXaxis()->GetXmax());
      haoe_energy[ich]->SetXTitle("Energy (keV)");
      haoe_energy[ich]->SetYTitle("A / E");
    }
    // populate the histograms
    while(reader.Next())
      for(int ich=0; ich<(int)channel->size(); ich++){
	if(use_fixedt){
	  double E = trappick->at(ich)*efscale[ich] + efoffset[ich];
	  hamps_energy[ich]->Fill(E, imax->at(ich)*E/trappick->at(ich));
	  haoe_energy[ich]->Fill(E, imax->at(ich)/trappick->at(ich));
	  hdcrs_energy[ich]->Fill(E, dcrslope->at(ich)-
				  dcre_slope[ich]*trappick->at(ich));
	}
	else{
	  double E = trapmax->at(ich)*escale[ich] + eoffset[ich];
	  hamps_energy[ich]->Fill(E, imax->at(ich)*E/trapmax->at(ich));
	  haoe_energy[ich]->Fill(E, imax->at(ich)/trapmax->at(ich));
	  hdcrs_energy[ich]->Fill(E, dcrslope->at(ich)-
				  dcre_slope[ich]*trapmax->at(ich)); 
	}
      }
    for(auto const& p : chan_map){
      int i = p.second;
      // get the AvsE profile by fitting in slices of energy
      int navsebins = (2550-50)/50;
      vector<double> avsex;
      vector<double> avsey;
      vector<double> avsexe;
      vector<double> avseye;
      for(int j=0; j<navsebins; j++){
	int bin0 = hamps_energy[i]->GetXaxis()->FindBin(25+j*50);
	int bin1 = hamps_energy[i]->GetXaxis()->FindBin(25+(j+1)*50);
	TH1D* h = hamps_energy[i]->ProjectionY("htmp", bin0, bin1);
	TF1* f = new TF1("ftmp", "gaus(0)+gaus(3)", 0.0, 1000.0);
	f->SetParameters(h->GetMaximum(),
			 h->GetBinCenter(h->GetMaximumBin()),
			 h->GetRMS()/2,
			 h->GetMaximum()*0.5,
			 h->GetBinCenter(h->GetMaximumBin())+h->GetStdDev(),
			 h->GetRMS()/2);
	h->Fit("ftmp", "QN");
	if(f->GetParameter(4)>0.0 && f->GetParameter(4)<10000 &&
	   f->GetParError(4)<10000){
	  avsex.push_back((1+j+0.5)*50);
	  avsexe.push_back(25);
	  if(f->GetParameter(0) > f->GetParameter(3)){
	    avsey.push_back(f->GetParameter(1));
	    avseye.push_back(f->GetParError(1));
	  }
	  else{
	    avsey.push_back(f->GetParameter(1));
	    avseye.push_back(f->GetParError(1));
	  }
	}
	delete h;
	delete f;
      }
      // fit to the AvsE profile with a quadratic
      TGraphErrors* g = new TGraphErrors(avsex.size(),
					 &avsex.front(),  &avsey.front(),
					 &avsexe.front(), &avseye.front());
      g->SetName(("gavse_"+to_string(i)).c_str());
      g->SetTitle("");
      g->GetHistogram()->SetXTitle("Energy (keV)");
      g->GetHistogram()->SetYTitle("Current Amplitude");
      g->SetMarkerColor(4);
      g->SetLineColor(4);
      TF1* f = new TF1(("favse_"+to_string(i)).c_str(),
		       "[0]+[1]*x+[2]*x*x",0.0, 20000.0);
      f->SetParameters(0.0, 1.0, 0.0);
      g->Fit(f, "QAF+", "", avsex[0], avsex[avsex.size()-1]);
      for(int j=0; j<3; j++){
	avse_param[i][j] = f->GetParameter(j);
	avse_uncert[i][j] = f->GetParError(j);
      }
      gavse[i] = g;
      delete f;
      // get AvsE distributions near the DEP, and do sideband subtraction
      TAxis* axis = hamps_energy[i]->GetXaxis();
      TH1D* havse0 = hamps_energy[i]->ProjectionY("havse0",
						  axis->FindBin(2615-2*511-10.),
						  axis->FindBin(2615-2*511+10.));
      TH1D* havse1 = hamps_energy[i]->ProjectionY("havse1",
						  axis->FindBin(2615-2*511-70.),
						  axis->FindBin(2615-2*511-50.));
      TH1D* havse2 = hamps_energy[i]->ProjectionY("havse2",
						  axis->FindBin(2615-2*511+70.),
						  axis->FindBin(2615-2*511+90.));
      for(int j=1; j<=(int)havse0->GetNbinsX(); j++){
	havse0->SetBinContent(j, havse0->GetBinContent(j) - 0.5 *
			      (havse1->GetBinContent(j)+havse2->GetBinContent(j)));
	if(j>1) havse0->SetBinContent(j, havse0->GetBinContent(j)+
				      havse0->GetBinContent(j-1));
      }
      // set the AvsE calibration parameters to 90% DEP acceptance
      double aint = havse0->GetMaximum();
      for(int j=1; j<=(int)havse0->GetNbinsX(); j++)
	if(havse0->GetBinContent(j)/aint > 0.1){
	  double x0 = havse0->GetBinCenter(j-1);
	  double y0 = havse0->GetBinContent(j-1)/aint;
	  double x1 = havse0->GetBinCenter(j);
	  double y1 = havse0->GetBinContent(j)/aint;
	  double E = (2615-2*511);
	  avse_param[i][3] = x0 + (x1-x0)*(0.1-y0)/(y1-y0);
	  avse_param[i][3] -= avse_param[i][0] + avse_param[i][1]*E;
	  avse_param[i][3] -= avse_param[i][2]*pow(E, 2);
	  break;
	}
      delete havse0;
      delete havse1;
      delete havse2;
      // find the upper and lower DCR cuts for 98% acceptance
      TH1D* hdcr_proj =
	hdcrs_energy[i]->ProjectionY("hdcr_proj",
				     hdcrs_energy[i]->GetXaxis()->FindBin(500),
				     hdcrs_energy[i]->GetXaxis()->FindBin(2500));
      double dintegral = hdcr_proj->Integral();
      hdcr_proj->SetBinContent(1, hdcr_proj->GetBinContent(1)/dintegral);
      for(int j=2; j<=(int)hdcr_proj->GetNbinsX(); j++){
	hdcr_proj->SetBinContent(j, hdcr_proj->GetBinContent(j)/dintegral +
				 hdcr_proj->GetBinContent(j-1));
	if(hdcr_proj->GetBinContent(j) > 0.01 && dcr_cut_lo[i] == 0.0)
	  dcr_cut_lo[i] = hdcr_proj->GetBinCenter(j);
	if(hdcr_proj->GetBinContent(j) > 0.99 && dcr_cut_hi[i] == 0.0)
	  dcr_cut_hi[i] = hdcr_proj->GetBinCenter(j);
      }
      delete hdcr_proj;
    }

    // compute the rise time correction (work in progress)
    cout << "computing rise time correction" << endl;
    reader.Restart();
    intree->SetBranchStatus("t1", true);
    intree->SetBranchStatus("t99", true);
    for(int ich=0; ich<(int)hrises_energy.size(); ich++){
      hrises_energy[ich] = new TH2D(("hrises_energy_"+to_string(ich)).c_str(), "",
				    hecal[ich]->GetXaxis()->GetNbins(),
				    hecal[ich]->GetXaxis()->GetXmin(),
				    hecal[ich]->GetXaxis()->GetXmax(),
				    hrise_energy[ich]->GetYaxis()->GetNbins(),
				    hrise_energy[ich]->GetYaxis()->GetXmin(),
				    1000);//hrise_energy[ich]->GetYaxis()->GetXmax());
      hrises_energy[ich]->SetXTitle("Energy (keV)");
      hrises_energy[ich]->SetYTitle("Rise Time (ns)");
    }
    while(reader.Next())
      for(int ich=0; ich<(int)channel->size(); ich++){
	double E = trappick->at(ich)*efscale[ich] + efoffset[ich];
	if(!use_fixedt) E = trapmax->at(ich)*escale[ich] + eoffset[ich];
	double avse = avse_param[ich][0] + avse_param[ich][1]*E;
	avse += avse_param[ich][2]*pow(E,2) - imax->at(ich);
	if(avse/avse_param[ich][3] <= -1.0) continue;
	hrises_energy[ich]->Fill(E, t99->at(ich)-t1->at(ich));
      }
    for(int ich=0; ich<(int)hrises_energy.size(); ich++){
      vector<double> rise_slope;
      vector<double> rise_error;
      vector<double> euncert(tl_peaks.size(), 1.0);
      TProfile* hproj =
	hrises_energy[ich]->ProfileX(("hrises_energy_proj_"+to_string(ich)).c_str(),
				     hrises_energy[ich]->GetYaxis()->FindBin(30.0),
				     hrises_energy[ich]->GetYaxis()->FindBin(60.0));
      TF1* fflat = new TF1(("fflat_"+to_string(ich)).c_str(), "[0]", 300, 3000);
      fflat->SetLineColor(8);
      hproj->Fit(fflat, "QFR+");
      trise_val[ich] = fflat->GetParameter(0);
      trise_uncert[ich] = fflat->GetParError(0);
      for(int ip=0; ip<(int)tl_peaks.size(); ip++){
	TF1* fn = new TF1(("frise_"+to_string(ich)+"_"+to_string(ip)).c_str(),
			  "[0]+[1]*x", tl_peaks[ip]*0.99, tl_peaks[ip]*1.01);
	hproj->Fit(fn, "QFR+");
	rise_slope.push_back(fn->GetParameter(1));
	rise_error.push_back(fn->GetParError(1));
      }
      hrises_energy_proj[ich] = hproj;
      hrise_slope[ich] = new TGraphErrors(tl_peaks.size(),
					  &tl_peaks.front(), &rise_slope.front(),
					  &euncert.front(), &rise_error.front());
      hrise_slope[ich]->SetName(("grise_slope_"+to_string(ich)).c_str());
      hrise_slope[ich]->SetTitle("");
      hrise_slope[ich]->GetHistogram()->SetXTitle("Energy (keV)");
      hrise_slope[ich]->GetHistogram()->SetYTitle("Rise Time Slope (ns/keV)");
      TF1* fconst = new TF1(("frise_const_"+to_string(ich)).c_str(),
			    "[0]", tl_peaks[0]-5, tl_peaks[tl_peaks.size()-1]+5);
      hrise_slope[ich]->Fit(fconst, "QFR+");
      trise_slope[ich] = fconst->GetParameter(0);
      trise_suncert[ich] = fconst->GetParError(0);
    }

    // print the calibration parameters
    for(auto const& p : chan_map){
      int i = p.second;
      cout << "channel 0:  " << endl;
      print_value("base    ", 3, base_mean[i],     base_rms[i]);
      print_value("brms    ", 3, brms_mean[i],     brms_rms[i]);
      print_value("pz      ", 3, pz_mean[i],       pz_rms[i]);
      print_value("avse p0 ", 3, avse_param[i][0], avse_uncert[i][0]);
      print_value("avse p1 ", 3, avse_param[i][1], avse_uncert[i][1]);
      print_value("avse p2 ", 3, avse_param[i][2], avse_uncert[i][2]);
      print_value("avse j  ", 3, avse_param[i][3], avse_uncert[i][3]);
      print_value("rise t  ", 3, trise_val[i],     trise_uncert[i]);
      print_value("rise m  ", 3, trise_slope[i],   trise_suncert[i]);
      print_value("dcr  m  ", 3, dcre_slope[i],    dcre_uncert[i]);
      print_value("dcr  lo ", 3, dcr_cut_lo[i],    0.0);
      print_value("dcr  hi ", 3, dcr_cut_hi[i],    0.0);
    }
  }

  intree->SetBranchStatus("*", true);
  reader.Restart();
  TTreeReaderValue<vector<int> > maxtime(reader, "maxtime");
  TTreeReaderValue<vector<int> > mintime(reader, "mintime");
  TTreeReaderValue<vector<int> > trapmaxtime(reader, "trapmaxtime");
  TTreeReaderValue<vector<double> > baseline(reader, "baseline");
  TTreeReaderValue<vector<double> > baserms(reader, "baserms");
  TTreeReaderValue<vector<double> > maxval(reader, "maxval");
  TTreeReaderValue<vector<double> > minval(reader, "minval");
  TTreeReaderValue<vector<double> > t0(reader, "t0");
  TTreeReaderValue<vector<double> > t10(reader, "t10");
  TTreeReaderValue<vector<double> > t50(reader, "t50");
  TTreeReaderValue<vector<double> > t90(reader, "t90");
  TTreeReaderValue<vector<double> > times(reader, "time");
  TTreeReaderValue<vector<double> > deltat(reader, "deltat");
  TTreeReaderValue<vector<vector<double> > > exp_param(reader, "exp_param");

  // apply the calibration parameters to the relevant input values
  cout << "applying calibration parameters to input trees" << endl;
  int iev = -1;
  while(reader.Next()){
    iev ++;
    if(iev % 100000 == 0 && iev > 0) outtree->AutoSave();
    int nwf = (int) channel->size();
    chan.assign(nwf, 0);
    baseSigma.assign(nwf, 0.0);
    brmsSigma.assign(nwf, 0.0);
    timestamp.assign(nwf, 0.0);
    dt.assign(nwf, 0.0);
    T0.assign(nwf, 0.0);
    trise.assign(nwf, 0.0);
    trapECal.assign(nwf, 0.0);
    trapEFCal.assign(nwf, 0.0);
    trapEFRCal.assign(nwf, 0.0);
    avse.assign(nwf, 0.0);
    aoe.assign(nwf, 0.0);
    dcr.assign(nwf, 0.0);
    for(int ich=0; ich<nwf; ich++){
      chan[ich] = channel->at(ich);
      baseSigma[ich] = (baseline->at(ich)-base_mean[ich])/base_rms[ich];
      brmsSigma[ich] = (baserms->at(ich)-brms_mean[ich])/brms_rms[ich];
      timestamp[ich] = times->at(ich);
      T0[ich] = t0->at(ich);
      dt[ich] = deltat->at(ich);
      trise[ich] = t99->at(ich)-t1->at(ich);
      trapECal[ich] = trapmax->at(ich)*escale[ich]+eoffset[ich];
      if(use_fixedt)
	trapEFCal[ich] = trappick->at(ich)*efscale[ich]+efoffset[ich];
      else
	trapEFCal[ich] = trapmax->at(ich)*escale[ich]+eoffset[ich];
      trapEFRCal[ich] = trapEFCal[ich] + ((trise_val[ich] -
					   trise[ich])) / trise_slope[ich];
      avse[ich] = avse_param[ich][0] + avse_param[ich][1]*trapEFCal[ich];
      avse[ich] += avse_param[ich][2]*pow(trapEFCal[ich], 2);
      avse[ich] -= imax->at(ich)*trapEFCal[ich]/trapEFCal[ich];
      avse[ich] /= avse_param[ich][3];
      aoe[ich] = imax->at(ich) / trapEFCal[ich];
      dcr[ich] = dcrslope->at(ich)-dcre_slope[ich]*trappick->at(ich);
      dcr[ich] = dcr[ich]*2/(dcr_cut_hi[ich]-dcr_cut_lo[ich]);
      dcr[ich] += (dcr_cut_lo[ich]+dcr_cut_hi[ich])/
	(dcr_cut_lo[ich]-dcr_cut_hi[ich]);
    }
    outtree->Fill();
  }

  cout << "writing output file" << endl;
  outfile->cd();
  outtree->AutoSave();
  for(auto const& p : chan_map){
    int i = p.second;
    if(hdeltat[i])
      hdeltat[i]->Write(("hdeltat_"+to_string(p.first)).c_str());
    if(henergy[i])
      henergy[i]->Write(("henergy_"+to_string(p.first)).c_str());
    if(henergyf[i])
      henergyf[i]->Write(("henergyf_"+to_string(p.first)).c_str());
    if(hbase[i])
      hbase[i]->Write(("hbase_"+to_string(p.first)).c_str());
    if(hbrms[i])
      hbrms[i]->Write(("hbrms_"+to_string(p.first)).c_str());
    if(hdecay[i])
      hdecay[i]->Write(("hdecay_"+to_string(p.first)).c_str());
    if(hbase_energy[i])
      hbase_energy[i]->Write(("hbase_energy_"+to_string(p.first)).c_str());
    if(hbrms_energy[i])
      hbrms_energy[i]->Write(("hbrms_energy_"+to_string(p.first)).c_str());
    if(hdecay_energy[i])
      hdecay_energy[i]->Write(("hdecay_energy_"+to_string(p.first)).c_str());
    if(hdecay_deltat[i])
      hdecay_deltat[i]->Write(("hdecay_deltat_"+to_string(p.first)).c_str());
    if(hamp_energy[i])
      hamp_energy[i]->Write(("hamp_energy_"+to_string(p.first)).c_str());
    if(hamps_energy[i])
      hamps_energy[i]->Write(("hamps_energy"+to_string(p.first)).c_str());
    if(hdcr_energy[i])
      hdcr_energy[i]->Write(("hdcr_energy_"+to_string(p.first)).c_str());
    if(hdcrs_energy[i]) hdcrs_energy[i]->Write();
    if(hrise_energy[i])
      hrise_energy[i]->Write(("hrise_energy_"+to_string(p.first)).c_str());
    if(hrises_energy[i]) hrises_energy[i]->Write();
    if(hrise_slope[i]) hrise_slope[i]->Write();
    if(hrises_energy_proj[i]) hrises_energy_proj[i]->Write();
    if(henergy[i]) henergy[i]->Write();
    if(hecal[i]) hecal[i]->Write();
    if(gavse[i]) gavse[i]->Write();
  }
  outfile->WriteObject(&chan_map,  "chan_map");
  outfile->WriteObject(&base_mean, "base_mean");
  outfile->WriteObject(&base_rms,  "base_rms");
  outfile->WriteObject(&brms_mean, "brms_mean");
  outfile->WriteObject(&brms_rms,  "brms_rms");
  outfile->WriteObject(&pz_mean,   "pz_mean");
  outfile->WriteObject(&pz_rms,    "pz_rms");
  outfile->Close();

  // write calibration parameters to a json file if not reading from one
  if(!json_config){
    cout << "writing calibration parameters to " << jfile << endl;
    for(int i=0; i<(int) base_mean.size(); i++){
      jvalue["base_mean"][i]  = base_mean[i];
      jvalue["base_rms"][i]   = base_rms[i];
      jvalue["pz_mean"][i]    = pz_mean[i];
      jvalue["escale"][i]     = escale[i];
      jvalue["eoffset"][i]    = eoffset[i];
      jvalue["efscale"][i]    = efscale[i];
      jvalue["efoffset"][i]   = efoffset[i];
      jvalue["avse_p0"][i]    = avse_param[i][0];
      jvalue["avse_p1"][i]    = avse_param[i][1];
      jvalue["avse_p2"][i]    = avse_param[i][2];
      jvalue["avse_j"][i]     = avse_param[i][3];
      jvalue["rise_t"][i]     = trise_val[i];
      jvalue["rise_m"][i]     = trise_slope[i];
      jvalue["dcre_slope"][i] = dcre_slope[i];
      jvalue["dcr_cut_lo"][i] = dcr_cut_lo[i];
      jvalue["dcr_cut_hi"][i] = dcr_cut_hi[i];
    }
    Json::StreamWriterBuilder builder;
    builder["commentStyle"] = "None";
    builder["indentation"]  = "    ";
    unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
    ofstream fstream(jfile);
    writer->write(jvalue, &fstream);
    fstream.close();
  }
      
  return 0;
}
