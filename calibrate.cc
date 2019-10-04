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
#include <algorithm>
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

// do a two point energy calibration to get the initial scale
void Th_scale(TH1D* h, double& offset, double& scale, int peak=0){
  // find the bin below which most of the spectrum lies
  int bin = 0;
  if(peak > 0) bin = h->GetXaxis()->FindBin(peak);
  else{
    for(int i=h->GetNbinsX(); i>200; i--){
      double c1 = h->Integral(i-199, i-100);
      double c0 = h->Integral(i-99,  i);
      if(c0 <= 10.0 || c1 <= 10.0) continue;
      if(c1-c0 > 20*sqrt(c0)){
	bool valid = true;
	for(int j=0; j<min(500, i-200); j++){
	  c1 = h->Integral(i-200-j, i-101-j);
	  if(c1 <= 10.0 || c1-c0 <= 10*sqrt(c0)){
	    valid = false;
	    break;
	  }
	}
	if(!valid) continue;
	bin = i+100;
	break;
      }
    }
  }
  // find the maximum in the upper third of that region, assume 2615 keV
  int bin1 = 0;
  double mval = 0.0;
  if(peak > 0){
    int b0 = h->GetXaxis()->FindBin(0.9*peak);
    int b1 = h->GetXaxis()->FindBin(1.1*peak);
    for(int i=b0; i<=b1; i++)
      if(h->GetBinContent(i) > mval){
	mval = h->GetBinContent(i);
	bin1 = i;
      }
  }
  else
    for(int i=max(1, (int)(2.*bin/3)); i<bin; i++)
      if(h->GetBinContent(i) > mval){
	mval = h->GetBinContent(i);
	bin1 = i;
      }
  double bk = bin / 2615.;
  TF1* f1 = new TF1("f1", "[0]+gaus(1)",
		    h->GetXaxis()->GetBinLowEdge(bin1-max(2, (int)(bk*4))),
		    h->GetXaxis()->GetBinUpEdge( bin1+max(2, (int)(bk*4))));
  double base = h->Integral((int)(bin1-65*bk), (int)(bin1-55*bk));
  base += h->Integral((int)(bin1+35*bk), (int)(bin1+45*bk));
  base /= 20;
  f1->SetParameters(base, h->GetBinContent(bin1),
		    h->GetBinCenter(bin1), 2*bk);
  f1->FixParameter(0, base);
  f1->SetParLimits(2, h->GetBinCenter(bin1)-bk/4, h->GetBinCenter(bin1)+bk/4);
  h->Fit(f1, "QMR+");
  double val1 = f1->GetParameter(2);
  // find the maximum at the approximate position of 583 keV
  int bin0 = (int) (bin1*583./2615);
  int mbin = 0;
  mval = 0.0;
  for(int i=(int)(bin0-0.005*bin1); i<(int)(bin0+0.005*bin1); i++)
    if(h->GetBinContent(i) > mval){
      mval = h->GetBinContent(i);
      mbin = i;
    }
  bin0 = mbin;
  TF1* f0 = new TF1("f0", "[0]+gaus(1)",
		    h->GetXaxis()->GetBinLowEdge(bin0-max(2, (int)(bk*2))),
		    h->GetXaxis()->GetBinUpEdge( bin0+max(2, (int)(bk*2))));
  base =  h->Integral((int)(bin0-30*bk), (int)(bin0-20*bk));
  base += h->Integral((int)(bin0+7*bk), (int)(bin0+17*bk));
  base /= (int) (22*bk);
  f0->SetParameters(base, h->GetBinContent(bin0)-base,
		    h->GetBinCenter(bin0), bk);
  f0->FixParameter(0, base);
  f0->SetParLimits(2, h->GetBinCenter(bin0)-bk/4, h->GetBinCenter(bin0)+bk/4);
  h->Fit(f0, "QRM+");
  double val0 = f0->GetParameter(2);
  scale = (2615-583) / (val1-val0);
  offset = 2615 - scale*val1;
}

// get a rough estimate of the resolution of the 2615 keV peak
pair<double, double> GetThRes(TH1D* h, pair<double, double> pr,
			      pair<double, double> sb0,
			      pair<double, double> sb1, double pos=0.0){
  double e0 = 0.0;
  double e1 = 0.0;
  if(pos == 0.0) Th_scale(h, e0, e1);
  else Th_scale(h, e0, e1, pos);
  int b0 = (int) h->GetXaxis()->FindBin((sb0.first-e0)/e1);
  int b1 = (int) h->GetXaxis()->FindBin((sb0.second-e0)/e1);
  int b2 = (int) h->GetXaxis()->FindBin((pr.first-e0)/e1);
  int b3 = (int) h->GetXaxis()->FindBin((pr.second-e0)/e1);
  int b4 = (int) h->GetXaxis()->FindBin((sb1.first-e0)/e1);
  int b5 = (int) h->GetXaxis()->FindBin((sb1.second-e0)/e1);
  TF1* f = new TF1("ftmp2", "pol0(0)+gaus(1)",
		   h->GetBinCenter(b1), h->GetBinCenter(b4));
  f->SetLineColor(8);
  f->FixParameter(0, (h->Integral(b0, b1)/(1+b1-b0)+
		      h->Integral(b4, b5)/(1+b5-b4))/2);
  f->SetParameter(1, h->GetBinContent((int)((b3+b2)*0.5))-
		  f->GetParameter(0));
  f->SetParameter(2, ((pr.first+pr.second)/2-e0)/e1);
  f->SetParameter(3, 3/e1);
  f->SetParLimits(2, h->GetBinCenter(b2), h->GetBinCenter(b3));
  h->Fit(f, "QMR+");
  pair<double, double> res;
  res.first = f->GetParameter(3);
  res.second = f->GetParError(3);
  delete f;
  return res;
}

// calibrate the AvsE cut from a scaled A vs E histogram
TGraphErrors* AvsECal(TH2D* ha, vector<double>& param, vector<double>& uncert){
  // values below from Walter, AvsE paper
  const vector<double> start({200,310,400,530,600,700,800,900,1050,1180,
	1250,1300,1465,1550,1680,1750,1800,1900,2000,2150,2220,2300});
  const double w = 25.0;
  const double dep = 2615 - 2*511;
  // get the AvsE profile by fitting in slices of energy
  vector<double> x;
  vector<double> y;
  vector<double> ex;
  vector<double> ey;
  for(int i=0; i<(int)start.size(); i++){
    int bin0 = ha->GetXaxis()->FindBin(start[i]);
    int bin1 = ha->GetXaxis()->FindBin(start[i]+w);
    TH1D* h = ha->ProjectionY("htmp", bin0, bin1);
    if(h->Integral() == 0.0){
      delete h;
      continue;
    }
    x.push_back(start[i]+0.5*w);
    ex.push_back(w/5);
    y.push_back(h->GetXaxis()->GetBinCenter(h->GetMaximumBin()));
    ey.push_back(h->GetMeanError());
    delete h;
  }
  if(x.size() == 0){
    cout << "AvsE calibration failed due to empty histogram" << endl;
    return NULL;
  }
  // fit to the AvsE profile with a quadratic
  TGraphErrors* g = new TGraphErrors(x.size(), &x.front(),  &y.front(),
				     &ex.front(), &ey.front());
  g->SetName(("gavse_"+string(ha->GetName())).c_str());
  g->SetTitle("");
  g->GetHistogram()->SetXTitle("Energy (keV)");
  g->GetHistogram()->SetYTitle("Current Amplitude");
  g->SetMarkerColor(4);
  g->SetLineColor(4);
  TF1* f = new TF1(("favse_"+string(ha->GetName())).c_str(),
		   "[0]+[1]*x+[2]*x*x",0.0, 20000.0);
  f->SetParameters(0.0, 1.0, 0.0);
  f->SetParLimits(2, -1e3, 0.0);
  g->Fit(f, "QAF+", "", x.front(), x.back());
  for(int i=0; i<3; i++){
    param[i] = f->GetParameter(i);
    uncert[i] = f->GetParError(i);
  }
  delete f;
  // get AvsE distributions near the DEP, and do sideband subtraction
  TAxis* axis = ha->GetXaxis();
  TH1D* havse0 = ha->ProjectionY("havse0",
				 axis->FindBin(dep-5),
				 axis->FindBin(dep+5));
  TH1D* havse1 = ha->ProjectionY("havse1",
				 axis->FindBin(dep-50),
				 axis->FindBin(dep-20));
  TH1D* havse2 = ha->ProjectionY("havse2",
				 axis->FindBin(dep+50),
				 axis->FindBin(dep+80));
  for(int i=0; i<=(int)havse0->GetNbinsX(); i++){
    double val = havse0->GetBinContent(i);
    val -= (havse1->GetBinContent(i)+havse2->GetBinContent(i)) * 0.5;
    havse0->SetBinContent(i, val);
    if(i>1) havse0->SetBinContent(i, havse0->GetBinContent(i)+
				  havse0->GetBinContent(i-1));
  }
  // set the AvsE calibration parameters to 90% DEP acceptance
  double aint = havse0->GetMaximum();
  for(int i=2; i<=(int)havse0->GetNbinsX(); i++)
    if(havse0->GetBinContent(i)/aint > 0.1){
      double x0 = havse0->GetBinCenter(i-1);
      double y0 = havse0->GetBinContent(i-1)/aint;
      double x1 = havse0->GetBinCenter(i);
      double y1 = havse0->GetBinContent(i)/aint;
      param[3] = x0 + (x1-x0)*(0.1-y0)/(y1-y0);
      param[3] -= param[0] + param[1]*dep + param[2]*pow(dep, 2);
      break;
    }
  delete havse0;
  delete havse1;
  delete havse2;
  return g;
}

// calibrate A/E from the scaled A/E vs energy histogram
void AoECal(TH2D* h, double& amin, double& auncert){
  const double dep = 2615. - 2*511;
  // get the A/E distribution of the DEP and sideband subtract
  TAxis* axis = h->GetXaxis();
  TH1D* havse0 = h->ProjectionY("havse0",
				 axis->FindBin(dep-5),
				 axis->FindBin(dep+5));
  TH1D* havse1 = h->ProjectionY("havse1",
				 axis->FindBin(dep-25),
				 axis->FindBin(dep-15));
  TH1D* havse2 = h->ProjectionY("havse2",
				 axis->FindBin(dep+15),
				 axis->FindBin(dep+25));
  for(int i=1; i<=(int)havse0->GetNbinsX(); i++){
    double val = havse0->GetBinContent(i);
    val -= (havse1->GetBinContent(i)+havse2->GetBinContent(i)) * 0.5;
    havse0->SetBinContent(i, val);
    if(i>1) havse0->SetBinContent(i, havse0->GetBinContent(i)+
				  havse0->GetBinContent(i-1));
  }
  double aint = havse0->GetMaximum();
  for(int i=2; i<=(int)havse0->GetNbinsX(); i++)
    if(havse0->GetBinContent(i)/aint > 0.1){
      double x0 = havse0->GetBinCenter(i-1);
      double y0 = havse0->GetBinContent(i-1)/aint;
      double x1 = havse0->GetBinCenter(i);
      double y1 = havse0->GetBinContent(i)/aint;
      amin = x0 + (x1-x0)*(0.1-y0)/(y1-y0);
      auncert = havse0->GetBinWidth(i)*0.5;
      break;
    }
}

// calibrate the DCR cut from the corrected DCR vs energy
void DCRCal(TH2D* h, double rejlo, double rejhi, double& cutlo, double& cuthi){
  TH1D* hp = h->ProjectionY("hdcr_proj",
			    h->GetXaxis()->FindBin(1850),
			    h->GetXaxis()->FindBin(2050));
  TH1D* hp2 = h->ProjectionY("hdcr_proj2",
			     h->GetXaxis()->FindBin(2150),
			     h->GetXaxis()->FindBin(2350));
  hp->Add(hp2);
  double dint = hp->Integral();
  hp->SetBinContent(1, hp->GetBinContent(1)/dint);
  for(int i=2; i<=(int)hp->GetNbinsX(); i++){
    hp->SetBinContent(i, hp->GetBinContent(i)/dint + hp->GetBinContent(i-1));
    if(hp->GetBinContent(i) > rejlo && cutlo == 0.0)
      cutlo = hp->GetBinCenter(i);
    if(hp->GetBinContent(i) > (1-rejhi) && cuthi == 0.0)
      cuthi = hp->GetBinCenter(i);
  }
  delete hp;
  delete hp2;
}


int main(int argc, char* argv[]){

  // root setup
  gStyle->SetOptStat(0);
  gStyle->SetLabelFont(132, "XYZ");
  gStyle->SetTitleFont(132, "XYZ");
  tdir = gROOT->CurrentDirectory();

  // Tl peaks to use for the rise time correction
  //const vector<double> tl_peaks({338.0, 510.5, 583.0, 726.5,
  //	794.0, 860.0, 911.0, 968.0, 2615.0});

  //option handling
  map<int, int> chan_map;
  vector<int> chan_cal;
  vector<int> chan_skip;
  vector<string> ifname;
  string ofname = "";
  bool json_config = false;
  string jfile = "calibrate.json";
  Json::Value jvalue;
  bool use_fixedt = true;
  map<int, int> thpeak;
  bool get_ct_decay = false;
  static struct option opts[]{
    {"help",             no_argument, NULL, 'h'},
    {"channel",    required_argument, NULL, 'c'},
    {"infile",     required_argument, NULL, 'i'},
    {"outfile",    required_argument, NULL, 'o'},
    {"skipchan",   required_argument, NULL, 's'},
    {"jsonconfig", required_argument, NULL, 'j'},
    {"jsonfile",   required_argument, NULL, 'J'},
    {"fixedtime",        no_argument, NULL, 'f'},
    {"thlocation", required_argument, NULL, 't'},
    {"chargetrap", no_argument, NULL, 'C'}
  };
  int opt = getopt_long(argc, argv, "hc:i:o:s:j:J:ft:C", opts, NULL);
  while(opt != -1){
    switch(opt){
    case 'h':
      cout << "options:"                               << endl;
      cout << "  -c channel to be calibrated"          << endl;
      cout << "  -i input filename"                    << endl;
      cout << "  -o output filename"                   << endl;
      cout << "  -s channel to be skipped"             << endl;
      cout << "  -j name of input json config file"    << endl;
      cout << "  -J name of output json config file"   << endl;
      cout << "  -f do not use the fixed time pickoff" << endl;
      cout << "  -C use charge trapping correction"    << endl;
      cout << "  -t channel,adc where adc is the estimated "
	   << " 2615 keV peak loaction"                << endl;
      return 0;
    case 'c': chan_cal.push_back(atoi(optarg));  break;
    case 'i': ifname.push_back(string(optarg));  break;
    case 'o': ofname = string(optarg);           break;
    case 's': chan_skip.push_back(atoi(optarg)); break;
    case 'j':{
      ifstream jfile(optarg);
      Json::CharReaderBuilder reader;
      string errors;
      if(!parseFromStream(reader, jfile, &jvalue, &errors)){
	cout << errors << endl;
	return 2;
      }
      json_config = true;
      jfile.close();
      break;
    }
    case 'J': jfile = string(optarg); break;
    case 'f': use_fixedt = false;     break;
    case 't':{ 
      string a = string(optarg);
      int chan = stoi(a.substr(0, a.find(",")));
      thpeak[chan] = stoi(a.substr(a.find(",")+1, a.size()-1));
      break;
    }
    case 'C': get_ct_decay = true; break;
    default: return 1;
    }
    opt = getopt_long(argc, argv, "hc:i:o:s:j:J:ft:C", opts, NULL);
  }
  assert(ofname != "" && ifname.size() > 0);
  
  // read in input files, setup input tree
  TChain* intree = new TChain("tree");
  vector<TH1D*> hdeltat, henergy, henergyf, henergyc1, henergyc2;
  vector<TH1D*> hbase, hbrms, hdecay;
  vector<TH2D*> hbase_energy, hbrms_energy, hdecay_energy, hdecay_deltat;
  vector<TH2D*> hamp_energy, hamps_energy, haoe_energy, haoes_energy;
  vector<TH2D*> hdcr_energy, hdcrs_energy, hrise_energy, hrises_energy;
  cout << "reading input histograms and configuration" << endl;
  string process_config = "";
  vector<string> serial_numbers;
  vector<double> pz_decay;
  vector<double> ct_decay;
  vector<double> ct_frac;
  vector<int> ct_method;
  for(auto const& n : ifname){
    intree->Add(n.c_str());
    TFile* f = TFile::Open(n.c_str());
    tdir->cd();
    string config(((TObjString*) f->Get("process_config"))->GetName());
    if(process_config == ""){
      process_config = config;
      // fixme - comments in the json will break this
      Json::Value value;
      Json::CharReaderBuilder builder;
      Json::CharReader* reader = builder.newCharReader();
      string jerrors;
      if(!reader->parse(config.c_str(),
			config.c_str()+config.size(), &value, &jerrors)){
	cout << jerrors << endl;
	return 2;
      }
      delete reader;
      vector<double> chans;
      SetJson(value, "channel_map", chans);
      for(auto const& i : chans){
	unsigned s = chan_map.size();
	chan_map[i] = s;
      }
      serial_numbers.resize(chan_map.size());
      pz_decay.resize(chan_map.size());
      ct_decay.resize(chan_map.size());
      ct_frac.resize(chan_map.size());
      ct_method.resize(chan_map.size());
      vector<int> channels;
      vector<string> detectors;
      SetJson(value, "channel_id", channels);
      SetJson(value, "det_serial", detectors);
      for(auto const& pr : chan_map){
	vector<int>::iterator i = find(channels.begin(),
				       channels.end(), pr.first);
	int index = distance(channels.begin(), i);
	Json::Value val = value[detectors[index]];
	serial_numbers[pr.second] = detectors[index];
	SetJson(val, "pz_decay", pz_decay[pr.second]);
	SetJson(val, "ct_decay", ct_decay[pr.second]);
	SetJson(val, "ct_frac",  ct_frac[pr.second]);
	SetJson(val, "ct_method", ct_method[pr.second]);
      }
      hdeltat.resize(chan_map.size(), NULL);
      henergy.resize(chan_map.size(), NULL);
      henergyf.resize(chan_map.size(), NULL);
      henergyc1.resize(chan_map.size(), NULL);
      henergyc2.resize(chan_map.size(), NULL);
      hbase.resize(chan_map.size(), NULL);
      hbrms.resize(chan_map.size(), NULL);
      hdecay.resize(chan_map.size(), NULL);
      hbase_energy.resize(chan_map.size(), NULL);
      hbrms_energy.resize(chan_map.size(), NULL);
      hdecay_energy.resize(chan_map.size(), NULL);
      hdecay_deltat.resize(chan_map.size(), NULL);
      hamp_energy.resize(chan_map.size(), NULL);
      hamps_energy.resize(chan_map.size(), NULL);
      haoe_energy.resize(chan_map.size(), NULL);
      haoes_energy.resize(chan_map.size(), NULL);
      hdcr_energy.resize(chan_map.size(), NULL);
      hdcrs_energy.resize(chan_map.size(), NULL);
      hrise_energy.resize(chan_map.size(), NULL);
      hrises_energy.resize(chan_map.size(), NULL);
    }
    else if(config != process_config){
      cout << "input file " << n << " processed with configuration:" << endl
	   << config << endl << "current configuration is:"
	   << process_config << endl;
      return 3;
    } 
    for(auto const& p : chan_map){
      string s = "_" + to_string(p.first);
      int i = p.second;
      hdeltat[i]      = Get1DHistOrSum(f, "hdeltat"      +s, hdeltat[i]);
      henergy[i]      = Get1DHistOrSum(f, "henergy"      +s, henergy[i]);
      henergyf[i]     = Get1DHistOrSum(f, "henergyf"     +s, henergyf[i]);
      henergyc1[i]    = Get1DHistOrSum(f, "henergyc1"    +s, henergyc1[i]);
      henergyc2[i]    = Get1DHistOrSum(f, "henergyc2"    +s, henergyc2[i]);
      hbase_energy[i] = Get2DHistOrSum(f, "hbase_energy" +s, hbase_energy[i]);
      hbrms_energy[i] = Get2DHistOrSum(f, "hbrms_energy" +s, hbrms_energy[i]);
      hdecay_energy[i]= Get2DHistOrSum(f, "hdecay_energy"+s, hdecay_energy[i]);
      hdecay_deltat[i]= Get2DHistOrSum(f, "hdecay_deltat"+s, hdecay_deltat[i]);
      if(get_ct_decay)      s = "c" + s;
      else if(use_fixedt)   s = "f" + s;
      hamp_energy[i]  = Get2DHistOrSum(f, "hamp_energy"  +s, hamp_energy[i]);
      haoe_energy[i]  = Get2DHistOrSum(f, "haoe_energy"  +s, haoe_energy[i]);
      hdcr_energy[i]  = Get2DHistOrSum(f, "hdcr_energy"  +s, hdcr_energy[i]);
      hrise_energy[i] = Get2DHistOrSum(f, "hrise_energy" +s, hrise_energy[i]);
    }
    f->Close();
  }
  intree->SetBranchStatus("*", false);
  intree->SetBranchStatus("channel", true);
  intree->SetBranchStatus("ct_decay", true);
  intree->SetBranchStatus("ct_value", true);
  intree->SetBranchStatus("ct_integral", true);
  intree->SetBranchStatus("trappick", true);
  TTreeReader reader(intree);
  TTreeReaderValue<vector<int> > channel(reader, "channel");
  TTreeReaderValue<vector<double> > trappick(reader, "trappick");
  TTreeReaderValue<vector<double> > ct1_trappick(reader, "ct1_trappick");
  TTreeReaderValue<vector<double> > ct2_trappick(reader, "ct2_trappick");
  TTreeReaderValue<vector<double> > ct_integral(reader, "ct_integral");
  TTreeReaderValue<vector<double> > trapmax(reader, "trapmax");
  TTreeReaderValue<vector<double> > imax(reader, "imax");
  TTreeReaderValue<vector<double> > dcrslope(reader, "dcrslope");
  TTreeReaderValue<vector<double> > t1(reader, "t1");
  TTreeReaderValue<vector<double> > t99(reader, "t99");
  TTreeReaderValue<vector<double> > sampling(reader, "sampling");
  TTreeReaderValue<vector<vector<double> > > ct_dec(reader, "ct_decay");
  TTreeReaderValue<vector<vector<double> > > ct_val(reader, "ct_value");
  
  // values and histograms to output
  vector<double> base_mean(chan_map.size(), 0.0);
  vector<double> brms_mean(chan_map.size(), 0.0);
  vector<double> base_uncert(chan_map.size(), 0.0);
  vector<double> brms_uncert(chan_map.size(), 0.0);
  vector<double> pz_mean(chan_map.size(), 0.0);
  vector<double> pz_uncert(chan_map.size(), 0.0);
  vector<double> ct1_mean(chan_map.size(), 0.0);
  vector<double> ct2_mean(chan_map.size(), 0.0);
  vector<double> ct1_uncert(chan_map.size(), 0.0);
  vector<double> ct2_uncert(chan_map.size(), 0.0);
  vector<double> escale(chan_map.size(), 0.0);
  vector<double> efscale(chan_map.size(), 0.0);
  vector<double> ecscale(chan_map.size(), 0.0);
  vector<double> eoffset(chan_map.size(), 0.0);
  vector<double> efoffset(chan_map.size(), 0.0);
  vector<double> ecoffset(chan_map.size(), 0.0);
  vector<TH1D*> hecal(chan_map.size(), NULL);
  vector<vector<double> > avse_param(chan_map.size(), vector<double>(4, 0.0));
  vector<vector<double> > avse_uncert(chan_map.size(), vector<double>(4, 0.0));
  vector<TGraphErrors*> gavse(chan_map.size(), NULL);
  vector<double> aoe_min(chan_map.size(), 0.0);
  vector<double> aoe_min_uncert(chan_map.size(), 0.0);
  vector<TGraphErrors*> gresct1(chan_map.size(), NULL);
  vector<TGraphErrors*> gresct2(chan_map.size(), NULL);
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

  // determine the optimal charge trapping time constant
  const int nct_steps = 100;
  if(get_ct_decay){
    vector<vector<double> > ct(chan_map.size());
    vector<vector<double> > dfrac(chan_map.size());
    vector<vector<TH1D*> > hE1(chan_map.size());
    vector<vector<TH1D*> > hE2(chan_map.size());
    bool first = true;
    bool valid = true;
    while(reader.Next()){
      if(!valid) break;
      for(int ich=0; ich<(int)channel->size(); ich++){
	if(find(chan_cal.begin(),
		chan_cal.end(), channel->at(ich)) == chan_cal.end())
	  continue;
	int i = chan_map[channel->at(ich)];
	if(first){
	  first = false;
	  if(ct_val->size() == 0){
	    valid = false;
	    break;
	  }
	  for(int j=0; j<(int)ct.size(); j++){
	    ct[j].resize(ct_dec->at(ich).size(), 0.0);
	    hE1[j].resize(ct_dec->at(ich).size(), NULL);
	    dfrac[j].resize(nct_steps, 0.0);
	    hE2[j].resize(nct_steps, NULL);
	    for(int k=0; k<(int)ct[j].size(); k++){
	      ct[j][k] = ct_dec->at(ich)[k];
	      string n = "hct1_" + to_string(j) + "_" + to_string(k);
	      hE1[j][k] = new TH1D(n.c_str(), "",
				   henergyc1[i]->GetXaxis()->GetNbins(),
				   henergyc1[i]->GetXaxis()->GetXmin(),
				   henergyc1[i]->GetXaxis()->GetXmax());
	    }
	    for(int k=0; k<nct_steps; k++){
	      dfrac[j][k] = (1 + (99*k/((double)nct_steps))) * 1e-8;
	      string n = "hct2_" + to_string(j) + "_" + to_string(k);
	      hE2[j][k] = new TH1D(n.c_str(), "",
				   henergyc1[i]->GetXaxis()->GetNbins()*4,
				   henergyc1[i]->GetXaxis()->GetXmin(),
				   henergyc1[i]->GetXaxis()->GetXmax());
	    }
	  }
	}
	for(int j=0; j<(int)ct_val->at(ich).size(); j++)
	  hE1[i][j]->Fill(ct_val->at(ich)[j]);
	for(int j=0; j<nct_steps; j++)
	  hE2[i][j]->Fill(trappick->at(ich)+dfrac[i][j]*ct_integral->at(ich));
      }
    }
    if(valid){
      map<int, int>::iterator it = chan_map.begin();
      for(int i=0; i<(int)hE1.size(); i++){
	int ch = (*it).first;
	it ++;
	vector<double> ct_uncert(hE1[i].size(), 1000.0);
	vector<double> dfrac_uncert(hE2[i].size(), 0.1*dfrac[i][0]);
	vector<double> res1(hE1[i].size(), 0.0);
	vector<double> res1_uncert(hE1[i].size(), 0.0);
	vector<double> res2(hE2[i].size(), 0.0);
	vector<double> res2_uncert(hE2[i].size(), 0.0);
	for(int j=0; j<(int)hE1[i].size(); j++){
	  pair<double, double> pr;
	  if(thpeak.find(ch) == thpeak.end())
	    pr = GetThRes(hE1[i][j], make_pair(2610., 2620.),
			  make_pair(2570., 2580.), make_pair(2650., 2660.));
	  else
	    pr = GetThRes(hE1[i][j], make_pair(2610., 2620.),
			  make_pair(2570., 2580.), make_pair(2650., 2660.),
			  thpeak[ch]);
	  res1[j] = pr.first;
	  res1_uncert[j] = pr.second;
	}
	for(int j=0; j<(int)res1.size(); j++)
	  if(abs(res1_uncert[j] / res1[j]) > 0.2){
	    res1.erase(res1.begin()+j);
	    ct[i].erase(ct[i].begin()+j);
	    res1_uncert.erase(res1_uncert.begin()+j);
	    ct_uncert.erase(ct_uncert.begin()+j);
	  }
	TGraphErrors* g1 = new TGraphErrors(res1.size(), &ct[i].front(),
					    &res1.front(), &ct_uncert.front(),
					    &res1_uncert.front());
	g1->SetName(("gresct1_"+to_string(i)).c_str());
	g1->SetTitle("");
	g1->GetHistogram()->SetXTitle("#tau (ns)");
	g1->GetHistogram()->SetYTitle("#sigma (ADC)");
	double minres = 1.e9;
	double mintau = 0.0;
	for(int j=0; j<2000; j++){
	  double tau = ct[i].front()+j*(ct[i].back()-ct[i].front())/2000;
	  double res = g1->Eval(tau, 0, "S");
	  if(res < minres){
	    minres = res;
	    mintau = tau;
	  }
	}
	ct1_mean[i] = mintau;
	gresct1[i] = g1;
	cout << "Ch " << ch << " CT method 1 min resolution of "
	     << minres << " ADC at tau " << mintau << " ns" << endl; 
	for(int j=0; j<(int)hE2[i].size(); j++){
	  pair<double, double> pr;
	  if(thpeak.find(ch) == thpeak.end())
	    pr = GetThRes(hE2[i][j], make_pair(2610., 2620.),
			  make_pair(2570., 2580.), make_pair(2650., 2660.));
	  else
	    pr = GetThRes(hE2[i][j], make_pair(2610., 2620.),
			  make_pair(2570., 2580.), make_pair(2650., 2660.),
			  thpeak[ch]);
	  res2[j] = pr.first;
	  res2_uncert[j] = pr.second;
	}
	for(int j=0; j<(int)res2.size(); j++)
	  if(abs(res2_uncert[j] / res2[j]) > 0.2){
	    res2.erase(res2.begin()+j);
	    dfrac[i].erase(dfrac[i].begin()+j);
	    res2_uncert.erase(res2_uncert.begin()+j);
	    dfrac_uncert.erase(dfrac_uncert.begin()+j);
	  }
	TGraphErrors* g2 = new TGraphErrors(res2.size(), &dfrac[i].front(),
					    &res2.front(),
					    &dfrac_uncert.front(),
					    &res2_uncert.front());
	g2->SetName(("gresct2_"+to_string(i)).c_str());
	g2->SetTitle("");
	g2->GetHistogram()->SetXTitle("D");
	g2->GetHistogram()->SetYTitle("#sigma (ADC)");
	minres = 1.e9;
	double mind = 0.0;
	for(int j=0; j<2000; j++){
	  double d =dfrac[i].front()+j*(dfrac[i].back()-dfrac[i].front())/2000;
	  double res = g2->Eval(d, 0, "S");
	  if(res < minres){
	    minres = res;
	    mind = d;
	  }
	}
	ct2_mean[i] = mind;
	gresct2[i] = g2;
	cout << "Ch " << ch << " CT method 2 min resolution of "
	     << minres << " ADC at c " << mind << endl;
      }
    }
    for(auto const& v : hE1) for(auto const& h : v) if(h) delete h;
    for(auto const& v : hE2) for(auto const& h : v) if(h) delete h;
  }
  
  // grab the calibration parameters if reading from an input configuration
  // otherwise compute them from the input data
  for(auto const& pr : chan_map){
    if(find(chan_skip.begin(), chan_skip.end(), pr.first) != chan_skip.end()){
      cout << "skipping channel " << pr.first << endl;
      continue;
    }
    else if(find(chan_cal.begin(), chan_cal.end(), pr.first)!=chan_cal.end()){
      if(!json_config){
	cout << "channel " << pr.first << " has not been specified "
	     << "for calibration, and no configuration file" << endl;
	return 4;
      }
      else{
	Json::Value value = jvalue[serial_numbers[pr.second]];
	SetJson(value, "base_mean",   base_mean[pr.second]);
	SetJson(value, "base_uncert", base_uncert[pr.second]);
	SetJson(value, "pz_mean",     pz_mean[pr.second]);
	SetJson(value, "ct1_mean",    ct1_mean[pr.second]);
	SetJson(value, "ct2_mean",    ct2_mean[pr.second]);
	SetJson(value, "ct_method",   ct_method[pr.second]);
	SetJson(value, "escale",      escale[pr.second]);
	SetJson(value, "efscale",     efscale[pr.second]);
	SetJson(value, "ecscale",     ecscale[pr.second]);
	SetJson(value, "eoffset",     eoffset[pr.second]);
	SetJson(value, "efoffset",    efoffset[pr.second]);
	SetJson(value, "ecoffset",    ecoffset[pr.second]);
	SetJson(value, "avse_p0",     avse_param[pr.second][0]);
	SetJson(value, "avse_p1",     avse_param[pr.second][1]);
	SetJson(value, "avse_p2",     avse_param[pr.second][2]);
	SetJson(value, "avse_j",      avse_param[pr.second][3]);
	SetJson(value, "aoe_min",     aoe_min[pr.second]);
	//SetJson(value, "rise_t",      trise_val[pr.second]);
	//SetJson(value, "rise_m",      trise_slope[pr.second]);
	SetJson(value, "dcre_slope",  dcre_slope[pr.second]);
	SetJson(value, "dcr_cut_lo",  dcr_cut_lo[pr.second]);
	SetJson(value, "dcr_cut_hi",  dcr_cut_hi[pr.second]);
	for(int i=0; i<(int)chan_map.size(); i++){
	  if(henergyf[i] == NULL) continue;
	  int nebins = henergyf[i]->FindBin((1.e4-efoffset[i])/efscale[i]);
	  double emin = efoffset[i];
	  double emax = henergyf[i]->GetXaxis()->GetBinUpEdge(nebins*efscale[i]+
							      efoffset[i]);
	  if(!use_fixedt){
	    nebins = henergy[i]->FindBin((1.e4-eoffset[i])/escale[i]);
	    emin = eoffset[i];
	    emax=henergy[i]->GetXaxis()->GetBinUpEdge(nebins*escale[i]+eoffset[i]);
	  }
	  hecal[i] = new TH1D(("hecal_"+to_string(i)).c_str(), "",
			      nebins, emin, emax);
	  for(int bin=1; bin<=nebins; bin++){
	    if(use_fixedt)
	      hecal[i]->SetBinContent(bin, henergyf[i]->GetBinContent(bin));
	    else
	      hecal[i]->SetBinContent(bin, henergy[i]->GetBinContent(bin));
	  }
	  hecal[i]->SetXTitle("Energy (keV)");
	  hecal[i]->SetYTitle("Entries");
	}
      }
    }
    else{
      int i = pr.second;
      cout << "getting initial energy scale calibration for channel "
	   << pr.first << endl;
      // initial energy scale correction
      if(thpeak.find(pr.first) == thpeak.end()){
	Th_scale(henergy[i], eoffset[i], escale[i]);
	Th_scale(henergyf[i], efoffset[i], efscale[i]);
	if(ct_method[i] == 1)
	  Th_scale(henergyc1[i], ecoffset[i], ecscale[i]);
	else
	  Th_scale(henergyc2[i], ecoffset[i], ecscale[i]);
      }
      else{
	Th_scale(henergy[i], eoffset[i], escale[i], thpeak[pr.first]);
	Th_scale(henergyf[i], efoffset[i], efscale[i], thpeak[pr.first]);
	if(ct_method[i] == 1)
	  Th_scale(henergyc1[i], ecoffset[i], ecscale[i], thpeak[pr.first]);
	else
	  Th_scale(henergyc2[i], ecoffset[i], ecscale[i], thpeak[pr.first]);
      }
      int nebins = henergyc1[i]->FindBin((1.e4-ecoffset[i])/ecscale[i]);
      double emin = ecoffset[i];
      double emax = henergyc1[i]->GetXaxis()->GetBinUpEdge(nebins*ecscale[i]+
							  ecoffset[i]);
      if(!get_ct_decay){
	nebins = henergyf[i]->FindBin((1.e4-efoffset[i])/escale[i]);
	emin = efoffset[i];
	emax=henergyf[i]->GetXaxis()->GetBinUpEdge(nebins*efscale[i]+
						   efoffset[i]);
      }
      if(!use_fixedt){
	nebins = henergy[i]->FindBin((1.e4-eoffset[i])/escale[i]);
	emin = eoffset[i];
	emax=henergy[i]->GetXaxis()->GetBinUpEdge(nebins*escale[i]+eoffset[i]);
      }
      hecal[i] = new TH1D(("hecal_"+to_string(i)).c_str(), "",
			  nebins, emin, emax);
      for(int bin=1; bin<=nebins; bin++){
	if(get_ct_decay){
	  if(ct_method[i] == 1)
	    hecal[i]->SetBinContent(bin, henergyc1[i]->GetBinContent(bin));
	  else
	    hecal[i]->SetBinContent(bin, henergyc2[i]->GetBinContent(bin));
	}
	else if(use_fixedt)
	  hecal[i]->SetBinContent(bin, henergyf[i]->GetBinContent(bin));
	else
	  hecal[i]->SetBinContent(bin, henergy[i]->GetBinContent(bin));
      }
      hecal[i]->SetXTitle("Energy (keV)");
      hecal[i]->SetYTitle("Entries");
      // baseline, baseline rms, and pz constant
      int ebin0 =
	hbase_energy[i]->GetXaxis()->FindBin((200-eoffset[i])/escale[i]);
      int ebin1 =
	hbase_energy[i]->GetXaxis()->FindBin((5000-eoffset[i])/escale[i]);
      hbase[i] =  hbase_energy[i]->ProjectionY(("hbase_"+to_string(i)).c_str(),
					       ebin0, ebin1);
      hbrms[i] =  hbrms_energy[i]->ProjectionY(("hbrms_"+to_string(i)).c_str(),
					       ebin0, ebin1);
      base_mean[i] = hbase[i]->GetMean();
      brms_mean[i] = hbrms[i]->GetMean();
      base_uncert[i]  = hbase[i]->GetMeanError();
      brms_uncert[i]  = hbrms[i]->GetMeanError();
      hdecay[i]=hdecay_energy[i]->ProjectionY(("hdecay_"+to_string(i)).c_str(),
						ebin0, ebin1);
      if(hdecay[i]->GetMaximum() > 0.9*hdecay[i]->Integral()){
	pz_mean[i] = hdecay[i]->GetBinCenter(hdecay[i]->GetMaximumBin());
	pz_uncert[i] = 0.0;
      }
      else{
	double pzm = hdecay[i]->GetBinCenter(hdecay[i]->GetMaximumBin());
	double pzr = hdecay[i]->GetRMS();
	TF1* fpz = new TF1("fpz", "gaus(0)", pzm-0.5*pzr, pzm+0.5*pzr);
	hdecay[i]->Fit(fpz, "QR+");
	pz_mean[i] = fpz->GetParameter(1);
	pz_uncert[i] = fpz->GetParError(1);
	delete fpz;
      }
      // dcr slope correction
      TProfile* dcr_prof = hdcr_energy[i]->ProfileX(("hdcr_prof_" +
						     to_string(i)).c_str());
      TF1* fdcr = new TF1("fdcr", "[0]+[1]*x",
			  0.0, dcr_prof->GetXaxis()->GetXmax());
      fdcr->FixParameter(0, 0.0);
      fdcr->SetParameter(1, -1.0e5);
      dcr_prof->Fit(fdcr, "QFR+");
      dcre_slope[i]   = fdcr->GetParameter(1);
      dcre_uncert[i] = fdcr->GetParError(1);
    }
  }

  // calibrate AvsE, A/E,  and DCR
  cout << "calibrating avse and dcr" << endl;
  for(auto const& ch : chan_cal){
    int ich = chan_map[ch];
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
    haoes_energy[ich] = new TH2D(("haoe_energy_"+to_string(ich)).c_str(), "",
				 hecal[ich]->GetXaxis()->GetNbins(),
				 hecal[ich]->GetXaxis()->GetXmin(),
				 hecal[ich]->GetXaxis()->GetXmax(),
				 haoe_energy[ich]->GetYaxis()->GetNbins(),
				 haoe_energy[ich]->GetYaxis()->GetXmin(),
				 haoe_energy[ich]->GetYaxis()->GetXmax());
    haoes_energy[ich]->SetXTitle("Energy (keV)");
    haoes_energy[ich]->SetYTitle("A/E");
    hdcrs_energy[ich] = new TH2D(("hdcrs_energy_"+to_string(ich)).c_str(),"",
				 hecal[ich]->GetXaxis()->GetNbins(),
				 hecal[ich]->GetXaxis()->GetXmin(),
				 hecal[ich]->GetXaxis()->GetXmax(),
				 hdcr_energy[ich]->GetYaxis()->GetNbins(),
				 hdcr_energy[ich]->GetYaxis()->GetXmin(),
				 hdcr_energy[ich]->GetYaxis()->GetXmax());
    hdcrs_energy[ich]->SetXTitle("Energy (keV)");
    hdcrs_energy[ich]->SetYTitle("DCR Slope");
  }
  
  // populate the histograms
  reader.Restart();
  intree->SetBranchStatus("ct_decay", false);
  intree->SetBranchStatus("ct_value", false);
  intree->SetBranchStatus("ct_integral", false);
  intree->SetBranchStatus("imax", true);
  intree->SetBranchStatus("ct1_trappick", true);
  intree->SetBranchStatus("ct2_trappick", true);
  intree->SetBranchStatus("trapmax", true);
  intree->SetBranchStatus("dcrslope", true);
  while(reader.Next()){
    for(int ich=0; ich<(int)channel->size(); ich++){
      if(find(chan_cal.begin(),
	      chan_cal.end(), channel->at(ich)) == chan_cal.end())
	continue;
      int i = chan_map[channel->at(ich)];
      if(get_ct_decay){
	double E = ct1_trappick->at(ich)*ecscale[i] + ecoffset[i];
	double ct = ct1_trappick->at(ich);
	if(ct_method[i] == 2){
	  E = ct2_trappick->at(ich)*ecscale[i] + ecoffset[i];
	  ct = ct2_trappick->at(ich);
	}
	hamps_energy[i]->Fill(E, imax->at(ich)*E/ct);
	haoes_energy[i]->Fill(E, imax->at(ich)/ct);
	hdcrs_energy[i]->Fill(E, dcrslope->at(ich)-dcre_slope[i]*ct);
      }
      else if(use_fixedt){
	double E = trappick->at(ich)*efscale[i] + efoffset[i];
	hamps_energy[i]->Fill(E, imax->at(ich)*E/trappick->at(ich));
	haoes_energy[i]->Fill(E, imax->at(ich)/trappick->at(ich));
	hdcrs_energy[i]->Fill(E, dcrslope->at(ich)-
			      dcre_slope[i]*trappick->at(ich));
      }
      else{
	double E = trapmax->at(ich)*escale[i] + eoffset[i];
	hamps_energy[i]->Fill(E, imax->at(ich)*E/trapmax->at(ich));
	haoes_energy[i]->Fill(E, imax->at(ich)/trapmax->at(ich));
	hdcrs_energy[i]->Fill(E, dcrslope->at(ich)-
			      dcre_slope[i]*trapmax->at(ich)); 
      }
    }
  }
  
  // calibrate AvsE and DCR
  for(auto const& ch : chan_cal){
    int i = chan_map[ch];
    gavse[i] = AvsECal(hamps_energy[i],  avse_param[i], avse_uncert[i]);
    AoECal(haoes_energy[i], aoe_min[i], aoe_min_uncert[i]);
    DCRCal(hdcrs_energy[i], 0.01, 0.01, dcr_cut_lo[i], dcr_cut_hi[i]);
  }
  
  // for now just plot the rise time, no energy correction based on it
  //cout << "computing rise time correction" << endl;
  reader.Restart();
  intree->SetBranchStatus("t1", true);
  intree->SetBranchStatus("t99", true);
  intree->SetBranchStatus("sampling", true);
  for(auto const& ch : chan_cal){
    int i= chan_map[ch];
    hrises_energy[i] = new TH2D(("hrises_energy_" +
				 to_string(i)).c_str(), "",
				hecal[i]->GetXaxis()->GetNbins(),
				hecal[i]->GetXaxis()->GetXmin(),
				hecal[i]->GetXaxis()->GetXmax(),
				hrise_energy[i]->GetYaxis()->GetNbins(),
				hrise_energy[i]->GetYaxis()->GetXmin(),
				hrise_energy[i]->GetYaxis()->GetXmax());
    hrises_energy[i]->SetXTitle("Energy (keV)");
    hrises_energy[i]->SetYTitle("Rise Time (ns)");
  }
  while(reader.Next())
    for(int ich=0; ich<(int)channel->size(); ich++){
      if(find(chan_cal.begin(),
	      chan_cal.end(), channel->at(ich)) == chan_cal.end())
	continue;
      int i = chan_map[channel->at(ich)];
      double E = ct1_trappick->at(ich)*ecscale[i] + ecoffset[i];
      double avse = -imax->at(ich)*E/ct1_trappick->at(ich);
      if(ct_method[i] == 2){
	E = ct2_trappick->at(ich)*ecscale[i] + ecoffset[i];
	avse = -imax->at(ich)*E/ct2_trappick->at(ich);
      }
      if(!get_ct_decay){
	E = trappick->at(ich)*efscale[i] + efoffset[i];
	avse = -imax->at(ich)*E/trappick->at(ich);
      }
      if(!use_fixedt){
	E = trapmax->at(ich)*escale[i] + eoffset[i];
	avse = -imax->at(ich)*E/trapmax->at(ich);
      }
      avse += avse_param[i][0] + avse_param[i][1]*E;
      avse += avse_param[i][2]*pow(E, 2);
      if(avse/avse_param[i][3] <= -1.0) continue;
      hrises_energy[i]->Fill(E, (t99->at(ich)-
				 t1->at(ich))*sampling->at(ich));
    }
  /*
  for(int ich=0; ich<(int)hrises_energy.size(); ich++){
    vector<double> rise_slope;
    vector<double> rise_error;
    vector<double> euncert(tl_peaks.size(), 1.0);
    TAxis* axis = hrises_energy[ich]->GetYaxis();
    TProfile* hproj = hrises_energy[ich]->ProfileX(("hrises_energy_proj_" +
						    to_string(ich)).c_str(),
						   axis->FindBin(30.0),
						   axis->FindBin(60.0));
    TF1* fflat = new TF1(("fflat_"+to_string(ich)).c_str(), "[0]", 300,3000);
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
  */
    
  // print the calibration parameters
  for(auto const& ch : chan_cal){
    int i = chan_map[ch];
    cout << "channel " << ch << ":  " << endl;
    print_value("base    ", 3, base_mean[i],     base_uncert[i]);
    print_value("brms    ", 3, brms_mean[i],     brms_uncert[i]);
    print_value("pz      ", 3, pz_mean[i],       pz_uncert[i]);
    print_value("ct tau  ", 3, ct1_mean[i],      ct1_uncert[i]);
    print_value("ct frac ", 6, ct2_mean[i],      ct2_uncert[i]);
    print_value("avse p0 ", 3, avse_param[i][0], avse_uncert[i][0]);
    print_value("avse p1 ", 3, avse_param[i][1], avse_uncert[i][1]);
    print_value("avse p2 ", 3, avse_param[i][2], avse_uncert[i][2]);
    print_value("avse j  ", 3, avse_param[i][3], avse_uncert[i][3]);
    print_value("aoe min ", 3, aoe_min[i],       aoe_min_uncert[i]);
    //print_value("rise t  ", 3, trise_val[i],     trise_uncert[i]);
    //print_value("rise m  ", 3, trise_slope[i],   trise_suncert[i]);
    print_value("dcr  m  ", 3, dcre_slope[i],    dcre_uncert[i]);
    print_value("dcr  lo ", 3, dcr_cut_lo[i],    0.0);
    print_value("dcr  hi ", 3, dcr_cut_hi[i],    0.0);
  }

  // fixme - copy over other useful values from the below
  intree->SetBranchStatus("*", true);
  reader.Restart();
  TTreeReaderValue<int> runIn(reader, "run");
  TTreeReaderValue<int> eventnumIn(reader, "eventnum");
  TTreeReaderValue<vector<string> >  detserialIn(reader, "detserial");
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
  //TTreeReaderValue<vector<vector<double> > > exp_param(reader, "exp_param");

  // output tree setup
  int run, eventnum;
  vector<string> detserial;
  vector<int> chan;
  vector<double> baseSigma, nbaserms, timestamp, dt, T0, trise;
  vector<double> trapECal, trapEFCal, trapEFCCal,/*trapEFRCal,*/avse, aoe, dcr;
  TFile* outfile = new TFile(ofname.c_str(), "recreate");
  tdir->cd();
  TTree* outtree = new TTree("tree", "tree");
  outtree->Branch("run", &run);
  outtree->Branch("eventnum", &eventnum);
  outtree->Branch("detserial", &detserial);
  outtree->Branch("channel", &chan);
  outtree->Branch("baseSigma", &baseSigma);
  outtree->Branch("nbaserms", &nbaserms);
  outtree->Branch("time", &timestamp);
  outtree->Branch("deltat", &dt);
  outtree->Branch("t0", &T0);
  outtree->Branch("trise", &trise);
  outtree->Branch("trapECal", &trapECal);
  outtree->Branch("trapEFCal", &trapEFCal);
  outtree->Branch("trapEFCCal", &trapEFCCal);
  //outtree->Branch("trapEFRCal", &trapEFRCal);
  outtree->Branch("avse", &avse);
  outtree->Branch("aoe", &aoe);
  outtree->Branch("dcr", &dcr);
  outtree->SetDirectory(outfile);

  // apply the calibration parameters to the relevant input values
  cout << "applying calibration parameters to input trees" << endl;
  int iev = -1;
  while(reader.Next()){
    iev ++;
    if(iev % 100000 == 0 && iev > 0) outtree->AutoSave();
    run = *runIn;
    eventnum = *eventnumIn;
    int nwf = (int) channel->size();
    detserial.assign(nwf, "");
    chan.assign(nwf, 0);
    baseSigma.assign(nwf, 0.0);
    nbaserms.assign(nwf, 0.0);
    timestamp.assign(nwf, 0.0);
    dt.assign(nwf, 0.0);
    T0.assign(nwf, 0.0);
    trise.assign(nwf, 0.0);
    trapECal.assign(nwf, 0.0);
    trapEFCal.assign(nwf, 0.0);
    trapEFCCal.assign(nwf, 0.0);
    //trapEFRCal.assign(nwf, 0.0);
    avse.assign(nwf, 0.0);
    aoe.assign(nwf, 0.0);
    dcr.assign(nwf, 0.0);
    for(int ich=0; ich<nwf; ich++){
      detserial[ich] = detserialIn->at(ich);
      chan[ich] = channel->at(ich);
      if(chan_map.find(chan[ich]) == chan_map.end() ||
	 find(chan_skip.begin(), chan_skip.end(), chan[ich])!=chan_skip.end())
	continue;
      int i = chan_map[chan[ich]];
      baseSigma[ich] = (baseline->at(ich)-base_mean[i])/brms_mean[i];
      nbaserms[ich] = baserms->at(ich)/brms_mean[i];
      timestamp[ich] = times->at(ich);
      T0[ich] = t0->at(ich);
      dt[ich] = deltat->at(ich);
      trise[ich] = t99->at(ich)-t1->at(ich);
      trapECal[ich] = trapmax->at(ich)*escale[i]+eoffset[i];
      trapEFCal[ich] = trappick->at(ich)*efscale[i]+efoffset[i];
      if(ct_method[i] == 1)
	trapEFCCal[ich] = ct1_trappick->at(ich)*ecscale[i]+ecoffset[i];
      else if(ct_method[i] == 2)
	trapEFCCal[ich] = ct2_trappick->at(ich)*ecscale[i]+ecoffset[i];
      //trapEFRCal[i] = trapEFCal[i] + ((trise_val[i] -
      //                                 trise[i])) / trise_slope[i];
      double E = trapEFCCal[ich];
      if(!get_ct_decay) E = trapEFCal[ich];
      if(!use_fixedt) E = trapECal[ich];
      
      avse[ich] = avse_param[i][0] + avse_param[i][1]*E;
      avse[ich] += avse_param[i][2]*pow(E, 2);
      if(get_ct_decay){
	if(ct_method[i] == 1){
	  avse[ich] -= imax->at(ich)*trapEFCCal[ich]/ct1_trappick->at(ich);
	  aoe[ich] = imax->at(ich) / ct1_trappick->at(ich);
	}
	else if(ct_method[i] == 2){
	  avse[ich] -= imax->at(ich)*trapEFCCal[ich]/ct2_trappick->at(ich);
	  aoe[ich] = imax->at(ich) / ct2_trappick->at(ich);
	}
      }
      else if(use_fixedt){
	avse[ich] -= imax->at(ich)*trapEFCal[ich]/trappick->at(ich);
	aoe[ich] = imax->at(ich) / trappick->at(ich);
      }
      else{
	avse[ich] -= imax->at(ich)*trapECal[ich]/trapmax->at(ich);
	aoe[ich] = imax->at(ich) / trapmax->at(ich);
      }
      avse[ich] /= avse_param[i][3];
      aoe[ich] -= aoe_min[i] - 1;
      // fixme - this places the cuts at -1 to 1, probably not what we want
      if(get_ct_decay)
	dcr[ich] = dcrslope->at(ich);//-dcre_slope[i]*ct_trappick->at(ich);
      else if(use_fixedt)
	dcr[ich] = dcrslope->at(ich);//-dcre_slope[i]*trappick->at(ich);
      else
	dcr[ich] = dcrslope->at(ich);//-dcre_slope[i]*trapmax->at(ich);
      //dcr[ich] = dcr[ich]*2/(dcr_cut_hi[i]-dcr_cut_lo[i]);
      //dcr[ich] += (dcr_cut_lo[i]+dcr_cut_hi[i])/
      //(dcr_cut_lo[i]-dcr_cut_hi[i]);
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
    if(henergyc1[i])
      henergyc1[i]->Write(("henergyc1_"+to_string(p.first)).c_str());
    if(henergyc2[i])
      henergyc2[i]->Write(("henergyc2_"+to_string(p.first)).c_str());
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
      hamps_energy[i]->Write(("hamps_energy_"+to_string(p.first)).c_str());
    if(haoe_energy[i])
      haoe_energy[i]->Write(("haoe_energy_"+to_string(p.first)).c_str());
    if(haoes_energy[i])
      haoes_energy[i]->Write(("haoes_energy_"+to_string(p.first)).c_str());
    if(hdcr_energy[i])
      hdcr_energy[i]->Write(("hdcr_energy_"+to_string(p.first)).c_str());
    if(hdcrs_energy[i])
      hdcrs_energy[i]->Write(("hdcrs_energy_"+to_string(p.first)).c_str());
    if(hrise_energy[i])
      hrise_energy[i]->Write(("hrise_energy_"+to_string(p.first)).c_str());
    if(hrises_energy[i])
      hrises_energy[i]->Write((("hrises_energy_"+to_string(p.first)).c_str()));
    if(hrise_slope[i]) hrise_slope[i]->Write();
    if(hrises_energy_proj[i]) hrises_energy_proj[i]->Write();
    if(henergy[i]) henergy[i]->Write();
    if(hecal[i]) hecal[i]->Write();
    if(gavse[i]) gavse[i]->Write();
    if(gresct1[i]) gresct1[i]->Write();
    if(gresct2[i]) gresct2[i]->Write();
  }
  outfile->WriteObject(&chan_map,    "chan_map");
  outfile->WriteObject(&base_mean,   "base_mean");
  outfile->WriteObject(&base_uncert, "base_uncert");
  outfile->WriteObject(&brms_mean,   "brms_mean");
  outfile->WriteObject(&brms_uncert, "brms_uncert");
  outfile->WriteObject(&pz_mean,     "pz_mean");
  outfile->WriteObject(&pz_uncert,   "pz_uncert");
  outfile->WriteObject(&ct1_mean,    "ct1_mean");
  outfile->WriteObject(&ct2_mean,    "ct2_mean");
  
  // write calibration parameters to a json file if not reading from one
  cout << "writing calibration parameters to " << jfile << endl;
  for(auto const& ch : chan_cal){
    int i = chan_map[ch];
    string s = serial_numbers[i];
    jvalue[s]["base_mean"]  = base_mean[i];
    jvalue[s]["base_rms"]   = brms_mean[i];
    jvalue[s]["pz_mean"]    = pz_mean[i];
    jvalue[s]["ct1_mean"]   = ct1_mean[i];
    jvalue[s]["ct2_mean"]   = ct2_mean[i];
    jvalue[s]["ct_method"]  = ct_method[i];
    jvalue[s]["escale"]     = escale[i];
    jvalue[s]["eoffset"]    = eoffset[i];
    jvalue[s]["efscale"]    = efscale[i];
    jvalue[s]["efoffset"]   = efoffset[i];
    jvalue[s]["ecscale"]    = ecscale[i];
    jvalue[s]["ecoffset"]   = ecscale[i];
    jvalue[s]["avse_p0"]    = avse_param[i][0];
    jvalue[s]["avse_p1"]    = avse_param[i][1];
    jvalue[s]["avse_p2"]    = avse_param[i][2];
    jvalue[s]["avse_j"]     = avse_param[i][3];
    jvalue[s]["aoe_min"]    = aoe_min[i];
    //jvalue[s]["rise_t"]     = trise_val[i];
    //jvalue[s]["rise_m"]     = trise_slope[i];
    jvalue[s]["dcre_slope"] = dcre_slope[i];
    jvalue[s]["dcr_cut_lo"] = dcr_cut_lo[i];
    jvalue[s]["dcr_cut_hi"] = dcr_cut_hi[i];
  }
  Json::StreamWriterBuilder builder;
  builder["indentation"] = "";
  TObjString* jout=new TObjString(Json::writeString(builder, jvalue).c_str());
  jout->Write("calibration_config");
  outfile->Close();
  builder["commentStyle"] = "None";
  builder["indentation"]  = "    ";
  unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
  ofstream fstream(jfile);
  writer->write(jvalue, &fstream);
  fstream.close();
  
  return 0;
}
