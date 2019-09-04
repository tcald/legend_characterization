#include "utils.hh"
#include <GATPeakShape.hh>
#include <GATHybridMonteCarlo.hh>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLine.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <json/json.h>
#include <json/value.h>
#include <json/reader.h>
#include <json/writer.h>
#include <utility>
#include <algorithm>
#include <getopt.h>
#include <assert.h>
#include <sstream>

using namespace std;

TDirectory* tdir;

int main(int argc, char* argv[]){

  // root setup
  gStyle->SetOptStat(0);
  gStyle->SetLabelFont(132, "XYZ");
  gStyle->SetTitleFont(132, "XYZ");
  tdir = gROOT->CurrentDirectory();

  // peak locations and fit ranges
  vector<double> th_peaks({277.0, 300.0, 583.0, 727.0, 861.0, 2615.0});
  vector<pair<double, double> > th_ranges;
  th_ranges.push_back(make_pair(260.0, 285.0));
  th_ranges.push_back(make_pair(290.0, 310.0));
  th_ranges.push_back(make_pair(570.0, 595.0));
  th_ranges.push_back(make_pair(710.0, 740.0));
  th_ranges.push_back(make_pair(845.0, 870.0));
  th_ranges.push_back(make_pair(2585.0, 2635.0));
  
  // default histogram binning
  int ebins1 = 12000;
  int ebins2 = 3000;
  double emin = 0.0;
  double emax = 3000.0;
  int dcrbins1 = 500;
  int dcrbins2 = 500;
  double dcrmin = -5.0;
  double dcrmax = 10.0;
  int aebins1 = 1000;
  int aebins2 = 1000;
  double avsemin = -15.0;
  double avsemax = 5.0;
  double aoemin = -2.0;
  double aoemax = -0.5;

  // option handling
  vector<string> infname;
  string outfname = "";
  map<int, int> chan_map;
  bool all_chan = true;
  vector<string> sources;
  vector<double> peaks;
  vector<pair<double, double> > ranges;
  bool do_psa = false;
  string write_dir = "";
  bool json_config = false;
  Json::Value jvalue;
  static struct option opts[]{
    {"help",           no_argument, NULL, 'h'},
    {"infile",   required_argument, NULL, 'i'},
    {"outfile",  required_argument, NULL, 'o'},
    {"channel",  required_argument, NULL, 'c'},
    {"source",   required_argument, NULL, 's'},
    {"writedir", required_argument, NULL, 'w'},
    {"jsonfile", required_argument, NULL, 'j'}
  };
  int opt = getopt_long(argc, argv, "hi:o:c:w:j:", opts, NULL);
  while(opt != -1){
    switch(opt){
    case 'h':
      cout << "options:"             << endl;
      cout << "  -i input filename"  << endl;
      cout << "  -o output filename" << endl;
      cout << "  -c channel number"  << endl;
      return 0;
    case 'i': infname.push_back(string(optarg)); break;
    case 'o': outfname = string(optarg);         break;
    case 'c':{
      unsigned i = chan_map.size();
      chan_map[atoi(optarg)] = (int) i;
      all_chan = false;
      break;
    }
    case 's':{
      if(find(sources.begin(), sources.end(),
	      string(optarg)) != sources.end()) break;
      sources.push_back(string(optarg));
      if(string(optarg) == "Th"){
	peaks.insert(peaks.end(), th_peaks.begin(), th_peaks.end());
	ranges.insert(ranges.end(), th_ranges.begin(), th_ranges.end());
	do_psa = true;
      }
      else cout << "unrecognized source type " << optarg << endl;
      break;
    }
    case 'w': write_dir = string(optarg); break;
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
    default: return 1;
    }
    opt = getopt_long(argc, argv, "hi:o:c:w:j:", opts, NULL);
  }
  assert(infname.size() != 0 && outfname != "");
  if(peaks.size() == 0){
    cout << "no source type specified, assuming 228Th" << endl;
    sources.push_back("Th");
    peaks.insert(peaks.end(), th_peaks.begin(), th_peaks.end());
    ranges.insert(ranges.end(), th_ranges.begin(), th_ranges.end());
    do_psa = true;
  }
  
  // read the histogram binning parameters from the json configuration
  if(json_config){
    SetJson(jvalue, "ebins1",   ebins1);
    SetJson(jvalue, "ebins2",   ebins2);
    SetJson(jvalue, "emin",     emin);
    SetJson(jvalue, "emax",     emax);
    SetJson(jvalue, "dcrbins1", dcrbins1);
    SetJson(jvalue, "dcrbins2", dcrbins2);
    SetJson(jvalue, "dcrmin",   dcrmin);
    SetJson(jvalue, "dcrmax",   dcrmax);
    SetJson(jvalue, "aebins1",  aebins1);
    SetJson(jvalue, "aebins2",  aebins2);
    SetJson(jvalue, "avsemin",  avsemin);
    SetJson(jvalue, "avsemax",  avsemax);
    SetJson(jvalue, "aoemin",   aoemin);
    SetJson(jvalue, "aoemax",   aoemax);
  }

  // template histograms
  stringstream ssebin;
  ssebin << fixed << setprecision(2) << (emax-emin)/ebins1;
  string selabel = "Events / " + ssebin.str() + " keV";
  TH1D* hE = new TH1D("hE", "", ebins1, emin, emax);
  hE->SetXTitle("Energy (keV)");
  hE->SetYTitle(selabel.c_str());
  TH1D* hdcr = new TH1D("hdcr", "", dcrbins1, dcrmin, dcrmax);
  hdcr->SetXTitle("DCR");
  hdcr->SetYTitle("Events");
  TH2D* hdcrEt = new TH2D("hdcrE", "",
			 ebins2, emin, emax,
			 dcrbins2, dcrmin, dcrmax);
  hdcrEt->SetXTitle("Energy (keV)");
  hdcrEt->SetYTitle("DCR");
  TH1D* havse = new TH1D("havse", "", aebins1, avsemin, avsemax);
  havse->SetXTitle("AvsE");
  havse->SetYTitle("Events");
  TH2D* havseEt = new TH2D("havseE", "",
			  ebins2, emin, emax,
			  aebins2, avsemin, avsemax);
  havseEt->SetXTitle("Energy (keV)");
  havseEt->SetYTitle("AvsE");
  TH1D* haoe = new TH1D("haoe", "", aebins1, aoemin, aoemax);
  haoe->SetXTitle("A/E");
  haoe->SetYTitle("Events");
  TH2D* haoeEt = new TH2D("haoeE", "",
			 ebins2, emin, emax,
			 aebins2, aoemin, aoemax);
  haoeEt->SetXTitle("Energy (keV)");
  haoeEt->SetYTitle("A/E");

  // histograms for each channel
  vector<TH1D*> hECal, hEFCal, hEFCCal, hEFCal_avse, hEFCal_aoe;
  vector<TH1D*> hdcr2615, havse2615, haoe2615;
  vector<TH2D*> hdcrE, havseE, haoeE;
  vector<vector<TH1D*> > hEpeak;
  for(auto const& pr : chan_map){
    string sch = to_string(pr.second);
    hECal.push_back((TH1D*) hE->Clone(("hECal_"+sch).c_str()));
    hEFCal.push_back((TH1D*) hE->Clone(("hEFCal_"+sch).c_str()));
    hEFCCal.push_back((TH1D*) hE->Clone(("heFCCal_"+sch).c_str()));
    hEFCal_avse.push_back((TH1D*) hE->Clone(("hEFCal_avse_"+sch).c_str()));
    hEFCal_aoe.push_back((TH1D*) hE->Clone(("hEFCal_aoe_"+sch).c_str()));
    hdcr2615.push_back((TH1D*) hdcr->Clone(("hdcr2615_"+sch).c_str()));
    havse2615.push_back((TH1D*) havse->Clone(("havse2615_"+sch).c_str()));
    haoe2615.push_back((TH1D*)  haoe->Clone(("haoe2615_"+sch).c_str()));
    hdcrE.push_back((TH2D*) hdcrEt->Clone(("hdcrE_"+sch).c_str()));
    havseE.push_back((TH2D*) havseEt->Clone(("havseE_"+sch).c_str()));
    haoeE.push_back((TH2D*) haoeEt->Clone(("haoeE_"+sch).c_str()));
    vector<TH1D*> vh;
    for(int i=0; i<(int)peaks.size(); i++)
      vh.push_back(new TH1D(("hEFCal_"+to_string((int)peaks[i])).c_str(), "",
			    8*((int)(ranges[i].second-ranges[i].first)),
			    ranges[i].first, ranges[i].second));
    hEpeak.push_back(vh);
  }

  // input tree setup
  TChain* tree = new TChain("tree", "tree");
  for(auto const& f : infname) tree->Add(f.c_str());
  tree->SetBranchStatus("*", false);
  tree->SetBranchStatus("channel", true);
  tree->SetBranchStatus("detserial", true);
  tree->SetBranchStatus("trapECal", true);
  tree->SetBranchStatus("trapEFCal", true);
  tree->SetBranchStatus("trapEFCCal", true);
  tree->SetBranchStatus("dcr", true);
  tree->SetBranchStatus("avse", true);
  tree->SetBranchStatus("aoe", true);
  tree->SetBranchStatus("baseSigma", true);
  tree->SetBranchStatus("nbaserms", true);
  TTreeReader reader(tree);
  TTreeReaderValue<vector<int> > channel(reader, "channel");
  TTreeReaderValue<vector<string> > detserial(reader, "detserial");
  TTreeReaderValue<vector<double> > trapECal(reader, "trapECal");
  TTreeReaderValue<vector<double> > trapEFCal(reader, "trapEFCal");
  TTreeReaderValue<vector<double> > trapEFCCal(reader, "trapEFCCal");
  TTreeReaderValue<vector<double> > dcr(reader, "dcr");
  TTreeReaderValue<vector<double> > avse(reader, "avse");
  TTreeReaderValue<vector<double> > aoe(reader, "aoe");
  TTreeReaderValue<vector<double> > baseSigma(reader, "baseSigma");
  TTreeReaderValue<vector<double> > nbaserms(reader, "nbaserms");

  // populate the histograms, initialize them here if using all channels
  vector<string> det_serial(chan_map.size());
  while(reader.Next())
    for(int ich=0; ich<(int)channel->size(); ich++){
      if(all_chan && chan_map.find(channel->at(ich)) == chan_map.end()){
	int i = (int) chan_map.size();
	chan_map[channel->at(ich)] = i;
	det_serial.resize(chan_map.size());
	string s = to_string(channel->at(ich));
	hECal.push_back((TH1D*) hE->Clone(("hECal_"+s).c_str()));
	hEFCal.push_back((TH1D*) hE->Clone(("hEFCal_"+s).c_str()));
	hEFCCal.push_back((TH1D*) hE->Clone(("heFCCal_"+s).c_str()));
	hEFCal_avse.push_back((TH1D*) hE->Clone(("hEFCal_avse_"+s).c_str()));
	hEFCal_aoe.push_back((TH1D*) hE->Clone(("hEFCal_aoe_"+s).c_str()));
	hdcr2615.push_back((TH1D*) hdcr->Clone(("hdcr2615_"+s).c_str()));
	havse2615.push_back((TH1D*)havse->Clone(("havse2615_"+s).c_str()));
	haoe2615.push_back((TH1D*) haoe->Clone(("haoe2615_"+s).c_str()));
	hdcrE.push_back((TH2D*) hdcrEt->Clone(("hdcrE_"+s).c_str()));
	havseE.push_back((TH2D*) havseEt->Clone(("havseE_"+s).c_str()));
	haoeE.push_back((TH2D*) haoeEt->Clone(("haoeE_"+s).c_str()));
	vector<TH1D*> vh;
	for(int j=0; j<(int)peaks.size(); j++)
	  vh.push_back(new TH1D(("hEFCal_"+to_string((int)peaks[j])).c_str(),
				"",8*((int)(ranges[j].second-ranges[j].first)),
				ranges[j].first, ranges[j].second));
	hEpeak.push_back(vh);
      }
      int i = chan_map[channel->at(ich)];
      if(find(det_serial.begin(), det_serial.end(), detserial->at(ich)) ==
	 det_serial.end())
	det_serial[i] = detserial->at(ich);
      if(abs(baseSigma->at(ich)) > 5) continue;
      if(abs(nbaserms->at(ich)) > 5) continue;
      hdcrE[i]->Fill(trapEFCal->at(ich), dcr->at(ich));
      //if(dcr->at(ich)<-1 || dcr->at(ich)>1) continue;
      hECal[i]->Fill(trapECal->at(ich));
      hEFCal[i]->Fill(trapEFCal->at(ich));
      hEFCCal[i]->Fill(trapEFCCal->at(ich));
      for(auto const& h : hEpeak[i]) h->Fill(trapEFCal->at(ich));
      if(trapEFCal->at(ich) > 2610 && trapEFCal->at(ich) < 2620){
	havse2615[i]->Fill(avse->at(ich));
	haoe2615[i]->Fill(aoe->at(ich));
	hdcr2615[i]->Fill(dcr->at(ich));
      }
      havseE[i]->Fill(trapEFCal->at(ich), avse->at(ich));
      haoeE[i]->Fill(trapEFCal->at(ich), aoe->at(ich));
      if(avse->at(ich) >= -1) hEFCal_avse[i]->Fill(trapEFCal->at(ich));
      if(aoe->at(ich) >= 0) hEFCal_aoe[i]->Fill(trapEFCal->at(ich));
    }

  // clean up template histograms
  delete hE;
  delete hdcr;
  delete havse;
  delete haoe;
  delete hdcrEt;
  delete havseEt;
  delete haoeEt;

  // plot the energy parameters
  vector<TCanvas*> cenergy(chan_map.size(), NULL);
  int ich = -1;
  for(auto const& pr : chan_map){
    ich ++;
    if(!hECal[ich]) continue;
    hECal[ich]->SetLineColor(2);
    hEFCal[ich]->SetLineColor(4);
    hEFCCal[ich]->SetLineColor(1);
    hEFCal_aoe[ich]->SetLineColor(8);
    TLegend* l = new TLegend(0.6, 0.6, 0.8, 0.8);
    l->SetLineColor(0);
    l->SetFillColor(0);
    l->AddEntry(hECal[ich],      "#font[132]{Trap Maximum}", "l");
    l->AddEntry(hEFCal[ich],     "#font[132]{Fixed Time Pickoff}", "l");
    l->AddEntry(hEFCCal[ich],    "#font[132]{CT Correction}", "l");
    l->AddEntry(hEFCal_aoe[ich], "#font[132]{A/E Cut}", "l");
    TCanvas* c = new TCanvas(("cenergy_"+to_string(pr.first)).c_str(), "",
			     0, 0, 1200, 800);
    c->cd(1)->SetLogy(true);
    c->cd(1);
    hEFCCal[ich]->Draw();
    hECal[ich]->Draw("same");
    hEFCal[ich]->Draw("same");
    hEFCal_aoe[ich]->Draw("same");
    l->Draw();
    cenergy[ich] = c;
    tdir->cd();
  }

  // fit the peaks using GATPeakShape
  vector<vector<TF1*> > fpeak;
  vector<vector<double> > detFWHM(hEpeak.size(),
				  vector<double>(peaks.size(), 0.0));
  vector<vector<double> > detFWHMuncert(hEpeak.size(),
					vector<double>(peaks.size(), 0.0));
  vector<vector<double> > detCent(hEpeak.size(),
				  vector<double>(peaks.size(), 0.0));
  vector<vector<double> > detCentUncert(hEpeak.size(),
					vector<double>(peaks.size(), 0.0));
  for(int i=0; i<(int)hEpeak.size(); i++){
    vector<TF1*> vf;
    int j = -1;
    for(auto const& h : hEpeak[i]){
      j ++;
      h->Sumw2();
      GATPeakShape* peak = new GATPeakShape();
      double p = h->GetBinCenter(h->GetMaximumBin());
      double r = sqrt(pow(0.38,2)+p*0.002+pow(p*0.0007,2));
      int bin0 = h->GetXaxis()->FindBin(p-r*5);
      int bin1 = h->GetXaxis()->FindBin(p+r*5);
      double bglo = h->Integral(1, bin0) / bin0;
      double bghi = h->Integral(bin1, h->GetNbinsX())/(h->GetNbinsX()-bin1+1);
      double norm = h->Integral(bin0, bin1, "width");
      norm -= (bglo+bghi)*(bin1-bin0)*h->GetBinWidth(1)/2;
      double s = 0.0;
      for(int bin=bin0; bin<bin1; bin++){
	double dE = h->GetBinCenter(bin) - p;
	if(dE < 0 && h->GetBinContent(bin) > bglo)
	  s += (h->GetBinContent(bin)-bglo) * pow(dE, 2);
	else if(dE >= 0 && h->GetBinContent(bin) > bghi)
	  s += (h->GetBinContent(bin)-bghi) * pow(dE, 2);
      }
      s = sqrt(s/norm)/3;
      peak->SetRange(h->GetXaxis()->GetBinLowEdge(1),
		     h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
      peak->SetAmplitude(norm);
      peak->SetGausPars(p, s);
      peak->SetStepHeight((bglo-bghi)/norm);
      peak->SetParameter(GATPeakShape::kFt, 0.2);
      peak->SetParameter(GATPeakShape::kTau, 0.001*peaks[i]);
      peak->SetParameter(GATPeakShape::kFht, 0.01);
      peak->SetParameter(GATPeakShape::kTauHT, 0.001*peaks[i]);
      peak->SetParameter(GATPeakShape::kB, bglo);
      peak->SetParameter(GATPeakShape::kM, 0.0);
      peak->SetParameter(GATPeakShape::kQ, 0.0);
      ROOT::Fit::FitResult result = peak->FitTo(h, 0.125, false);
      TF1* f = (TF1*) peak->GetTF1().Clone(("fn_"+to_string(i)+"_"+
					    to_string(vf.size())).c_str());
      vf.push_back(f);
      peak->GetFWHM(detFWHM[i][j], detFWHMuncert[i][j]);
      peak->GetCentroid(detCent[i][j], detCentUncert[i][j]);
      delete peak;
    }
    fpeak.push_back(vf);
  }

  // plot the peaks and fit results on a single canvas
  vector<TCanvas*> cpeak;
  vector<double> speaks;
  speaks.insert(speaks.begin(), peaks.begin(), peaks.end());
  sort(speaks.begin(), speaks.end());
  ich = -1;
  for(auto const& pr : chan_map){
    ich ++;
    int nx = min(3, (int)peaks.size());
    int ny = ceil(((double)peaks.size())/3);
    TCanvas* c = new TCanvas(("cpeak_"+to_string(pr.first)).c_str(),
			     "", 800*nx, 600*ny);
    c->Divide(nx, ny);
    for(int i=0; i<(int)speaks.size(); i++){
      int j;
      for(j=0; j<(int)speaks.size(); j++)
	if(abs(peaks[j]-speaks[i]) < 1.0e-6) break;
      c->cd(i+1)->SetLogy(true);
      c->cd(i+1);
      hEpeak[ich][j]->SetXTitle("Energy (keV)");
      stringstream ss;
      ss << fixed << setprecision(3) << hEpeak[ich][j]->GetBinWidth(1);
      string sl = "Events / " + ss.str() + " keV";
      hEpeak[ich][j]->SetYTitle(sl.c_str());
      hEpeak[ich][j]->SetLineColor(4);
      hEpeak[ich][j]->SetMarkerColor(4);
      hEpeak[ich][j]->Draw();
      fpeak[ich][j]->SetLineWidth(1);
      fpeak[ich][j]->Draw("same");
      tdir->cd();
    }
    cpeak.push_back(c);
  }

  // plot the energy linearity
  vector<TCanvas*> celin;
  TLine* line = new TLine();
  line->SetLineWidth(1);
  line->SetLineStyle(2);
  line->SetLineColor(8);
  ich = -1;
  for(auto const& pr : chan_map){
    ich ++;
    vector<double> xuncert(peaks.size(), 0.0);
    TGraphErrors* g0 = new TGraphErrors(peaks.size(),
					&peaks.front(), &detCent[ich].front(),
					&xuncert.front(),
					&detCentUncert[ich].front());
    g0->SetTitle("");
    TH1F* h0 = new TH1F(("hcent0_"+to_string(pr.first)).c_str(), "",
			100, 0.0, 3000.0);
    h0->GetXaxis()->SetTitle("True Energy (keV)");
    h0->GetYaxis()->SetTitle("Peak Position (keV)");
    h0->GetYaxis()->SetRangeUser(0.0, 3000.0);
    h0->GetXaxis()->SetTitleOffset(1.2);
    h0->GetYaxis()->SetTitleOffset(0.5*7./3);
    g0->SetHistogram(h0);
    g0->SetLineColor(4);
    g0->SetMarkerColor(4);
    g0->SetMarkerStyle(20);
    g0->SetMarkerSize(0.8);
    TF1* f = new TF1(("flin_"+to_string(pr.first)).c_str(),
		     "[0]+[1]*x", 0.0, 3000.0);
    f->SetParameters(0.0, 1.0);
    g0->Fit(f, "QEMR+");
    vector<double> res(peaks.size(), 0.0);
    for(int i=0; i<(int)res.size(); i++)
      res[i] = detCent[ich][i] - f->Eval(peaks[i]);
    double maxres = *max_element(res.begin(), res.end());
    TGraphErrors* g1 = new TGraphErrors(peaks.size(),
					&peaks.front(), &res.front(),
					&xuncert.front(),
					&detCentUncert[ich].front());
    g1->SetTitle("");
    TH1F* h1 = new TH1F(("hcent_"+to_string(pr.first)).c_str(), "",
			100, 0.0, 3000.0);
    h1->GetXaxis()->SetTitle("True Energy (keV)");
    h1->GetYaxis()->SetTitle("Residual (keV)");
    h1->GetYaxis()->SetRangeUser(-1.1*maxres, 1.1*maxres);
    double lsize = g1->GetXaxis()->GetLabelSize();
    double tsize = g1->GetXaxis()->GetTitleSize();
    h1->GetXaxis()->SetLabelSize(lsize*7./3);
    h1->GetXaxis()->SetTitleSize(tsize*7./3);
    h1->GetYaxis()->SetLabelSize(lsize*7./3);
    h1->GetYaxis()->SetTitleSize(tsize*7./3);
    h1->GetYaxis()->SetTitleOffset(0.5);
    g1->SetHistogram(h1);
    g1->SetLineColor(4);
    g1->SetMarkerColor(4);
    g1->SetMarkerStyle(20);
    g1->SetMarkerSize(0.8);
    TCanvas* c = new TCanvas(("celin_"+to_string(pr.first)).c_str(),
			     "", 1200, 1000);
    TPad* p0 = new TPad(("pl0_"+to_string(pr.second)).c_str(), "",
			0.0, 0.3, 1.0, 1.0);
    p0->SetTopMargin(0.125);
    p0->SetBottomMargin(0.0);
    p0->SetLeftMargin(0.1);
    p0->SetRightMargin(0.027);
    p0->Draw();
    TPad* p1 = new TPad(("pl1_"+to_string(pr.second)).c_str(), "",
			0.0, 0.0, 1.0, 0.3);
    p1->SetTopMargin(0.0);
    p1->SetBottomMargin(0.25);
    p1->SetLeftMargin(0.1);
    p1->SetRightMargin(0.027);
    p1->Draw();
    p0->cd();
    h0->Draw("hist X+");
    g0->Draw("same P");
    p1->cd();
    h1->Draw("hist");
    g1->Draw("same P");
    line->DrawLine(0.0, 0.0, 3000.0, 0.0);
    c->Update();
    celin.push_back(c);
    tdir->cd();
  }

  // plot the FWHM vs energy and the residuals
  vector<TGraphErrors*> fwhmE;
  vector<TGraphErrors*> fwhmResid;
  vector<TF1*> ffwhmE;
  vector<TCanvas*> cfwhm;
  vector<double> roiFWHM;
  ich=-1;
  for(auto const& pr : chan_map){
    ich ++;
    vector<double> xuncert(peaks.size(), 0.0);
    TGraphErrors* g = new TGraphErrors(peaks.size(),
				       &peaks.front(), &detFWHM[ich].front(),
				       &xuncert.front(),
				       &detFWHMuncert[ich].front());
    g->SetTitle("");
    TH1F* h0 = new TH1F(("htmp0_"+to_string(pr.first)).c_str(), "",
			100, 0.0, 3000.0);
    h0->GetXaxis()->SetTitle("Energy (keV)");
    h0->GetYaxis()->SetTitle("FWHM (keV)");
    h0->GetYaxis()->SetRangeUser(0.0, 4.0);
    h0->GetXaxis()->SetTitleOffset(1.2);
    h0->GetYaxis()->SetTitleOffset(0.5*7./3);
    g->SetHistogram(h0);
    g->SetLineColor(4);
    g->SetMarkerColor(4);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.8);
    TF1* f = new TF1(("fres_"+to_string(pr.first)).c_str(),
		     "sqrt([0]+[1]*x+[2]*x*x)", 0.0, 3000.0);
    f->SetParameters(0.38, 5e-3, 5e-6);
    g->Fit(f, "QEMR+");
    roiFWHM.push_back(f->Eval(2039.0));
    fwhmE.push_back(g);
    ffwhmE.push_back(f);
    vector<double> res(peaks.size(), 0.0);
    for(int i=0; i<(int)res.size(); i++)
      res[i] = detFWHM[ich][i] - f->Eval(peaks[i]);
    double maxres = *max_element(res.begin(), res.end());
    g = new TGraphErrors(peaks.size(), &peaks.front(), &res.front(),
			 &xuncert.front(), &detFWHMuncert[ich].front());
    g->SetTitle("");
    TH1F* h1 = new TH1F(("htmp1_"+to_string(pr.first)).c_str(), "",
			100, 0.0, 3000.0);
    h1->GetXaxis()->SetTitle("Energy (keV)");
    h1->GetYaxis()->SetTitle("Residual (keV)");
    h1->GetYaxis()->SetRangeUser(-1.1*maxres, 1.1*maxres);
    double lsize = g->GetXaxis()->GetLabelSize();
    double tsize = g->GetXaxis()->GetTitleSize();
    h1->GetXaxis()->SetLabelSize(lsize*7./3);
    h1->GetXaxis()->SetTitleSize(tsize*7./3);
    h1->GetYaxis()->SetLabelSize(lsize*7./3);
    h1->GetYaxis()->SetTitleSize(tsize*7./3);
    h1->GetYaxis()->SetTitleOffset(0.5);
    g->SetHistogram(h1);
    g->SetLineColor(4);
    g->SetMarkerColor(4);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.8);
    fwhmResid.push_back(g);
    TCanvas* c = new TCanvas(("cfwhm_"+to_string(pr.first)).c_str(),
			     "", 1200, 1000);
    TPad* p0 = new TPad(("p0_"+to_string(pr.second)).c_str(), "",
			0.0, 0.3, 1.0, 1.0);
    p0->SetTopMargin(0.125);
    p0->SetBottomMargin(0.0);
    p0->SetLeftMargin(0.1);
    p0->SetRightMargin(0.027);
    p0->Draw();
    TPad* p1 = new TPad(("p1_"+to_string(pr.second)).c_str(), "",
			0.0, 0.0, 1.0, 0.3);
    p1->SetTopMargin(0.0);
    p1->SetBottomMargin(0.25);
    p1->SetLeftMargin(0.1);
    p1->SetRightMargin(0.027);
    p1->Draw();
    p0->cd();
    h0->Draw("hist X+");
    fwhmE[ich]->Draw("same P");
    line->DrawLine(2039.0, 0.0, 2039.0, roiFWHM[ich]);
    line->DrawLine(0.0, roiFWHM[ich], 2039.0, roiFWHM[ich]);
    p1->cd();
    h1->Draw("hist");
    fwhmResid[ich]->Draw("same P");
    line->DrawLine(0.0, 0.0, 3000.0, 0.0);
    c->Update();
    cfwhm.push_back(c);
    tdir->cd();
  }
  delete line;

  // determine the PSA efficiencies
  if(do_psa){
    int i = -1;
    for(auto const& pr : chan_map){
      i ++;
      double dep = 2615-2*511.;
      double dres = ffwhmE[i]->Eval(dep);
      pair<double, double> d0(max(1525.0, dep-8*dres-50), dep-8*dres);
      pair<double, double> d1(dep-3*dres, dep+3*dres);
      pair<double, double> d2(1650.0, 1700.0);
      double sep = 2615-2*511.;
      double sres = ffwhmE[i]->Eval(sep);
      pair<double, double> s0(sep-8*sres-50, sep-8*sres);
      pair<double, double> s1(sep-3*sres, sep-3*sres);
      pair<double, double> s2(sep+8*sres, sep+8*sres+50);
      TAxis* xa = havseE[i]->GetXaxis();
      TH1D* h0 = havseE[i]->ProjectionY("havse0_tmp",
					xa->FindBin(d0.first),
					xa->FindBin(d0.second));
      TH1D* h1 = havseE[i]->ProjectionY("havse1_tmp",
					xa->FindBin(d1.first),
					xa->FindBin(d1.second));
      TH1D* h2 = havseE[i]->ProjectionY("havse2_tmp",
					xa->FindBin(d2.first),
					xa->FindBin(d2.second));
      for(int j=1; j<=h1->GetNbinsX(); j++){
	double b0 = h0->GetBinContent(j) / (d0.second-d0.first);
	double b2 = h2->GetBinContent(j) / (d2.second-d2.first);
	double b1 = h1->GetBinContent(j)-(b0+b2)*0.5*(d1.second-d1.first);
	h1->SetBinContent(j, max(0.0, b1));
      }
      double dacc = h1->Integral(h1->GetXaxis()->FindBin(-1.0),
				 h1->GetNbinsX());
      dacc /= h1->Integral(0, h1->GetNbinsX());

      delete h0;
      delete h1;
      delete h2;
      h0 = havseE[i]->ProjectionY("havse0_tmp",
				  xa->FindBin(s0.first),
				  xa->FindBin(s0.second));
      h1 = havseE[i]->ProjectionY("havse1_tmp",
				  xa->FindBin(s1.first),
				  xa->FindBin(s1.second));
      h2 = havseE[i]->ProjectionY("havse2_tmp",
				  xa->FindBin(s2.first),
				  xa->FindBin(s2.second));
      for(int j=1; j<=h1->GetNbinsX(); j++){
	double b0 = h0->GetBinContent(j) / (s0.second-s0.first);
	double b2 = h2->GetBinContent(j) / (s2.second-s2.first);
	double b1 = h1->GetBinContent(j)-(b0+b2)*0.5*(s1.second-s1.first);
	h1->SetBinContent(j, max(0.0, b1));
      }
      double sacc = h1->Integral(h1->GetXaxis()->FindBin(-1.0),
				 h1->GetNbinsX());
      sacc /= h1->Integral(0, h1->GetNbinsX());
      delete h0;
      delete h1;
      delete h2;
      cout << "Channel " << pr.first << ":" << endl;
      cout << "  AvsE DEP acceptance: " << dacc << endl;
      cout << "  AvsE SEP acceptance: " << sacc << endl;
      stringstream ssavse;
      ssavse << "DEP: " << fixed << setprecision(1) << 100*dacc << "%  "
	     << "SEP: " << fixed << setprecision(1) << 100*sacc << "%";
      havseE[i]->SetTitle(("#font[132]{"+ssavse.str()+"}").c_str());
      h0 = haoeE[i]->ProjectionY("haoe0_tmp",
				 xa->FindBin(d0.first),
				 xa->FindBin(d0.second));
      h1 = haoeE[i]->ProjectionY("haoe1_tmp",
				 xa->FindBin(d1.first),
				 xa->FindBin(d1.second));
      h2 = haoeE[i]->ProjectionY("haoe2_tmp",
				 xa->FindBin(d2.first),
				 xa->FindBin(d2.second));
      for(int j=1; j<=h1->GetNbinsX(); j++){
	double b0 = h0->GetBinContent(j) / (d0.second-d0.first);
	double b2 = h2->GetBinContent(j) / (d2.second-d2.first);
	double b1 = h1->GetBinContent(j)-(b0+b2)*0.5*(d1.second-d1.first);
	h1->SetBinContent(j, max(0.0, b1));
      }
      dacc = h1->Integral(h1->GetXaxis()->FindBin(-1.0),
				 h1->GetNbinsX());
      dacc /= h1->Integral(0, h1->GetNbinsX());
      delete h0;
      delete h1;
      delete h2;
      h0 = haoeE[i]->ProjectionY("haoe0_tmp",
				 xa->FindBin(s0.first),
				 xa->FindBin(s0.second));
      h1 = haoeE[i]->ProjectionY("haoe1_tmp",
				 xa->FindBin(s1.first),
				 xa->FindBin(s1.second));
      h2 = haoeE[i]->ProjectionY("haoe2_tmp",
				 xa->FindBin(s2.first),
				 xa->FindBin(s2.second));
      for(int j=1; j<=h1->GetNbinsX(); j++){
	double b0 = h0->GetBinContent(j) / (s0.second-s0.first);
	double b2 = h2->GetBinContent(j) / (s2.second-s2.first);
	double b1 = h1->GetBinContent(j)-(b0+b2)*0.5*(s1.second-s1.first);
	h1->SetBinContent(j, max(0.0, b1));
      }
      sacc = h1->Integral(h1->GetXaxis()->FindBin(-1.0),
				 h1->GetNbinsX());
      sacc /= h1->Integral(0, h1->GetNbinsX());
      delete h0;
      delete h1;
      delete h2;
      cout << "  A/E DEP acceptance: " << dacc << endl;
      cout << "  A/E SEP acceptance: " << sacc << endl;
      stringstream ssaoe;
      ssaoe << "DEP: " << fixed << setprecision(1) << 100*dacc << "%  "
	    << "SEP: " << fixed << setprecision(1) << 100*sacc << "%";
      haoeE[i]->SetTitle(("#font[132]{"+ssaoe.str()+"}").c_str());
    }
  }

  // plot the PSA parameters
  vector<TCanvas*> cpsa;
  vector<TCanvas*> cpsaE;
  ich = -1;
  for(auto const& pr : chan_map){
    ich ++;
    TCanvas* c = new TCanvas(("cpsa_"+to_string(pr.first)).c_str(), "",
			     0, 0, 1400, 400);
    c->Divide(3, 1);
    c->cd(1)->SetLogy(true);
    c->cd(1);
    havse2615[ich]->Draw();
    c->cd(2)->SetLogy(true);
    c->cd(2);
    haoe2615[ich]->Draw();
    c->cd(3)->SetLogy(true);
    c->cd(3);
    hdcr2615[ich]->Draw();
    cpsa.push_back(c);
    c = new TCanvas(("cpsaE_"+to_string(pr.first)).c_str(), "",
		    0, 0, 1400, 400);
    c->Divide(3, 1);
    c->cd(1)->SetLogz(true);
    c->cd(1);
    havseE[ich]->Draw("colz");
    c->cd(2)->SetLogz(true);
    c->cd(2);
    haoeE[ich]->Draw("colz");
    c->cd(3)->SetLogz(true);
    c->cd(3);
    hdcrE[ich]->Draw("colz");
    cpsaE.push_back(c);
  }

  if(write_dir != ""){
    ich = -1;
    for(auto const& pr : chan_map){
      ich ++;
      string s = "_" + det_serial[pr.second] + "_";
      for(int i=0; i<(int)sources.size(); i++){
	if(i<(int)sources.size()-1) s += sources[i] + "_";
	else s += sources[i] + ".pdf";
      }
      if(cenergy[ich]) cenergy[ich]->Print((write_dir+"/energy"+s).c_str());
      if(cpeak[ich])     cpeak[ich]->Print((write_dir+"/peakfit"+s).c_str());
      if(celin[ich])     celin[ich]->Print((write_dir+"/elin"+s).c_str());
      if(cfwhm[ich])     cfwhm[ich]->Print((write_dir+"/fwhm"+s).c_str());
      if(cpsa[ich])       cpsa[ich]->Print((write_dir+"/psa"+s).c_str());
      if(cpsaE[ich])     cpsaE[ich]->Print((write_dir+"/psaE"+s).c_str());
    }
  }
    
  // write histograms to output file
  TFile* outfile = new TFile(outfname.c_str(), "recreate");
  outfile->cd();
  for(int i=0; i<(int)chan_map.size(); i++){
    if(hECal[i])             hECal[i]->Write();
    if(hEFCal[i])           hEFCal[i]->Write();
    if(hEFCCal[i])         hEFCCal[i]->Write();
    if(hEFCal_avse[i]) hEFCal_avse[i]->Write();
    if(hEFCal_aoe[i])   hEFCal_aoe[i]->Write();
    if(hdcr2615[i])       hdcr2615[i]->Write();
    if(havse2615[i])     havse2615[i]->Write();
    if(haoe2615[i])       haoe2615[i]->Write();
    if(hdcrE[i])             hdcrE[i]->Write();
    if(havseE[i])           havseE[i]->Write();
    if(haoeE[i])             haoeE[i]->Write();
    if(cenergy[i])         cenergy[i]->Write();
    if(cpeak[i])             cpeak[i]->Write();
    if(celin[i])             celin[i]->Write();
    if(cfwhm[i])             cfwhm[i]->Write();
    if(cpsa[i])               cpsa[i]->Write();
    if(cpsaE[i])             cpsaE[i]->Write();
  }
  outfile->Close();
  
  return 0;
}
