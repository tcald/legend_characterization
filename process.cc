#include "utils.hh"
#include <MGTEvent.hh>
#include <MGTWaveform.hh>
#include <MGWFPoleZeroCorrection.hh>
#include <MGWFTrapezoidalFilter.hh>
#include <MGWFAsymTrapezoidalFilter.hh>
#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLine.h>
#include <TPolyMarker.h>
#include <json/json.h>
#include <json/value.h>
#include <json/reader.h>
#include <utility>
#include <numeric>
#include <functional>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <assert.h>

using namespace std;

TDirectory* tdir;

int main(int argc, char* argv[]){

  // default parameters, later set in json configuration
  double slow_ramp =     2500.0;
  double slow_flat =     1500.0;
  double fast_ramp =     40.0;
  double fast_flat =     100.0;
  double fast_fall =     2000.0;
  double avse_ramp =     100.0;
  double avse_flat =     300.0;
  double t0_thresh =     2.0;
  double pz_decay  =     49000.0;
  double pzd_decay =     0.0;
  double os_amplitude =  0.0;
  double os_decay     =  0.0;
  int    nbase_samples = 200;
  int    nefit_samples = 300;
  int    ndcr_samples  = 100;

  // root setup
  gStyle->SetOptStat(0);
  gStyle->SetLabelFont(132, "XYZ");
  gStyle->SetTitleFont(132, "XYZ");
  tdir = gROOT->CurrentDirectory();
  TPolyMarker* marker = new TPolyMarker(8);
  marker->SetMarkerSize(0.8);
  marker->SetMarkerStyle(4);
  marker->SetMarkerColor(2);

  // option handling
  map<int, int> chan_map;
  string base_dir = "";
  int max_wf = 0;
  string outfname = "";
  string fnamebase = "";
  int run_start = -1;
  int run_end = -1;
  string infname = "";
  bool fit_tail = true;
  bool json_config = false;
  bool float_pz = true;
  int nbits = 16;
  int write_wf = 0;
  double rc_decay = 0.0;
  int update_percentage = 5;
  static struct option opts[]{
    {"help",             no_argument, NULL, 'h'},
    {"datadir",    required_argument, NULL, 'd'},
    {"channel",    required_argument, NULL, 'c'},
    {"nwaveforms", required_argument, NULL, 'n'},
    {"outfile",    required_argument, NULL, 'o'},
    {"fnamebase",  required_argument, NULL, 'f'},
    {"runstart",   required_argument, NULL, 'r'},
    {"runend",     required_argument, NULL, 'R'},
    {"infile",     required_argument, NULL, 'i'},
    {"fittail",          no_argument, NULL, 'F'},
    {"bits",       required_argument, NULL, 'b'},
    {"writewf",    required_argument, NULL, 'w'},
    {"decayconst", required_argument, NULL, 'D'},
    {"jsoncofig",  required_argument, NULL, 'j'},
    {"update",     required_argument, NULL, 'u'}
  };
  int opt = getopt_long(argc, argv,
			"hd:c:n:o:f:r:R:i:Fb:w:D:j:u:", opts, NULL);
  while(opt != -1){
    switch(opt){
    case 'h':
      cout << "options:"                             << endl;
      cout << "  -i input filename"                  << endl;
      cout << "  -o output filename"                 << endl;
      cout << "  -r starting run number"             << endl;
      cout << "  -R ending run number"               << endl;
      cout << "  -d base directory to data"          << endl;
      cout << "  -f base of filenames"               << endl;
      cout << "  -c channel to be analyzed"          << endl;
      cout << "  -n number of events to analyze"     << endl;
      cout << "  -F do not fit exponential tail"     << endl;
      cout << "  -b number of ADC bits"              << endl;
      cout << "  -w plot every n'th wf"              << endl;
      cout << "  -D decay constant in ns"            << endl;
      cout << "  -j name of json configuration file" << endl;
      cout << "  -u update interval (percentage)"    << endl;
      return 0;
    case 'd': base_dir  = string(optarg); break;
    case 'c': chan_map[atoi(optarg)] = (int) chan_map.size()-1; break;
    case 'n': max_wf    = atoi(optarg);   break;
    case 'o': outfname  = string(optarg); break;
    case 'f': fnamebase = string(optarg); break;
    case 'r': run_start = atoi(optarg);   break;
    case 'R': run_end   = atoi(optarg);   break;
    case 'i': infname   = string(optarg); break;
    case 'F': fit_tail  = false;          break;
    case 'b': nbits     = atoi(optarg);   break;
    case 'w': write_wf  = atoi(optarg);   break;
    case 'D': rc_decay  = atof(optarg);   break;
    case 'j':{
      ifstream jfile(optarg);
      Json::CharReaderBuilder reader;
      Json::Value value;
      string errors;
      if(!parseFromStream(reader, jfile, &value, &errors)){
	cout << errors << endl;
	return 2;
      }
      cout << "reading parameters for " << value["manufacturer"].asString()
	   << " " << value["serial"].asString() << endl;
      SetJson(value, "slow_ramp",     slow_ramp);
      SetJson(value, "slow_flat",     slow_flat);
      SetJson(value, "fast_ramp",     fast_ramp);
      SetJson(value, "fast_flat",     fast_flat);
      SetJson(value, "fast_fall",     fast_fall);
      SetJson(value, "avse_ramp",     avse_ramp);
      SetJson(value, "t0_thresh",     t0_thresh);
      SetJson(value, "pzd_decay",     pzd_decay);
      SetJson(value, "os_amplitude",  os_amplitude);
      SetJson(value, "nbase_samples", nbase_samples);
      SetJson(value, "nefit_samples", nefit_samples);
      SetJson(value, "ndcr_samples",  ndcr_samples);
      if(SetJson(value, "pz_decay", pz_decay)) float_pz = false;
      json_config = true;
      break;
    }
    case 'u': update_percentage = atoi(optarg); break;
    default: return 1;
    }
    opt = getopt_long(argc, argv,
		      "hd:c:n:o:f:r:R:i:Fb:w:D:j:u:", opts, NULL);
  }
  assert(infname != "" || (run_start > 0 && run_end > 0));
  if(rc_decay != 0.0 && json_config){
    pz_decay = rc_decay;
    float_pz = false;
    cout << "over-riding default pz constant with user-specified value "
	 << pz_decay << " ns" << endl;
  }
    
  // output tree
  vector<int>    channel, maxtime, mintime, trapmaxtime;
  vector<double> baseline, baserms, maxval, minval, trappick;
  vector<double> t0, t1, t10, t50, t90, t99, imax, dcrslope;
  vector<double> trapmax, time, deltat;
  vector<vector<double> > exp_param;
  TFile* outfile = new TFile(outfname.c_str(), "recreate");
  if(write_wf)
    for(auto const& p : chan_map)
      outfile->mkdir(("ch"+to_string(p.first)).c_str());
  tdir->cd();
  TTree* outtree = new TTree("tree", "tree");
  outtree->Branch("channel", &channel);
  outtree->Branch("maxtime", &maxtime);
  outtree->Branch("mintime", &mintime);
  outtree->Branch("trapmaxtime", &trapmaxtime);
  outtree->Branch("baseline", &baseline);
  outtree->Branch("baserms", &baserms);
  outtree->Branch("maxval", &maxval);
  outtree->Branch("minval", &minval);
  outtree->Branch("t0", &t0);
  outtree->Branch("t1", &t1);
  outtree->Branch("t10", &t10);
  outtree->Branch("t50", &t50);
  outtree->Branch("t90", &t90);
  outtree->Branch("t99", &t99);
  outtree->Branch("imax", &imax);
  outtree->Branch("dcrslope", &dcrslope);
  outtree->Branch("trapmax", &trapmax);
  outtree->Branch("trappick", &trappick);
  outtree->Branch("time", &time);
  outtree->Branch("deltat", &deltat);
  outtree->Branch("exp_param", &exp_param);
  outtree->SetDirectory(outfile);

  // output histograms
  vector<double> tlast(chan_map.size(), 0.0);
  vector<int> wf_count(chan_map.size(), -1);
  vector<double> xfit, yfit;
  vector<TH1D*> hdeltat(chan_map.size(), NULL);
  vector<TH1D*> henergy(chan_map.size(), NULL);
  vector<TH1D*> henergyf(chan_map.size(), NULL);
  vector<TH2D*> hbase_energy(chan_map.size(), NULL);
  vector<TH2D*> hbrms_energy(chan_map.size(), NULL);
  vector<TH2D*> hdecay_energy(chan_map.size(), NULL);
  vector<TH2D*> hdecay_deltat(chan_map.size(), NULL);
  vector<TH2D*> hamp_energy(chan_map.size(), NULL);
  vector<TH2D*> hdcr_energy(chan_map.size(), NULL);
  vector<TH2D*> hrise_energy(chan_map.size(), NULL);
  int adcbins = (1<<nbits) / 16;
  for(int i=0; i<(int)hbase_energy.size(); i++){
    hdeltat[i] = new TH1D(("hdeltat_"+to_string(chan_map[i])).c_str(), "",
			  1e5, 0.0, 10000.0);
    hdeltat[i]->SetXTitle("Delta t (#mu s)");
    hdeltat[i]->SetYTitle("Entries");
    henergy[i] = new TH1D(("henergy_"+to_string(chan_map[i])).c_str(), "",
			  1<<nbits, 0.0, 1<<nbits);
    henergy[i]->SetXTitle("Trap Maximum (ADC)");
    henergy[i]->SetYTitle("Entries");
    henergyf[i] = new TH1D(("henergyf_"+to_string(chan_map[i])).c_str(), "",
			   1<<nbits, 0.0, 1<<nbits);
    henergyf[i]->SetXTitle("Fixed Time Pickoff (ADC)");
    henergyf[i]->SetYTitle("Entries");
    hbase_energy[i] = new TH2D(("hbase_energy_"+to_string(chan_map[i])).c_str(), "",
			       adcbins, 0.0, 10000.0, 400, -20.0, 20.0);
    hbase_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hbase_energy[i]->SetYTitle("Baseline (ADC)");
    hbrms_energy[i] = new TH2D(("hbrms_energy_"+to_string(chan_map[i])).c_str(), "",
			       adcbins, 0.0, 10000.0, 400, 0.0, 40.0);
    hbrms_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hbrms_energy[i]->SetYTitle("Baseline RMS (ADC)");
    if(fit_tail){
      hdecay_energy[i] = new TH2D(("hdecay_energy_"+to_string(chan_map[i])).c_str(),
				  "", adcbins, 0.0, 10000.0,
				  400, 0.0, 100.0);
      hdecay_energy[i]->SetXTitle("Trap Maximum (ADC)");
      hdecay_energy[i]->SetYTitle("Decay Constant (#mu s)");
      hdecay_deltat[i] = new TH2D(("hdecay_deltat_"+to_string(chan_map[i])).c_str(),
				  "", 1000, 0.0, 1000.0, 400, 0.0, 100.0);
      hdecay_deltat[i]->SetYTitle("Decay Constant (#mu s)");
      hdecay_deltat[i]->SetXTitle("#Delta t (#mu s)");
    }
    hamp_energy[i] = new TH2D(("hamp_energy_"+to_string(chan_map[i])).c_str(), "",
			      adcbins, 0.0, 10000.0, 500, 0.0, 5000.0);
    hamp_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hamp_energy[i]->SetYTitle("Maximum Current");
    hdcr_energy[i] = new TH2D(("hdcr_energy_"+to_string(chan_map[i])).c_str(), "",
			      adcbins, 0.0, 10000.0,
			      400, -100.0, 100.0);
    hdcr_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hdcr_energy[i]->SetYTitle("Raw DCR Slope");
    hrise_energy[i] = new TH2D(("hrise_energy_"+to_string(chan_map[i])).c_str(), "",
			       2*adcbins, 0.0, 10000.0,
			       1000, 0.0, 5000.0);
    hrise_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hrise_energy[i]->SetYTitle("Rise Time (ns)");
  }

  // input reader
  TChain tree("MGTree", "MGTree");                                  
  string fname = base_dir + "/" + fnamebase;                
  if(infname != "") tree.Add(infname.c_str());
  else
    for(int irun=run_start; irun<=run_end; irun++)
      tree.Add((fname + to_string(irun) + ".root").c_str());
  int nentries = (int) tree.GetEntries();
  if(max_wf > 0) nentries = min(nentries, max_wf);
  TTreeReader reader(&tree);
  TTreeReaderValue<MGTEvent> event(reader, "event");

  MGWFPoleZeroCorrection* pole_zero=new MGWFPoleZeroCorrection();
  MGWFTrapezoidalFilter* slow_trap=new MGWFTrapezoidalFilter(slow_ramp,slow_flat);
  MGWFAsymTrapezoidalFilter* fast_trap = new MGWFAsymTrapezoidalFilter();
  fast_trap->SetRampTime(fast_ramp);
  fast_trap->SetFlatTime(fast_flat);
  fast_trap->SetFallTime(fast_fall);
  MGWFTrapezoidalFilter* avse_trap=new MGWFTrapezoidalFilter(avse_ramp,avse_flat);
  
  int iev = -1;
  int cpercent = -1;
  int lpercent = -1;
  while(reader.Next()){
    iev ++;
    cpercent = (int) (100.0*iev/nentries);
    if(cpercent % update_percentage == 0 && cpercent != lpercent)
      cout << cpercent << "%" << endl;
    lpercent = cpercent;
    if(iev >= nentries) break;
    if(iev % 100000 == 0 && iev > 0) outtree->AutoSave();
    TClonesArray* wfs = event->GetWaveforms();
    int nwf = (int) chan_map.size();
    // clear previous values
    channel.assign(nwf, 0);
    maxtime.assign(nwf, 0);
    mintime.assign(nwf, 0);
    trapmaxtime.assign(nwf, 0);
    baseline.assign(nwf, 0.0);
    baserms.assign(nwf, 0.0);
    maxval.assign(nwf, 0.0);
    minval.assign(nwf, 0.0);
    t0.assign(nwf, 0.0);
    t1.assign(nwf, 0.0);
    t10.assign(nwf, 0.0);
    t50.assign(nwf, 0.0);
    t90.assign(nwf, 0.0);
    t99.assign(nwf, 0.0);
    imax.assign(nwf, 0.0);
    dcrslope.assign(nwf, 0.0);
    trapmax.assign(nwf, 0.0);
    time.assign(nwf, 0.0);
    deltat.assign(nwf, 0.0);
    trappick.assign(nwf, 0.0);
    exp_param.assign(nwf, vector<double>(4, 0.0));
    // start analyzing waveforms for channels selected by user
    for(int iwf=0; iwf<(int)wfs->GetEntriesFast(); iwf++){
      if(chan_map.find(iwf) == chan_map.end()) continue;
      int index = chan_map[iwf];
      wf_count[index] ++;
      bool write = false;
      if(write_wf>0) if(wf_count[index] % write_wf == 0) write = true;
      channel[index] = iwf;
      MGTWaveform* wf = (MGTWaveform*) wfs->At(iwf);
      TH1D* hwf_orig = NULL;
      if(write) hwf_orig = wf->GimmeUniqueHist();
      time[index] = event->GetTime() + wf->GetTOffset();
      deltat[index] = time[index] - tlast[index];
      hdeltat[index]->Fill(deltat[index] * 1e6);
      double sampling = wf->GetSamplingPeriod();
      vector<double> vwf = wf->GetVectorData();
      assert((int)vwf.size()>nbase_samples && (int)vwf.size()>=2*nefit_samples);
      baseline[index] = accumulate(vwf.begin(), vwf.begin()+nbase_samples, 0);
      baseline[index] /= nbase_samples;
      double base_orig = baseline[index];
      // exponential baseline correction from last event, fit exponential tail
      if(fit_tail){
	for_each(xfit.begin(), xfit.end(), [&](double& s){s-=deltat[index]*1e9;});
	for(int i=0; i<nefit_samples; i++){
	  xfit.push_back(i*sampling);
	  yfit.push_back(vwf[i]);
	}
	TF1* fn = new TF1("ftmp", "[0]+[1]*exp(-x/[2])",
			  xfit[0], xfit[xfit.size()-1]);
	fn->SetParameter(0, baseline[index]);
	fn->SetParameter(1, vwf[0]-baseline[index]);
	if(float_pz) fn->SetParameter(2, pz_decay);
	else fn->FixParameter(2, pz_decay);
	TGraph* gr = new TGraph(xfit.size(), &xfit.front(), &yfit.front());
	gr->Fit(fn, "QRNW");
	for(int i=0; i<3; i++) exp_param[index][i] = fn->GetParameter(i);
	exp_param[index][3] = fn->GetChisquare() / (xfit.size()-3);
	xfit.resize(nefit_samples, 0.0);
	yfit.resize(nefit_samples, 0.0);
	for(int i=(int)vwf.size()-nefit_samples; i<(int)vwf.size(); i++){
	  int j = i - (int)vwf.size() + nefit_samples;
	  xfit[j] = i*sampling;
	  yfit[j] = vwf[i];
	}
	for(int i=0; i<(int)vwf.size(); i++)
	  vwf[i] = vwf[i] - fn->Eval(i*sampling);
	delete gr;
	delete fn;
	baseline[index] = accumulate(vwf.begin(), vwf.begin()+nbase_samples,0);
	baseline[index] /= nbase_samples;
      }
      // offset baseline subtraction, baseline rms, extrema, and decay constant
      tlast[index] = time[index];
      for_each(vwf.begin(), vwf.end(), [&](double& s){s-=baseline[index];});
      baserms[index]=sqrt(inner_product(vwf.begin(), vwf.begin()+nbase_samples,
					vwf.begin(), 0.0) / nbase_samples);
      wf->SetData(vwf);
      TH1D* hwf_cor = NULL;
      if(write) hwf_cor = wf->GimmeUniqueHist();
      vector<double>::iterator itmax = max_element(vwf.begin(), vwf.end());
      vector<double>::iterator itmin = min_element(vwf.begin(), vwf.end());
      maxtime[index] = distance(vwf.begin(), itmax);
      mintime[index] = distance(vwf.begin(), itmin);
      maxval[index]  = *itmax;
      minval[index]  = *itmin;
      if(itmax == vwf.begin() || itmax == vwf.end()) continue;
      // second pass pz correction and adc overshoot correction
      MGTWaveform* pz = new MGTWaveform();
      pole_zero->SetDecayConstant(exp_param[index][2]);
      pole_zero->TransformOutOfPlace(*wf, *pz);
      vector<double> vwff = pz->GetVectorData();
      MGTWaveform* wff = pz;
      MGTWaveform* pzd = NULL;
      if(pzd_decay > 0.0){
	pzd = new MGTWaveform();
	PZDiff(*pz, *pzd, pzd_decay);
	vwff = pzd->GetVectorData();
	wff = pzd;
      }
      MGTWaveform* os = NULL;
      if(os_decay > 0.0 && os_amplitude > 0.0){
	os = new MGTWaveform();
	OvershootCorrection(*pzd, *os, os_amplitude, os_decay);
	vwff = os->GetVectorData();
	wff = os;
      }
      // get time points
      MGTWaveform* ft = new MGTWaveform();
      ft->SetTOffset(fast_fall+fast_flat);
      MGTWaveform* st = new MGTWaveform();
      MGTWaveform* at = new MGTWaveform();
      fast_trap->TransformOutOfPlace(*wff, *ft);
      slow_trap->TransformOutOfPlace(*wff, *st);
      avse_trap->TransformOutOfPlace(*wff, *at);
      vector<double> vft = ft->GetVectorData();
      vector<double> vst = st->GetVectorData();
      vector<double> vat = at->GetVectorData();
      vector<double>::iterator itstmax = max_element(vst.begin(), vst.end());
      vector<double>::iterator itftmax = max_element(vft.begin(), vft.end()-1);
      imax[index] = *max_element(vat.begin(), vat.end());
      trapmaxtime[index] = distance(vst.begin(), itstmax);
      trapmax[index] = *itstmax;
      hbase_energy[index]->Fill(trapmax[index], baseline[index]);
      hbrms_energy[index]->Fill(trapmax[index], baserms[index]);
      if(fit_tail){
	hdecay_energy[index]->Fill(trapmax[index], exp_param[index][2]/1000);
	hdecay_deltat[index]->Fill(deltat[index]*1e6,exp_param[index][2]/1000);
      }
      double ftbase = accumulate(vft.begin(), vft.begin()+nbase_samples, 0.0);
      ftbase /= nbase_samples;
      for_each(vft.begin(), vft.end(), [&](double& s){s-=ftbase;});
      double ftrms = sqrt(inner_product(vft.begin(), vft.begin()+nbase_samples,
					vft.begin(), 0.0) / nbase_samples);
      double ftoffset = ft->GetTOffset()/sampling;
      for(int i=distance(vft.begin(), itftmax); i>=0; i--){
	double val = vft[i] / (*itftmax);
	double diff = vft[i+1] - vft[i];
	if(diff == 0.0) diff = 1.0e12;
	if(val < 0.01 && t1[index] == 0.0)
	  t1[index] = i + (vft[i+1]-0.01*(*itftmax))/diff + ftoffset;
	if(val < 0.1 && t10[index] == 0.0)
	  t10[index] = i + (vft[i+1]-0.1*(*itftmax))/diff + ftoffset;
	if(val < 0.5 && t50[index] == 0.0)
	  t50[index] = i + (vft[i+1]-0.5*(*itftmax))/diff + ftoffset;
	if(val < 0.9 && t90[index] == 0.0)
	  t90[index] = i + (vft[i+1]-0.9*(*itftmax))/diff + ftoffset;
	if(val < 0.99 && t99[index] == 0.0)
	  t99[index] = i + (vft[i+1]-0.99*(*itftmax))/diff + ftoffset;
	if(vft[i] < t0_thresh*ftrms){
	  t0[index] = i + ftoffset;
	  break;
	}
      }
      t0[index] = min(max(0., t0[index]), (double)vft.size());
      t1[index] = min(max(t0[index],   t1[index]), (double)vft.size());
      t10[index]= min(max(t1[index],  t10[index]), (double)vft.size());
      t50[index]= min(max(t10[index], t50[index]), (double)vft.size());
      t90[index]= min(max(t50[index], t90[index]), (double)vft.size());
      t99[index]= min(max(t90[index], t99[index]), (double)vft.size());
      // fixed time energy pickoff
      int pickoff = t0[index] + (int) (slow_ramp+slow_flat-400)/sampling;
      if(pickoff < 0 || pickoff >= (int) vst.size()) trappick[index] = 0;
      else trappick[index] = vst[pickoff];
      henergy[index]->Fill(trapmax[index]);
      henergyf[index]->Fill(trappick[index]);
      hamp_energy[index]->Fill(trappick[index], imax[index]);
      hrise_energy[index]->Fill(trappick[index],
				(t99[index]-t1[index])*sampling-fast_ramp);
      // compute dcr
      int dcr_start = min((int)t99[index], (int)vwff.size()-2*ndcr_samples-1);
      dcrslope[index]  = accumulate(vwff.end()-ndcr_samples, vwff.end(), 0.0);
      dcrslope[index] -= accumulate(vwff.begin()+dcr_start,
				    vwff.begin()+dcr_start+ndcr_samples, 0.0);
      dcrslope[index] /= vwff.size() - dcr_start - ndcr_samples;
      hdcr_energy[index]->Fill(trappick[index], dcrslope[index]);
      // plot waveform and write to output file
      if(write){
	for(int i=1; i<=(int)hwf_orig->GetNbinsX(); i++)
	  hwf_orig->SetBinContent(i, hwf_orig->GetBinContent(i) - base_orig);
	TH1D* hpz = pz->GimmeUniqueHist();
	TH1D* hpzd= NULL;
	if(pzd) hpzd = pzd->GimmeUniqueHist();
	TH1D* hos = NULL;
	if(os) hos = os->GimmeUniqueHist();
	TH1D* hft = ft->GimmeUniqueHist();
	TH1D* hst = st->GimmeUniqueHist();
	TH1D* hat = at->GimmeUniqueHist();
	string cname = "c_"+to_string(iev)+"_"+to_string(chan_map[iwf]);
	TCanvas* c = new TCanvas(cname.c_str(), "", 0, 0, 1200, 800);
	hwf_orig->SetTitle("");
	hwf_orig->SetLineColor(2);
	hwf_orig->SetMarkerColor(2);
	hwf_orig->Draw();
	hwf_cor->SetLineColor(4);
	hwf_cor->SetMarkerColor(4);
	hwf_cor->Draw("same");
	hpz->SetLineColor(28);
	hpz->SetMarkerColor(28);
	hpz->Draw("same");
	if(hpzd){
	  hpzd->SetLineColor(40);
	  hpzd->SetMarkerColor(40);
	  hpzd->Draw("same");
	}
	if(hos){
	  hos->SetLineColor(9);
	  hos->SetMarkerColor(9);
	  hos->Draw("same");
	}
	hft->SetLineColor(6);
	hft->SetMarkerColor(6);
	hft->Draw("same c");
	hst->SetLineColor(8);
	hst->SetMarkerColor(8);
	hst->Draw("same");
	hat->SetLineColor(7);
	hat->SetMarkerColor(7);
	hat->Draw("same");
	double fto = ftoffset;
	marker->SetPoint(0, t0[index]*sampling, ft->At((int)(t0[index]-fto)));
	marker->SetPoint(1, t1[index]*sampling, ft->At((int)(t1[index]-fto)));
	marker->SetPoint(2,t10[index]*sampling, ft->At((int)(t10[index]-fto)));
	marker->SetPoint(3,t50[index]*sampling, ft->At((int)(t50[index]-fto)));
	marker->SetPoint(4,t90[index]*sampling, ft->At((int)(t90[index]-fto)));
	marker->SetPoint(5,t99[index]*sampling, ft->At((int)(t99[index]-fto)));
	marker->SetPoint(6, trapmaxtime[index]*sampling, trapmax[index]);
	marker->SetPoint(7, pickoff*sampling, trappick[index]);
	marker->Draw();
	c->Update();
	outfile->cd(("ch"+to_string(iwf)).c_str());
	c->Write();
	tdir->cd();
      }
      if(pz)  delete pz;
      if(pzd) delete pzd;
      if(os)  delete os;
      if(ft)  delete ft;
      if(st)  delete st;
      if(at)  delete at;
    }
    outtree->Fill();
    if(iev < nentries-1) delete wfs;
    if(iev < nentries-1) delete event->GetDigitizerData();
    //event->ClearEventData();
  }
  if(pole_zero) delete pole_zero;
  if(slow_trap) delete slow_trap;
  if(fast_trap) delete fast_trap;
  if(marker)    delete marker;

  // write outputs
  outfile->cd();
  outtree->AutoSave();
  for(int i=0; i<(int)chan_map.size(); i++){
    hdeltat[i]->Write();
    henergy[i]->Write();
    henergyf[i]->Write();
    hbase_energy[i]->Write();
    hbrms_energy[i]->Write();
    if(hdecay_energy[i]) hdecay_energy[i]->Write();
    if(hdecay_deltat[i]) hdecay_deltat[i]->Write();
    hamp_energy[i]->Write();
    hdcr_energy[i]->Write();
    hrise_energy[i]->Write();
  }
  outfile->Close();
    
  return 0;
}
