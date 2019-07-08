#include "utils.hh"
#include <MGTRun.hh>
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
  vector<string> serial_numbers;
  string base_dir = "";
  int max_wf = 0;
  string outfname = "";
  string fnamebase = "";
  string jsonfname = "";
  int run_start = -1;
  int run_end = -1;
  string infname = "";
  bool fit_tail = true;
  bool float_pz = false;
  int nbits = 16;
  int write_wf = 0;
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
			"hd:c:n:o:f:r:R:i:Fb:w:Dj:u:", opts, NULL);
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
      cout << "  -D fit decay constants"             << endl;
      cout << "  -j name of json configuration file" << endl;
      cout << "  -u update interval (percentage)"    << endl;
      return 0;
    case 'd': base_dir  = string(optarg); break;
    case 'c':{
      unsigned i = chan_map.size();
      chan_map[atoi(optarg)] = (int) i;     break;
    }
    case 'n': max_wf      = atoi(optarg);   break;
    case 'o': outfname    = string(optarg); break;
    case 'f': fnamebase   = string(optarg); break;
    case 'r': run_start   = atoi(optarg);   break;
    case 'R': run_end     = atoi(optarg);   break;
    case 'i': infname     = string(optarg); break;
    case 'F': fit_tail    = false;          break;
    case 'b': nbits       = atoi(optarg);   break;
    case 'w': write_wf    = atoi(optarg);   break;
    case 'D': float_pz    = true;           break;
    case 'j': jsonfname   = string(optarg); break;
    case 'u': update_percentage = atoi(optarg); break;
    default: return 1;
    }
    opt = getopt_long(argc, argv,
		      "hd:c:n:o:f:r:R:i:Fb:w:Dj:u:", opts, NULL);
  }
  assert(infname != "" || (run_start > 0 && run_end > 0));

  // default parameters
  vector<double> slow_ramp(chan_map.size(), 2500.0);
  vector<double> slow_flat(chan_map.size(), 1500.0);
  vector<double> fast_ramp(chan_map.size(), 40.0);
  vector<double> fast_flat(chan_map.size(), 100.0);
  vector<double> fast_fall(chan_map.size(), 2000.0);
  vector<double> avse_ramp(chan_map.size(), 100.0);
  vector<double> avse_flat(chan_map.size(), 300.0);
  vector<double> t0_thresh(chan_map.size(), 2.0);
  vector<double> pz_decay(chan_map.size(), 49000.0);
  vector<double> pzd_decay(chan_map.size(), 0.0);
  vector<double> os_amplitude(chan_map.size(), 0.0);
  vector<double> os_decay(chan_map.size(), 0.0);
  vector<int>    nbase_samples(chan_map.size(), 200);
  vector<int>    nefit_samples(chan_map.size(), 300);
  vector<int>    ndcr_samples(chan_map.size(), 100);

  // read the json configuration file if provided
  Json::Value jvalue;
  if(jsonfname != ""){
    ifstream jfile(jsonfname);
    Json::CharReaderBuilder jreader;
    string jerrors;
    if(!Json::parseFromStream(jreader, jfile, &jvalue, &jerrors)){
      cout << jerrors << endl;
      return 2;
    }
    jfile.close();
    vector<int> channels;
    vector<string> detectors;
    SetJson(jvalue, "channel_id", channels);
    SetJson(jvalue, "det_serial", detectors);
    assert(channels.size() == detectors.size());
    if(chan_map.size() == 0){
      cout << "no user specified channel(s), attempting to use channel_id "
	   << "list from configuration file" << endl;
      for(auto const& i : channels){
	unsigned s = chan_map.size();
	chan_map[i] = s;
      }
    }
    for(auto const& pr : chan_map){
      vector<int>::iterator i = find(channels.begin(),
				     channels.end(), pr.first);
      if(i == channels.end()){
	cout <<"channel "<<pr.first<<" not found in "<<jsonfname<<endl;
	return 3;
      }
      int index = distance(channels.begin(), i);
      if(!jvalue.isMember(detectors[index])){
	cout << "configuration not found for " << detectors[index] << " in "
	     << jsonfname << endl;
	return 3;
      }
      serial_numbers.push_back(detectors[index]);
      Json::Value value = jvalue[detectors[index]];
      cout << "reading parameters for channel " << pr.first << " detector "
	   << detectors[index] << endl;
      SetJson(value, "slow_ramp",     slow_ramp[pr.second]);
      SetJson(value, "slow_flat",     slow_flat[pr.second]);
      SetJson(value, "fast_ramp",     fast_ramp[pr.second]);
      SetJson(value, "fast_flat",     fast_flat[pr.second]);
      SetJson(value, "fast_fall",     fast_fall[pr.second]);
      SetJson(value, "avse_ramp",     avse_ramp[pr.second]);
      SetJson(value, "t0_thresh",     t0_thresh[pr.second]);
      SetJson(value, "pzd_decay",     pzd_decay[pr.second]);
      SetJson(value, "os_amplitude",  os_amplitude[pr.second]);
      SetJson(value, "os_decay",      os_decay[pr.second]);
      SetJson(value, "nbase_samples", nbase_samples[pr.second]);
      SetJson(value, "nefit_samples", nefit_samples[pr.second]);
      SetJson(value, "ndcr_samples",  ndcr_samples[pr.second]);
      SetJson(value, "pz_decay",      pz_decay[pr.second]);
    }
  }
  if(chan_map.size() == 0){
    cout << "no user-specified channels or configuration channel list" << endl;
    return 4;
  }
  
  // output tree
  int run, eventnum;
  vector<string> detserial;
  vector<int>    channel, maxtime, mintime, trapmaxtime;
  vector<double> baseline, baserms, maxval, minval, trappick;
  vector<double> t0, t1, t10, t50, t90, t99, imax, dcrslope;
  vector<double> trapmax, time, deltat, stime;
  vector<vector<double> > exp_param;
  TFile* outfile = new TFile(outfname.c_str(), "recreate");
  if(write_wf)
    for(auto const& p : chan_map)
      outfile->mkdir(("ch"+to_string(p.first)).c_str());
  tdir->cd();
  TTree* outtree = new TTree("tree", "tree");
  outtree->Branch("run", &run);
  outtree->Branch("eventnum", &eventnum);
  outtree->Branch("detserial", &detserial);
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
  outtree->Branch("sampling", &stime);
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
  vector<TH2D*> hamp_energyf(chan_map.size(), NULL);
  vector<TH2D*> hdcr_energyf(chan_map.size(), NULL);
  vector<TH2D*> hrise_energyf(chan_map.size(), NULL);
  int adcbins = (1<<nbits) / 16;
  for(auto const& pr : chan_map){
    int ch = pr.first;
    int i = pr.second;
    hdeltat[i] = new TH1D(("hdeltat_"+to_string(ch)).c_str(), "",
			  1e5, 0.0, 10000.0);
    hdeltat[i]->SetXTitle("Delta t (#mu s)");
    hdeltat[i]->SetYTitle("Entries");
    henergy[i] = new TH1D(("henergy_"+to_string(ch)).c_str(), "",
			  1<<nbits, 0.0, 1<<nbits);
    henergy[i]->SetXTitle("Trap Maximum (ADC)");
    henergy[i]->SetYTitle("Entries");
    henergyf[i] = new TH1D(("henergyf_"+to_string(ch)).c_str(), "",
			   1<<nbits, 0.0, 1<<nbits);
    henergyf[i]->SetXTitle("Fixed Time Pickoff (ADC)");
    henergyf[i]->SetYTitle("Entries");
    hbase_energy[i] = new TH2D(("hbase_energy_"+to_string(ch)).c_str(), "",
			       adcbins, 0.0, 10000.0, 400, -20.0, 20.0);
    hbase_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hbase_energy[i]->SetYTitle("Baseline (ADC)");
    hbrms_energy[i] = new TH2D(("hbrms_energy_"+to_string(ch)).c_str(), "",
			       adcbins, 0.0, 10000.0, 400, 0.0, 40.0);
    hbrms_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hbrms_energy[i]->SetYTitle("Baseline RMS (ADC)");
    hdecay_energy[i] = new TH2D(("hdecay_energy_"+to_string(ch)).c_str(),
				"", adcbins, 0.0, 10000.0,
				400, 0.0, 100.0);
    hdecay_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hdecay_energy[i]->SetYTitle("Decay Constant (#mu s)");
    hdecay_deltat[i] = new TH2D(("hdecay_deltat_"+to_string(ch)).c_str(),
				"", 1000, 0.0, 1000.0, 400, 0.0, 100.0);
    hdecay_deltat[i]->SetYTitle("Decay Constant (#mu s)");
    hdecay_deltat[i]->SetXTitle("#Delta t (#mu s)");
    hamp_energy[i] = new TH2D(("hamp_energy_"+to_string(ch)).c_str(), "",
			      adcbins, 0.0, 10000.0, 500, 0.0, 5000.0);
    hamp_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hamp_energy[i]->SetYTitle("Maximum Current");
    hdcr_energy[i] = new TH2D(("hdcr_energy_"+to_string(ch)).c_str(), "",
			      adcbins, 0.0, 10000.0,
			      400, -100.0, 100.0);
    hdcr_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hdcr_energy[i]->SetYTitle("Raw DCR Slope");
    hrise_energy[i] = new TH2D(("hrise_energy_"+to_string(ch)).c_str(), "",
			       2*adcbins, 0.0, 10000.0,
			       1000, 0.0, 5000.0);
    hrise_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hrise_energy[i]->SetYTitle("Rise Time (ns)");
    hamp_energyf[i] = (TH2D*) hamp_energy[i]->Clone(("hamp_energyf_" +
						     to_string(ch)).c_str());
    hdcr_energyf[i] = (TH2D*) hdcr_energy[i]->Clone(("hdcr_energyf_" +
						     to_string(ch)).c_str());
    hrise_energyf[i] =(TH2D*) hrise_energy[i]->Clone(("hrise_energyf_" +
						      to_string(ch)).c_str());
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
  TTreeReaderValue<MGTRun> mgtrun(reader, "run");
  TTreeReaderValue<MGTEvent> event(reader, "event");

  MGWFPoleZeroCorrection* pole_zero=new MGWFPoleZeroCorrection();
  MGWFTrapezoidalFilter* slow_trap=new MGWFTrapezoidalFilter(slow_ramp[0],
							     slow_flat[0]);
  MGWFAsymTrapezoidalFilter* fast_trap = new MGWFAsymTrapezoidalFilter();
  fast_trap->SetRampTime(fast_ramp[0]);
  fast_trap->SetFlatTime(fast_flat[0]);
  fast_trap->SetFallTime(fast_fall[0]);
  MGWFTrapezoidalFilter* avse_trap=new MGWFTrapezoidalFilter(avse_ramp[0],
							     avse_flat[0]);

  int lrun = -1;
  int iev = -1;
  int cpercent = -1;
  int lpercent = -1;
  vector<int> skip_chan;
  while(reader.Next()){
    iev ++;
    cpercent = (int) (100.0*iev/nentries);
    if(cpercent % update_percentage == 0 && cpercent != lpercent)
      cout << cpercent << "%" << endl;
    lpercent = cpercent;
    if(iev >= nentries) break;
    if(iev % 100000 == 0 && iev > 0) outtree->AutoSave();
    run = mgtrun->GetRunNumber();
    // event->GetEventNumber() is 0 for the basic event builder
    if(run != lrun){
      lrun = run;
      eventnum = -1;
    }
    eventnum ++;
    TClonesArray* wfs = event->GetWaveforms();
    int nwf = (int) wfs->GetEntries();
    // clear previous values
    channel.assign(nwf, 0);
    detserial.assign(nwf, "");
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
    stime.assign(nwf, 0.0);
    imax.assign(nwf, 0.0);
    dcrslope.assign(nwf, 0.0);
    trapmax.assign(nwf, 0.0);
    time.assign(nwf, 0.0);
    deltat.assign(nwf, 0.0);
    trappick.assign(nwf, 0.0);
    exp_param.assign(nwf, vector<double>(4, 0.0));
    // start analyzing waveforms for channels selected by user
    for(int iwf=0; iwf<(int)wfs->GetEntriesFast(); iwf++){
      MGTWaveform* wf = (MGTWaveform*) wfs->At(iwf);
      channel[iwf] = (int) wf->GetID();
      if(chan_map.find(channel[iwf]) == chan_map.end()){
	if(find(skip_chan.begin(),
		skip_chan.end(), channel[iwf]) == skip_chan.end()){
	  skip_chan.push_back(channel[iwf]);
	  continue;
	}
      }
      int index = chan_map[channel[iwf]];
      detserial[iwf] = serial_numbers[index];
      wf_count[index] ++;
      bool write = false;
      if(write_wf>0) if(wf_count[index] % write_wf == 0) write = true;
      TH1D* hwf_orig = NULL;
      if(write) hwf_orig = wf->GimmeUniqueHist();
      time[iwf] = (event->GetTime() + wf->GetTOffset())/1e9;
      deltat[iwf] = time[iwf] - tlast[index];
      hdeltat[index]->Fill(deltat[iwf] * 1e6);
      double sampling = wf->GetSamplingPeriod();
      stime[iwf] = sampling;
      vector<double> vwf = wf->GetVectorData();
      assert((int)vwf.size()>nbase_samples[index] &&
	     (int)vwf.size()>=2*nefit_samples[index]);
      baseline[iwf] = accumulate(vwf.begin(),
				 vwf.begin()+nbase_samples[index], 0);
      baseline[iwf] /= nbase_samples[index];
      double base_orig = baseline[iwf];
      // exponential baseline correction from last event, fit exponential tail
      if(fit_tail){
	for_each(xfit.begin(), xfit.end(),
		 [&](double& s){s-=deltat[iwf]*1e9;});
	for(int i=0; i<nefit_samples[index]; i++){
	  xfit.push_back(i*sampling);
	  yfit.push_back(vwf[i]);
	}
	TF1* fn = new TF1("ftmp", "[0]+[1]*exp(-x/[2])",
			  xfit[0], xfit[xfit.size()-1]);
	fn->SetParameter(0, baseline[iwf]);
	fn->SetParameter(1, vwf[0]-baseline[iwf]);
	if(float_pz) fn->SetParameter(2, pz_decay[index]);
	else fn->FixParameter(2, pz_decay[index]);
	TGraph* gr = new TGraph(xfit.size(), &xfit.front(), &yfit.front());
	gr->Fit(fn, "QRNW");
	for(int i=0; i<3; i++) exp_param[iwf][i] = fn->GetParameter(i);
	exp_param[iwf][3] = fn->GetChisquare() / (xfit.size()-3);
	xfit.resize(nefit_samples[index], 0.0);
	yfit.resize(nefit_samples[index], 0.0);
	for(int i=(int)vwf.size()-nefit_samples[index]; i<(int)vwf.size();i++){
	  int j = i - (int)vwf.size() + nefit_samples[index];
	  xfit[j] = i*sampling;
	  yfit[j] = vwf[i];
	}
	for(int i=0; i<(int)vwf.size(); i++)
	  vwf[i] = vwf[i] - fn->Eval(i*sampling);
	delete gr;
	delete fn;
	baseline[iwf] = accumulate(vwf.begin(),
				   vwf.begin()+nbase_samples[index],0);
	baseline[iwf] /= nbase_samples[index];
      }
      // offset baseline subtraction, baseline rms, extrema, and decay constant
      tlast[index] = time[iwf];
      for_each(vwf.begin(), vwf.end(), [&](double& s){s-=baseline[iwf];});
      baserms[iwf]=sqrt(inner_product(vwf.begin(),
				      vwf.begin()+nbase_samples[index],
				      vwf.begin(), 0.0)/nbase_samples[index]);
      wf->SetData(vwf);
      TH1D* hwf_cor = NULL;
      if(write) hwf_cor = wf->GimmeUniqueHist();
      vector<double>::iterator itmax = max_element(vwf.begin(), vwf.end());
      vector<double>::iterator itmin = min_element(vwf.begin(), vwf.end());
      maxtime[iwf] = distance(vwf.begin(), itmax);
      mintime[iwf] = distance(vwf.begin(), itmin);
      maxval[iwf]  = *itmax;
      minval[iwf]  = *itmin;
      if(itmax == vwf.begin() || itmax == vwf.end()) continue;
      // second pass pz correction and adc overshoot correction
      MGTWaveform* pz = new MGTWaveform();
      if(fit_tail) pole_zero->SetDecayConstant(exp_param[iwf][2]);
      else pole_zero->SetDecayConstant(pz_decay[index]);
      pole_zero->TransformOutOfPlace(*wf, *pz);
      vector<double> vwff = pz->GetVectorData();
      MGTWaveform* wff = pz;
      MGTWaveform* pzd = NULL;
      if(pzd_decay[index] > 0.0){
	pzd = new MGTWaveform();
	PZDiff(*pz, *pzd, pzd_decay[index]);
	vwff = pzd->GetVectorData();
	wff = pzd;
      }
      MGTWaveform* os = NULL;
      if(os_decay[index] > 0.0 && os_amplitude[index] > 0.0){
	os = new MGTWaveform();
	OvershootCorrection(*pzd, *os, os_amplitude[index], os_decay[index]);
	vwff = os->GetVectorData();
	wff = os;
      }
      // get time points
      MGTWaveform* ft = new MGTWaveform();
      ft->SetTOffset(fast_fall[index]+fast_flat[index]);
      MGTWaveform* st = new MGTWaveform();
      MGTWaveform* at = new MGTWaveform();
      fast_trap->SetFlatTime(fast_flat[index]);
      fast_trap->SetRampTime(fast_ramp[index]);
      fast_trap->SetFallTime(fast_fall[index]);
      fast_trap->TransformOutOfPlace(*wff, *ft);
      slow_trap->SetFlatTime(slow_flat[index]);
      slow_trap->SetRampTime(slow_ramp[index]);
      slow_trap->TransformOutOfPlace(*wff, *st);
      avse_trap->SetFlatTime(avse_flat[index]);
      avse_trap->SetRampTime(avse_ramp[index]);
      avse_trap->TransformOutOfPlace(*wff, *at);
      vector<double> vft = ft->GetVectorData();
      vector<double> vst = st->GetVectorData();
      vector<double> vat = at->GetVectorData();
      vector<double>::iterator itstmax = max_element(vst.begin(), vst.end());
      vector<double>::iterator itftmax = max_element(vft.begin(), vft.end()-1);
      imax[iwf] = *max_element(vat.begin(), vat.end());
      trapmaxtime[iwf] = distance(vst.begin(), itstmax);
      trapmax[iwf] = *itstmax;
      hbase_energy[index]->Fill(trapmax[iwf], baseline[iwf]);
      hbrms_energy[index]->Fill(trapmax[iwf], baserms[iwf]);
      if(fit_tail){
	if(deltat[iwf]*1e6 < 250)
	  hdecay_energy[index]->Fill(trapmax[iwf], exp_param[iwf][2]/1000);
	hdecay_deltat[index]->Fill(deltat[iwf]*1e6,exp_param[iwf][2]/1000);
      }
      double ftbase = accumulate(vft.begin(),
				 vft.begin()+nbase_samples[index], 0.0);
      ftbase /= nbase_samples[index];
      for_each(vft.begin(), vft.end(), [&](double& s){s-=ftbase;});
      double ftrms = sqrt(inner_product(vft.begin(),
					vft.begin()+nbase_samples[index],
					vft.begin(),0.0)/nbase_samples[index]);
      double ftoffset = ft->GetTOffset()/sampling;
      for(int i=distance(vft.begin(), itftmax); i>=0; i--){
	double val = vft[i] / (*itftmax);
	double diff = vft[i+1] - vft[i];
	if(diff == 0.0) diff = 1.0e12;
	if(val < 0.01 && t1[iwf] == 0.0)
	  t1[iwf] = i + (vft[i+1]-0.01*(*itftmax))/diff + ftoffset;
	if(val < 0.1 && t10[iwf] == 0.0)
	  t10[iwf] = i + (vft[i+1]-0.1*(*itftmax))/diff + ftoffset;
	if(val < 0.5 && t50[iwf] == 0.0)
	  t50[iwf] = i + (vft[i+1]-0.5*(*itftmax))/diff + ftoffset;
	if(val < 0.9 && t90[iwf] == 0.0)
	  t90[iwf] = i + (vft[i+1]-0.9*(*itftmax))/diff + ftoffset;
	if(val < 0.99 && t99[iwf] == 0.0)
	  t99[iwf] = i + (vft[i+1]-0.99*(*itftmax))/diff + ftoffset;
	if(vft[i] < t0_thresh[index]*ftrms){
	  t0[iwf] = i + ftoffset;
	  break;
	}
      }
      t0[iwf] = min(max(0., t0[iwf]), (double)vft.size());
      t1[iwf] = min(max(t0[iwf],   t1[iwf]), (double)vft.size());
      t10[iwf]= min(max(t1[iwf],  t10[iwf]), (double)vft.size());
      t50[iwf]= min(max(t10[iwf], t50[iwf]), (double)vft.size());
      t90[iwf]= min(max(t50[iwf], t90[iwf]), (double)vft.size());
      t99[iwf]= min(max(t90[iwf], t99[iwf]), (double)vft.size());
      // fixed time energy pickoff
      int pickoff = t0[iwf] +
	(int) (slow_ramp[index] + 0.75*slow_flat[index])/sampling;
      if(pickoff < 0 || pickoff >= (int) vst.size()) trappick[iwf] = 0;
      else trappick[iwf] = vst[pickoff];
      henergy[index]->Fill(trapmax[iwf]);
      henergyf[index]->Fill(trappick[iwf]);
      hamp_energy[index]->Fill(trapmax[iwf], imax[iwf]);
      hamp_energyf[index]->Fill(trappick[iwf], imax[iwf]);
      double tr = (t99[iwf] - t1[iwf])*sampling - fast_ramp[index];
      hrise_energy[index]->Fill(trapmax[iwf], tr);
      hrise_energyf[index]->Fill(trappick[iwf], tr);
      // compute dcr
      int dcr_start = min((int)t99[iwf],
			  (int)vwff.size()-2*ndcr_samples[index]-1);
      dcrslope[iwf] = accumulate(vwff.end()-ndcr_samples[index], vwff.end(),0);
      dcrslope[iwf]-= accumulate(vwff.begin()+dcr_start,vwff.begin()+dcr_start+
				 ndcr_samples[index], 0.0);
      dcrslope[iwf] /= vwff.size() - dcr_start - ndcr_samples[index];
      hdcr_energy[index]->Fill(trapmax[iwf], dcrslope[iwf]);
      hdcr_energyf[index]->Fill(trappick[iwf], dcrslope[iwf]);
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
	string cname = "c_"+to_string(iev)+"_"+to_string(channel[iwf]);
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
	marker->SetPoint(0, t0[iwf]*sampling, ft->At((int)(t0[iwf]-fto)));
	marker->SetPoint(1, t1[iwf]*sampling, ft->At((int)(t1[iwf]-fto)));
	marker->SetPoint(2,t10[iwf]*sampling, ft->At((int)(t10[iwf]-fto)));
	marker->SetPoint(3,t50[iwf]*sampling, ft->At((int)(t50[iwf]-fto)));
	marker->SetPoint(4,t90[iwf]*sampling, ft->At((int)(t90[iwf]-fto)));
	marker->SetPoint(5,t99[iwf]*sampling, ft->At((int)(t99[iwf]-fto)));
	marker->SetPoint(6, trapmaxtime[iwf]*sampling, trapmax[iwf]);
	marker->SetPoint(7, pickoff*sampling, trappick[iwf]);
	marker->Draw();
	c->Update();
	outfile->cd(("ch"+to_string(channel[iwf])).c_str());
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
    //if(iev < nentries-1) delete wfs;
    //if(iev < nentries-1) delete event->GetDigitizerData();
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
    hamp_energyf[i]->Write();
    hdcr_energyf[i]->Write();
    hrise_energyf[i]->Write();
  }
  int ch_count = 0;
  for(auto const& p : chan_map){
    jvalue["channel_map"][ch_count] = p.first;
    ch_count ++;
  }
  Json::StreamWriterBuilder builder;
  builder["indentation"] = "";
  TObjString* jout=new TObjString(Json::writeString(builder, jvalue).c_str());
  jout->Write("process_config");
  outfile->Close();

  if(skip_chan.size() != 0){
    cout << "did not process channel(s) ";
    for(auto const& ch : skip_chan)
      cout << ch << " ";
    cout << endl;
  }
  
  return 0;
}
