#include "MultiWaveform.hh"
#include "utils.hh"
#ifdef __CUDA
#include "gutils.hh"
#endif
#include <MGTRun.hh>
#include <MGTEvent.hh>
#include <MGTWaveform.hh>
#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLine.h>
#include <TPolyMarker.h>
#include <TStopwatch.h>
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
#include <future>
#include <thread>
#include <chrono>

using namespace std;
using namespace std::chrono;

TDirectory* tdir;

enum MultiWaveformState{
  EMPTY,
  POPULATING,
  WAITING,
  CPU_PROCESSING,
  GPU_PROCESSING,
  COMPLETE,
  WRITTEN
};

int main(int argc, char* argv[]){

  vector<double> times(7, 0.0);
  TStopwatch twatch, pwatch;
  twatch.Start();
  pwatch.Start();
  
  // root setup
  gStyle->SetOptStat(0);
  gStyle->SetLabelFont(132, "XYZ");
  gStyle->SetTitleFont(132, "XYZ");
  tdir = gROOT->CurrentDirectory();
  
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
  bool float_ct = false;
  const int nct_steps = 20;
  int nbits = 16;
  int write_wf = 0;
  int update_percentage = 5;
  int process_block = 100000;
  int max_threads = 4;
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
    {"bits",       required_argument, NULL, 'b'},
    {"writewf",    required_argument, NULL, 'w'},
    {"decayconst",       no_argument, NULL, 'D'},
    {"chargetrap",       no_argument, NULL, 'C'},
    {"jsoncofig",  required_argument, NULL, 'j'},
    {"update",     required_argument, NULL, 'u'},
    {"process_block", required_argument, NULL, 'p'},
    {"threads",       required_argument, NULL, 't'}
  };
  int opt = getopt_long(argc, argv,
			"hd:c:n:o:f:r:R:i:b:w:Cj:u:p:t:", opts, NULL);
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
      cout << "  -b number of ADC bits"              << endl;
      cout << "  -w plot every n'th wf"              << endl;
      cout << "  -C vary charge trapping correction" << endl;
      cout << "  -j name of json configuration file" << endl;
      cout << "  -u update interval (percentage)"    << endl;
      cout << "  -p number of wf to process in block"<< endl;
      cout << "  -t max cpu threads to use"          << endl;
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
    case 'b': nbits       = atoi(optarg);   break;
    case 'w': write_wf    = atoi(optarg);   break;
    case 'C': float_ct    = true;           break;
    case 'j': jsonfname   = string(optarg); break;
    case 'u': update_percentage = atoi(optarg); break;
    case 'p': process_block     = atoi(optarg); break;
    case 't': max_threads       = atoi(optarg); break;
    default: return 1;
    }
    opt = getopt_long(argc, argv,
		      "hd:c:n:o:f:r:R:i:b:w:Cj:u:p:t:", opts, NULL);
  }
  assert(infname != "" || (run_start > 0 && run_end > 0));

  // multi-waveforms for processing and default analysis parameters
  vector<vector<MultiWaveform*> >
    mwf(chan_map.size(), vector<MultiWaveform*>(max_threads, NULL));
  vector<vector<int> > wfcount(max_threads);
  vector<future<vector<map<string, MultiWaveform*> > > > threads(max_threads);
  vector<MultiWaveformState> sthread(max_threads, EMPTY);
  vector<int> bthread(max_threads, -1);
  int cthread = 0;
  int block_count = 0;
  sthread[cthread] = POPULATING;
  bthread[cthread] = block_count;
  for(auto const& pr : chan_map)
    for(auto& wf : mwf[pr.second]){
      wf = new MultiWaveform(pr.first);
      wf->SetID(pr.first);
      wf->SetParam("slow_ramp", 2500.0);
      wf->SetParam("slow_flat", 1500.0);
      wf->SetParam("fast_ramp", 40.0);
      wf->SetParam("fast_flat", 100.0);
      wf->SetParam("fast_fall", 2000.0);
      wf->SetParam("avse_ramp", 100.0);
      wf->SetParam("avse_flat", 300.0);
      wf->SetParam("t0_thresh", 2.0);
      wf->SetParam("pz_decay", 58000.0);
      wf->SetParam("ct_offset", 2000.0);
      wf->SetParam("ct_time", 100000.0);
      wf->SetParam("ct_frac", 1.0e-8);
      wf->SetParam("resting_base", 0.0);
      wf->SetParam("nbase_samples", 200);
      wf->SetParam("nefit_samples", 300);
      wf->SetParam("ndcr_samples", 100);
      wf->SetParam("ct_method", 2);
    }
  
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
      vector<int>::iterator i = std::find(channels.begin(),
				     channels.end(), pr.first);
      if(i == channels.end()){
	cout <<"channel "<<pr.first<<" not found in "<<jsonfname<<endl;
	return 3;
      }
      int index = std::distance(channels.begin(), i);
      if(!jvalue.isMember(detectors[index])){
	cout << "configuration not found for " << detectors[index] << " in "
	     << jsonfname << endl;
	return 3;
      }
      serial_numbers.push_back(detectors[index]);
      Json::Value value = jvalue[detectors[index]];
      cout << "reading parameters for channel " << pr.first << " detector "
	   << detectors[index] << endl;
      bool v = true;
      for(auto& wf : mwf[index]){
	wf->SetParam("slow_rise", GetJsonD(value, "slow_ramp", v));
	wf->SetParam("slow_flat", GetJsonD(value, "slow_flat", v));
	wf->SetParam("fast_rise", GetJsonD(value, "fast_ramp", v));
	wf->SetParam("fast_flat", GetJsonD(value, "fast_flat", v));
	wf->SetParam("fast_fall", GetJsonD(value, "fast_fall", v));
	wf->SetParam("avse_rise", GetJsonD(value, "avse_ramp", v));
	wf->SetParam("avse_flat", GetJsonD(value, "avse_flat", v));
	wf->SetParam("pz_decay",  GetJsonD(value, "pz_decay", v));
	wf->SetParam("t0_thresh", GetJsonD(value, "t0_thresh", v));
	wf->SetParam("ct_offset", GetJsonD(value, "ct_offset", v));
	wf->SetParam("ct_method", GetJsonI(value, "ct_method", v));
	wf->SetParam("ct_decay",  GetJsonD(value, "ct_decay", v));
	wf->SetParam("ct_frac",   GetJsonD(value, "ct_frac", v));
	wf->SetParam("resting_base",  GetJsonD(value, "resting_base", v));
	wf->SetParam("nbase_samples", GetJsonI(value, "nbase_samples", v));
	wf->SetParam("nefit_samples", GetJsonI(value, "nefit_samples", v));
	wf->SetParam("ndcr_samples",  GetJsonI(value, "ndcr_samples", v));
	v = false;
      }
    }
  }
  if(chan_map.size() == 0){
    cout << "no user-specified channels or configuration channel list" << endl;
    return 4;
  }
  for(auto const& vwf : mwf)
    for(auto const& w : vwf)
      if(w->GetParam("ct_method")<1 || w->GetParam("ct_method")>2){
	cout << "CT method must be 1 or 2" << endl;
	return 5;
      }
  for(auto& vwf : mwf)
    for(auto& wf : vwf){
      if(float_ct){
	wf->SetParam("nct_steps", nct_steps);
	for(int ict=0; ict<nct_steps; ict++)
	  wf->SetParam("ct_decay_"+to_string(ict), wf->GetParam("pz_decay") *
		      (0.5+1.5*((float)ict)/nct_steps));
      }
      else wf->SetParam("nct_steps", 0.0);
    }
    
  // output tree
  int run, eventnum;
  vector<string> detserial;
  vector<int>    channel;
  vector<double> baseline, baserms, trappick, ct1_trappick, ct2_trappick;
  vector<double> t0, t1, t10, t50, t90, t99, imax, dcrslope, pickoff;
  vector<double> trapmax, time, deltat, stime, ct_integral, decay_const;
  vector<vector<double> > ct_decay, ct_value;
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
  outtree->Branch("baseline", &baseline);
  outtree->Branch("baserms", &baserms);
  outtree->Branch("t0", &t0);
  outtree->Branch("t1", &t1);
  outtree->Branch("t10", &t10);
  outtree->Branch("t50", &t50);
  outtree->Branch("t90", &t90);
  outtree->Branch("t99", &t99);
  outtree->Branch("pickoff", &pickoff);
  outtree->Branch("imax", &imax);
  outtree->Branch("dcrslope", &dcrslope);
  outtree->Branch("trapmax", &trapmax);
  outtree->Branch("trappick", &trappick);
  outtree->Branch("ct1_trappick", &ct1_trappick);
  outtree->Branch("ct2_trappick", &ct2_trappick);
  outtree->Branch("time", &time);
  outtree->Branch("deltat", &deltat);
  outtree->Branch("sampling", &stime);
  outtree->Branch("decay_const", &decay_const);
  outtree->Branch("ct_integral", &ct_integral);
  outtree->Branch("ct_decay", &ct_decay);
  outtree->Branch("ct_value", &ct_value);
  outtree->SetDirectory(outfile);

  // output histograms
  vector<TH1D*> hdeltat(chan_map.size(), NULL);
  vector<TH1D*> henergy(chan_map.size(), NULL);
  vector<TH1D*> henergyf(chan_map.size(), NULL);
  vector<TH1D*> henergyc1(chan_map.size(), NULL);
  vector<TH1D*> henergyc2(chan_map.size(), NULL);
  vector<TH2D*> hbase_energy(chan_map.size(), NULL);
  vector<TH2D*> hbase_deltat(chan_map.size(), NULL);
  vector<TH2D*> hbrms_deltat(chan_map.size(), NULL);
  vector<TH2D*> hbrms_energy(chan_map.size(), NULL);
  vector<TH2D*> hdecay_energy(chan_map.size(), NULL);
  vector<TH2D*> hdecay_deltat(chan_map.size(), NULL);
  vector<TH2D*> hamp_energy(chan_map.size(), NULL);
  vector<TH2D*> haoe_energy(chan_map.size(), NULL);
  vector<TH2D*> hdcr_energy(chan_map.size(), NULL);
  vector<TH2D*> hrise_energy(chan_map.size(), NULL);
  vector<TH2D*> hamp_energyf(chan_map.size(), NULL);
  vector<TH2D*> haoe_energyf(chan_map.size(), NULL);
  vector<TH2D*> hdcr_energyf(chan_map.size(), NULL);
  vector<TH2D*> hrise_energyf(chan_map.size(), NULL);
  vector<TH2D*> hamp_energyc(chan_map.size(), NULL);
  vector<TH2D*> haoe_energyc(chan_map.size(), NULL);
  vector<TH2D*> hdcr_energyc(chan_map.size(), NULL);
  vector<TH2D*> hrise_energyc(chan_map.size(), NULL);
  int adcbins = (1<<nbits) / 16;
  for(auto const& pr : chan_map){
    int ch = pr.first;
    int i = pr.second;
    string s = to_string(ch);
    hdeltat[i] = new TH1D(("hdeltat_"+s).c_str(), "", 1e5, 0.0, 10000.0);
    hdeltat[i]->SetXTitle("Delta t (#mu s)");
    hdeltat[i]->SetYTitle("Entries");
    henergy[i] = new TH1D(("henergy_"+s).c_str(), "",
			  2*(1<<nbits), 0.0, 1<<nbits);
    henergy[i]->SetXTitle("Trap Maximum (ADC)");
    henergy[i]->SetYTitle("Entries");
    henergyf[i] = new TH1D(("henergyf_"+s).c_str(), "",
			   2*(1<<nbits), 0.0, 1<<nbits);
    henergyf[i]->SetXTitle("Fixed Time Pickoff (ADC)");
    henergyf[i]->SetYTitle("Entries");
    henergyc1[i] = new TH1D(("henergyc1_"+s).c_str(), "",
			    2*(1<<nbits), 0.0, 1<<nbits);
    henergyc1[i]->SetXTitle("CT1 Fixed Time Pickoff (ADC)");
    henergyc1[i]->SetYTitle("Entries");
    henergyc2[i] = new TH1D(("henergyc2_"+s).c_str(), "",
			    2*(1<<nbits), 0.0, 1<<nbits);
    henergyc2[i]->SetXTitle("CT2 Fixed Time Pickoff (ADC)");
    henergyc2[i]->SetYTitle("Entries");
    hbase_energy[i] = new TH2D(("hbase_energy_"+s).c_str(), "",
			       adcbins, 0.0, 10000.0, 400, -20.0, 20.0);
    hbase_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hbase_energy[i]->SetYTitle("Baseline (ADC)");
    hbase_deltat[i] = new TH2D(("hbase_deltat_"+s).c_str(), "",
			       1000, 0.0, 1000.0, 1000, -20.0, 20.0);
    hbase_deltat[i]->SetXTitle("#Deltat t (#mu s)");
    hbase_deltat[i]->SetYTitle("Baseline (ADC)");
    hbrms_deltat[i] = new TH2D(("hbrms_deltat_"+s).c_str(), "",
				1000, 0.0, 1000.0, 100, 0.0, 40.0);
    hbrms_deltat[i]->SetXTitle("#Deltat t (#mu s)");
    hbrms_deltat[i]->SetYTitle("Baseline RMS (ADC)");
    hbrms_energy[i] = new TH2D(("hbrms_energy_"+s).c_str(), "",
			       adcbins, 0.0, 10000.0, 400, 0.0, 40.0);
    hbrms_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hbrms_energy[i]->SetYTitle("Baseline RMS (ADC)");
    hdecay_energy[i] = new TH2D(("hdecay_energy_"+s).c_str(), "",
				adcbins, 0.0, 1<<nbits, 400, 0.0, 100.0);
    hdecay_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hdecay_energy[i]->SetYTitle("Decay Constant (#mu s)");
    hdecay_deltat[i] = new TH2D(("hdecay_deltat_"+s).c_str(), "",
				1000, 0.0, 1000.0, 400, 0.0, 100.0);
    hdecay_deltat[i]->SetYTitle("Decay Constant (#mu s)");
    hdecay_deltat[i]->SetXTitle("#Delta t (#mu s)");
    hamp_energy[i] = new TH2D(("hamp_energy_"+s).c_str(), "",
			      adcbins, 0.0, 10000.0, 500, 0.0, 5000.0);
    hamp_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hamp_energy[i]->SetYTitle("Maximum Current");
    haoe_energy[i] = new TH2D(("haoe_energy_"+s).c_str(), "",
			      adcbins, 0.0, 10000.0, 500, 0.0, 2.0);
    haoe_energy[i]->SetXTitle("Trap Maximum (ADC)");
    haoe_energy[i]->SetYTitle("Maximum Current / Trap Maximum");
    hdcr_energy[i] = new TH2D(("hdcr_energy_"+s).c_str(), "",
			      adcbins, 0.0, 10000.0, 400, -100.0, 100.0);
    hdcr_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hdcr_energy[i]->SetYTitle("Raw DCR Slope");
    hrise_energy[i] = new TH2D(("hrise_energy_"+s).c_str(), "",
			       2*adcbins, 0.0, 10000.0, 1000, 0.0, 5000.0);
    hrise_energy[i]->SetXTitle("Trap Maximum (ADC)");
    hrise_energy[i]->SetYTitle("Rise Time (ns)");
    s = "_energyf_" + to_string(ch);
    hamp_energyf[i] = (TH2D*) hamp_energy[i]->Clone(("hamp"+s).c_str());
    haoe_energyf[i] = (TH2D*) haoe_energy[i]->Clone(("haoe"+s).c_str());
    hdcr_energyf[i] = (TH2D*) hdcr_energy[i]->Clone(("hdcr"+s).c_str());
    hrise_energyf[i] =(TH2D*) hrise_energy[i]->Clone(("hrise"+s).c_str());
    s = "_energyc_" + to_string(ch);
    hamp_energyc[i] = (TH2D*) hamp_energy[i]->Clone(("hamp"+s).c_str());
    haoe_energyc[i] = (TH2D*) haoe_energy[i]->Clone(("haoe"+s).c_str());
    hdcr_energyc[i] = (TH2D*) hdcr_energy[i]->Clone(("hdcr"+s).c_str());
    hrise_energyc[i] =(TH2D*) hrise_energy[i]->Clone(("hrise"+s).c_str());
  }

  // input reader
  TChain tree("MGTree", "MGTree");                                  
  string fname = base_dir + "/" + fnamebase;                
  if(infname != "") tree.Add(infname.c_str());
  else
    for(int irun=run_start; irun<=run_end; irun++)
      tree.Add((fname + to_string(irun) + ".root").c_str());
  int nentries = (int) tree.GetEntries();
  if(max_wf > 0) nentries = std::min(nentries, max_wf);
  TTreeReader reader(&tree);
  TTreeReaderValue<MGTRun> mjtrun(reader, "run");
  TTreeReaderValue<MGTEvent> event(reader, "event");

  // start reading/processing events
  times[0] += pwatch.RealTime();
  pwatch.Start();
  int lrun = -1;
  int iev = -1;
  int cpercent = -1;
  int lpercent = -1;
  int lprocess = 0;
  int event_count = -1;
  vector<int> wf_count(chan_map.size(), -1);
  vector<int> skip_chan;
  vector<vector<double> > baselines(chan_map.size());
  vector<bool> rebin_base(chan_map.size(), false);
  bool read_events = true;
  bool processing_complete = false;
  map<int, int> ev_offset;
  while(!processing_complete){
    if(read_events && iev < nentries-1){
      reader.Next();
      iev ++;
      cpercent = (int) (100.0*iev/nentries);
      if(cpercent % update_percentage == 0 && cpercent != lpercent)
	cout << cpercent << "%" << endl;
      lpercent = cpercent;
      if(iev % 10000 == 0 && iev > 0) outtree->AutoSave();
      run = mjtrun->GetRunNumber();
      if(run != lrun){
	lrun = run;
	event_count = -1;
	ev_offset[run] = iev;
      }
      event_count ++;
      TClonesArray* wfs = event->GetWaveforms();
      wfcount[cthread].push_back((int)wfs->GetEntriesFast());
      // get waveforms for channels selected by user
      for(int iwf=0; iwf<(int)wfs->GetEntriesFast(); iwf++){
	MGTWaveform* wf = (MGTWaveform*) wfs->At(iwf);
	int chan = (int) wf->GetID();
	if(chan_map.find(chan) == chan_map.end()){
	  if(std::find(skip_chan.begin(),
		       skip_chan.end(), chan) == skip_chan.end())
	    skip_chan.push_back(chan);
	  continue;
	}
	int i = chan_map[chan];
	int j = cthread;
	int wfi = mwf[i][j]->GetNWaveforms();
	mwf[i][j]->AddWaveform(wf);
	mwf[i][j]->SetWFParam(wfi,  "run",      (float) run);
	mwf[i][j]->SetWFParam(wfi,  "eventnum", (float) event_count);
	mwf[i][j]->SetWFParam(wfi,  "wfindex",  (float) iwf);
	mwf[i][j]->SetWFDParam(wfi, "time", event->GetTime()+wf->GetTOffset());
	if(wfi > 0)
	  mwf[i][j]->SetWFDParam(wfi, "tlast",
				 mwf[i][j]->GetWFDParam(wfi-1, "time"));
	else{
	  float s = wf->GetSamplingPeriod();
	  mwf[i][j]->SetParam("sampling", s);
	  mwf[i][j]->SetParam("slow_nrise",mwf[i][j]->GetParam("slow_rise")/s);
	  mwf[i][j]->SetParam("slow_nflat",mwf[i][j]->GetParam("slow_flat")/s);
	  mwf[i][j]->SetParam("fast_nrise",mwf[i][j]->GetParam("fast_rise")/s);
	  mwf[i][j]->SetParam("fast_nflat",mwf[i][j]->GetParam("fast_flat")/s);
	  mwf[i][j]->SetParam("fast_nfall",mwf[i][j]->GetParam("fast_fall")/s);
	  mwf[i][j]->SetParam("avse_nrise",mwf[i][j]->GetParam("avse_rise")/s);
	  mwf[i][j]->SetParam("avse_nflat",mwf[i][j]->GetParam("avse_flat")/s);
	  mwf[i][j]->SetWFDParam(wfi, "tlast", 0.0);
	}
      }    
      // check if we should process waveforms, if so, process and fit
      // if the gpu has enough memory, launch there, otherwise on cpu
      bool process = false;
      if(iev >= nentries-1) process = true;
      else{
	for(int i=0; i<(int)mwf.size(); i++)
	  if(mwf[i][cthread]->GetNWaveforms() >= process_block) process = true;
      }
      if(!process) continue;
      times[1] += pwatch.RealTime();
      pwatch.Start();
      if(block_count == 0)
	cout<<"approximately "<<1+nentries/iev<<" wf blocks to process"<<endl;
      #ifdef __CUDA
      size_t max_size = 0;
      for(unsigned i=0; i<chan_map.size(); i++){
	size_t s = 3 * sizeof(float)*mwf[i][cthread]->wfy.size();
	if(write_wf) s *= 2;
	max_size = max(max_size, s);
      }
      if(max_size < GPUMemory()){
	int c = cthread;
	threads[cthread] = async(launch::async, [&]{
	    vector<map<string, MultiWaveform*> > twf(chan_map.size());
	    for(unsigned i=0; i<chan_map.size(); i++)
	    twf[i] = ProcessMultiWaveformGPU(mwf[i][c], write_wf>0);
	    return twf;});
	bthread[cthread] = block_count;
	sthread[cthread] = GPU_PROCESSING;
	cout << "processing block " << block_count
	     << " on GPU with thread " << cthread << endl;
      }
      #endif
      if(sthread[cthread] != GPU_PROCESSING){
	int c = cthread;
	threads[cthread] = async(launch::async, [&]{
	    vector<map<string, MultiWaveform*> > twf(chan_map.size());
	    for(unsigned i=0; i<chan_map.size(); i++)
	      twf[i] = ProcessMultiWaveform(mwf[i][c], write_wf>0, bthread[c]);
	    return twf;});
	bthread[cthread] = block_count;
	sthread[cthread] = CPU_PROCESSING;
	cout << "processing block " << block_count
	     << " on CPU with thread " << cthread << endl;
	this_thread::sleep_for(milliseconds(10));
      }
    }
    // check for finished blocks, wait for threads to finish if last event
    int min_proc = 1e9;
    processing_complete = true;
    for(unsigned i=0; i<sthread.size(); i++){
      if(sthread[i] == CPU_PROCESSING || sthread[i] == GPU_PROCESSING){
	processing_complete = false;
	if(threads[i].wait_for(milliseconds(0)) == future_status::ready)
	  sthread[i] = COMPLETE;
	else{
	  if(iev < nentries-1) min_proc = min(min_proc, bthread[i]);
	  else{
	    cout << "waiting for thread " << i << " to finish..." << endl;
	    TStopwatch wtmp;
	    wtmp.Start();
	    threads[i].wait();
	    sthread[i] = COMPLETE;
	    times[5] -= wtmp.RealTime();
	  }
	}
      }
      else if(sthread[i] == COMPLETE) processing_complete = false;
    }
    times[5] += pwatch.RealTime();
    pwatch.Start();
    // check for completed threads to write
    int mw_index = 0;
    int mw_block = 1e9;
    while(mw_index >= 0){
      mw_index = -1;
      mw_block = 1e9;
      for(unsigned ith=0; ith<threads.size(); ith++)
	if(sthread[ith] == COMPLETE && bthread[ith] < min_proc)
	  if(bthread[ith] < mw_block){
	    mw_block = bthread[ith];
	    mw_index = ith;
	  }
      if(mw_index < 0) continue;
      sthread[mw_index] = WRITTEN;
      cout << "writing events from block " << bthread[mw_index]
	   << ", thread " << mw_index << endl;
      vector<map<string, MultiWaveform*> > twf = threads[mw_index].get();
      // fill histograms, tree, and plot
      vector<int> chindex(chan_map.size(), 0);
      for(int jev=0; jev<(int)wfcount[mw_index].size(); jev++){
	int nwf = wfcount[mw_index][jev];
	if(nwf > 1) cout << nwf << endl;
	// clear previous values
	channel.assign(nwf, 0);
	detserial.assign(nwf, "");
	baseline.assign(nwf, 0.0);
	baserms.assign(nwf, 0.0);
	t0.assign(nwf, 0.0);
	t1.assign(nwf, 0.0);
	t10.assign(nwf, 0.0);
	t50.assign(nwf, 0.0);
	t90.assign(nwf, 0.0);
	t99.assign(nwf, 0.0);
	pickoff.assign(nwf, 0.0);
	stime.assign(nwf, 0.0);
	imax.assign(nwf, 0.0);
	dcrslope.assign(nwf, 0.0);
	trapmax.assign(nwf, 0.0);
	time.assign(nwf, 0.0);
	deltat.assign(nwf, 0.0);
	decay_const.assign(nwf, 0.0);
	trappick.assign(nwf, 0.0);
	ct1_trappick.assign(nwf, 0.0);
	ct2_trappick.assign(nwf, 0.0);
	ct_integral.assign(nwf, 0.0);
	if(float_ct){
	  ct_decay.assign(nwf, vector<double>(nct_steps, 0.0));
	  ct_value.assign(nwf, vector<double>(nct_steps, 0.0));
	}
	else{
	  ct_decay.resize(0);
	  ct_value.resize(0);
	  }
	int tmprun = run;
	for(int ich=0; ich<(int)mwf.size(); ich++){
	  int ith = mw_index;
	  if(jev == 0){
	    int ptype = (int) mwf[ich][ith]->GetParam("proc_type");
	    if(ptype == 0) times[2] += mwf[ich][ith]->GetParam("proc_time");
	    else times[6] += mwf[ich][ith]->GetParam("proc_time");
	  }
	  int ev = (int) mwf[ich][ith]->GetWFParam(chindex[ich], "eventnum");
	  ev += ev_offset[mwf[ich][ith]->GetWFParam(chindex[ich], "run")];
	  while(chindex[ich]<mwf[ich][ith]->GetNWaveforms()){
	    if(ev >= lprocess + jev) break;
	    chindex[ich] ++;
	    ev = (int) mwf[ich][ith]->GetWFParam(chindex[ich], "eventnum");
	    ev += ev_offset[mwf[ich][ith]->GetWFParam(chindex[ich], "run")];
	  }
	  if(chindex[ich] >= mwf[ich][ith]->GetNWaveforms()) continue;
	  if(ev != lprocess+jev) continue;
	  int iwf = (int) mwf[ich][ith]->GetWFParam(chindex[ich], "wfindex");
	  for(auto const& pr : chan_map)
	    if(pr.second == ich){
	      channel[iwf] = pr.first;
	      break;
	    }
	  detserial[iwf] = serial_numbers[ich];
	  wf_count[ich] ++;
	  bool write = false;
	  if(write_wf>0) if(wf_count[ich] % write_wf == 0) write = true;
	  int jwf = chindex[ich];
	  eventnum          = (int) mwf[ich][ith]->GetWFParam(jwf, "eventnum");
	  run               = (int) mwf[ich][ith]->GetWFParam(jwf, "run");
	  channel[iwf]      = mwf[ich][ith]->GetID();
	  stime[iwf]        = mwf[ich][ith]->GetParam("sampling");
	  time[iwf]         = mwf[ich][ith]->GetWFDParam(jwf, "time");
	  deltat[iwf]     = time[iwf]-mwf[ich][ith]->GetWFDParam(jwf, "tlast");
	  baseline[iwf]     = mwf[ich][ith]->GetWFParam(jwf, "base");
	  baserms[iwf]      = mwf[ich][ith]->GetWFParam(jwf, "base_rms");
	  trapmax[iwf]      = mwf[ich][ith]->GetWFParam(jwf, "trap_max");
	  trappick[iwf]     = mwf[ich][ith]->GetWFParam(jwf, "trappick");
	  imax[iwf]         = mwf[ich][ith]->GetWFParam(jwf, "imax");
	  t0[iwf]           = mwf[ich][ith]->GetWFParam(jwf, "t0");
	  t1[iwf]           = mwf[ich][ith]->GetWFParam(jwf, "t1");
	  t10[iwf]          = mwf[ich][ith]->GetWFParam(jwf, "t10");
	  t50[iwf]          = mwf[ich][ith]->GetWFParam(jwf, "t50");
	  t90[iwf]          = mwf[ich][ith]->GetWFParam(jwf, "t90");
	  t99[iwf]          = mwf[ich][ith]->GetWFParam(jwf, "t99");
	  pickoff[iwf]      = mwf[ich][ith]->GetWFParam(jwf, "pickoff");
	  ct_integral[iwf]  = mwf[ich][ith]->GetWFParam(jwf, "ct_integral");
	  ct1_trappick[iwf] = mwf[ich][ith]->GetWFParam(jwf, "ct1_trappick");
	  ct2_trappick[iwf] = mwf[ich][ith]->GetWFParam(jwf, "ct2_trappick");
	  dcrslope[iwf]     = mwf[ich][ith]->GetWFParam(jwf, "dcrslope");
	  decay_const[iwf]  = mwf[ich][ith]->GetWFParam(jwf, "decay_const");
	  if(!rebin_base[ich]){
	    // rebin baseline histograms based on the first 100 events
	    if(baselines[ich].size() < 100)
	      baselines[ich].push_back(baseline[iwf]);
	    else{
	      rebin_base[ich] = true;
	      delete hbase_energy[ich];
	      delete hbase_deltat[ich];
	      vector<double> v = baselines[ich];
	      double mean = accumulate(v.begin(), v.end(), 0.0)/v.size();
	      for_each(v.begin(), v.end(), [&](double&s){s-=mean;});
	      double rms = sqrt(inner_product(v.begin(), v.end(),
					      v.begin(), 0.0)/v.size());
	      string s = "hbase_energy_" + to_string(channel[iwf]);
	      hbase_energy[ich] = new TH2D(s.c_str(), "",
					   adcbins, 0.0, 10000.0,
					   400, mean-10*rms, mean+10*rms);
	      hbase_energy[ich]->SetXTitle("Trap Maximum (ADC)");
	      hbase_energy[ich]->SetYTitle("Baseline (ADC)");
	      s = "hbase_deltat_" + to_string(channel[iwf]);
	      hbase_deltat[ich] = new TH2D(s.c_str(), "",
					   1000, 0.0, 1000.0,
					   1000, mean-10*rms, mean+10*rms);
	      hbase_deltat[ich]->SetXTitle("#Delta t (#mu s)");
	      hbase_deltat[ich]->SetYTitle("Baseline (ADC)");
	      baselines[ich].clear();
	    }
	  }
	  hdeltat[ich]->Fill(deltat[iwf]/1.e3);
	  hbase_energy[ich]->Fill(trapmax[iwf],     baseline[iwf]);
	  hbase_deltat[ich]->Fill(deltat[iwf]/1.e3, baseline[iwf]);
	  hbrms_deltat[ich]->Fill(deltat[iwf]/1.e3, baserms[iwf]);
	  hbrms_energy[ich]->Fill(trapmax[iwf],     baserms[iwf]);
	  hdecay_energy[iwf]->Fill(trapmax[iwf],    decay_const[iwf]/1000);
	  hdecay_deltat[iwf]->Fill(deltat[iwf]/1e3, decay_const[iwf]/1000);
	  double ctE = ct1_trappick[iwf];
	  if((int) mwf[ich][ith]->GetParam("ct_method")==2)
	    ctE = ct2_trappick[iwf];
	  henergy[ich]->Fill(trapmax[iwf]);
	  henergyf[ich]->Fill(trappick[iwf]);
	  henergyc1[ich]->Fill(ct1_trappick[iwf]);
	  henergyc2[ich]->Fill(ct2_trappick[iwf]);
	  hamp_energy[ich]->Fill(trapmax[iwf], imax[iwf]);
	  hamp_energyf[ich]->Fill(trappick[iwf], imax[iwf]);
	  hamp_energyc[ich]->Fill(ctE, imax[iwf]);
	  haoe_energy[ich]->Fill(trapmax[iwf], imax[iwf]/trapmax[iwf]);
	  haoe_energyf[ich]->Fill(trappick[iwf], imax[iwf]/trappick[iwf]);
	  haoe_energyc[ich]->Fill(ctE, imax[iwf]/ctE);
	  double tr = (t99[iwf] - t1[iwf])*stime[iwf];
	  tr -= mwf[ich][ith]->GetParam("fast_ramp");
	  hrise_energy[ich]->Fill(trapmax[iwf], tr);
	  hrise_energyf[ich]->Fill(trappick[iwf], tr);
	  hrise_energyc[ich]->Fill(ctE, tr);
	  hdcr_energy[ich]->Fill(trapmax[iwf], dcrslope[iwf]);
	  hdcr_energyf[ich]->Fill(trappick[iwf], dcrslope[iwf]);
	  hdcr_energyc[ich]->Fill(ctE, dcrslope[iwf]);
	  if(float_ct)
	    for(int i=0; i<=nct_steps; i++){
	      ct_decay[iwf][i] = mwf[ich][ith]->GetParam("ct_decay_" +
							to_string(i));
	      ct_value[iwf][i] =
		mwf[ich][ith]->GetWFParam(jwf, "ct1_trappick_"+to_string(i));
	    }
	  // plot waveform and write to output file
	  if(write){
	    string s = "wf_"+to_string(lprocess+jev)+"_"+to_string(ich)+"_";
	    double sampling = mwf[ich][ith]->GetParam("sampling");
	    TH1D* hwf = GetWFHist(twf[ich]["base_sub"],  jwf, s+"wf");
	    TH1D* hpz = GetWFHist(twf[ich]["pz_cor"],    jwf, s+"pz");
	    TH1D* hft = GetWFHist(twf[ich]["fast_trap"], jwf, s+"ft");
	    TH1D* hst = GetWFHist(twf[ich]["slow_trap"], jwf, s+"st");
	    TH1D* hat = GetWFHist(twf[ich]["avse_trap"], jwf, s+"at");
	    TH1D* hct = GetWFHist(twf[ich]["pz_ct"],     jwf, s+"ct");
	    TPolyMarker* marker = new TPolyMarker(7);
	    marker->SetMarkerSize(0.8);
	    marker->SetMarkerStyle(4);
	    marker->SetMarkerColor(2);
	    TPolyMarker* marker2 = new TPolyMarker(1);
	    marker2->SetMarkerSize(0.8);
	    marker2->SetMarkerStyle(29);
	    marker2->SetMarkerColor(2);
	    string cname = "c"+to_string(lprocess+jev);
	    cname += "_"+to_string(channel[iwf]);
	    TCanvas* c = new TCanvas(cname.c_str(), "", 0, 0, 1200, 800);
	    hwf->SetLineColor(4);
	    hwf->SetMarkerColor(4);
	    hwf->Draw();
	    hpz->SetLineColor(2);
	    hpz->SetMarkerColor(2);
	    hpz->Draw("same");
	    hft->SetLineColor(6);
	    hft->SetMarkerColor(6);
	    hft->Draw("same c");
	    hst->SetLineColor(8);
	    hst->SetMarkerColor(8);
	    hst->Draw("same");
	    hat->SetLineColor(7);
	    hat->SetMarkerColor(7);
	    hat->Draw("same");
	    hct->SetLineColor(1);
	    hct->SetMarkerColor(1);
	    hct->Draw("same");
	    int foff =
	      (int)(twf[ich]["fast_trap"]->GetParam("fast_flat") +
		    twf[ich]["fast_trap"]->GetParam("fast_fall"))/sampling;
	    marker2->SetPoint(0,t0[iwf]*sampling,
			      hft->GetBinContent(1+(int)(t0[iwf]-foff)));
	    marker->SetPoint(0, t1[iwf]*sampling,
			     hft->GetBinContent(1+(int)(t1[iwf]-foff)));
	    marker->SetPoint(1,t10[iwf]*sampling,
			     hft->GetBinContent(1+(int)(t10[iwf]-foff)));
	    marker->SetPoint(2,t50[iwf]*sampling,
			     hft->GetBinContent(1+(int)(t50[iwf]-foff)));
	    marker->SetPoint(3,t90[iwf]*sampling,
			     hft->GetBinContent(1+(int)(t90[iwf]-foff)));
	    marker->SetPoint(4,t99[iwf]*sampling,
			     hft->GetBinContent(1+(int)(t99[iwf]-foff)));
	    marker->SetPoint(5, hst->GetBinCenter(hst->GetMaximumBin()),
			     trapmax[iwf]);
	    marker->SetPoint(6, pickoff[iwf]*sampling, trappick[iwf]);
	    marker->Draw();
	    marker2->Draw();
	    c->Update();
	    outfile->cd(("ch"+to_string(channel[iwf])).c_str());
	    c->Write();
	    tdir->cd();
	    delete hpz;
	    delete hft;
	    delete hst;
	    delete hat;
	    delete hct;
	    delete marker;
	    delete marker2;
	    delete c;
	  }
	}
	outtree->Fill();
	run = tmprun;
      }
      lprocess += (int) wfcount[mw_index].size();
      wfcount[mw_index].resize(0);
      for(auto& m : twf) for(auto& p : m) if(p.second) delete p.second;
      twf.clear();
      for(auto& vwf : mwf) vwf[mw_index]->ClearWaveforms();
      times[3] += pwatch.RealTime();
      pwatch.Start();
    }
    // set the next thread index
    block_count ++;
    cthread = -1;
    int wcount = -1;
    if(!processing_complete){
      while(cthread < 0){
	wcount ++;
	if(wcount == 1) cout << "waiting for available thread..." << endl;
	for(int ith=0; ith<(int)threads.size(); ith++)
	  if(sthread[ith] == EMPTY || sthread[ith] == WRITTEN){
	    cthread = ith;
	    bthread[cthread] = block_count;
	  sthread[cthread] = POPULATING;
	  if(iev < nentries - 1)
	    cout << "reading waveforms into thread " << ith << endl;
	  read_events = true;
	  break;
	  }
	  else if(sthread[ith] == CPU_PROCESSING ||
		  sthread[ith] == GPU_PROCESSING){
	    if(threads[ith].wait_for(milliseconds(1)) == future_status::ready){
	      cthread = ith;
	      sthread[ith] = COMPLETE;
	      read_events = false;
	      block_count --;
	      break;
	    }
	  }
      }
    }
    times[5] += pwatch.RealTime();
    pwatch.Start();
  }

  // write outputs
  outfile->cd();
  outtree->AutoSave();
  for(int i=0; i<(int)chan_map.size(); i++){
    if(hdeltat[i]) hdeltat[i]->Write();
    if(henergy[i]) henergy[i]->Write();
    if(henergyf[i]) henergyf[i]->Write();
    if(henergyc1[i]) henergyc1[i]->Write();
    if(henergyc2[i]) henergyc2[i]->Write();
    if(hbase_energy[i]) hbase_energy[i]->Write(hbase_energy[i]->GetName());
    if(hbase_deltat[i]) hbase_deltat[i]->Write(hbase_deltat[i]->GetName());
    if(hbrms_deltat[i]) hbrms_deltat[i]->Write(hbrms_deltat[i]->GetName());
    if(hbrms_energy[i]) hbrms_energy[i]->Write();
    if(hdecay_energy[i]) hdecay_energy[i]->Write();
    if(hdecay_deltat[i]) hdecay_deltat[i]->Write();
    if(hamp_energy[i]) hamp_energy[i]->Write();
    if(haoe_energy[i]) haoe_energy[i]->Write();
    if(hdcr_energy[i]) hdcr_energy[i]->Write();
    if(hrise_energy[i]) hrise_energy[i]->Write();
    if(hamp_energyf[i]) hamp_energyf[i]->Write();
    if(haoe_energyf[i]) haoe_energyf[i]->Write();
    if(hdcr_energyf[i]) hdcr_energyf[i]->Write();
    if(hrise_energyf[i]) hrise_energyf[i]->Write();
    if(hamp_energyc[i]) hamp_energyc[i]->Write();
    if(haoe_energyc[i]) haoe_energyc[i]->Write();
    if(hdcr_energyc[i]) hdcr_energyc[i]->Write();
    if(hrise_energyc[i]) hrise_energyc[i]->Write();
  }
  tdir->cd();
  int ch_count = 0;
  for(auto const& p : chan_map){
    jvalue["channel_map"][ch_count] = p.first;
    ch_count ++;
  }
  Json::StreamWriterBuilder builder;
  builder["indentation"] = "";
  TObjString* jout=new TObjString(Json::writeString(builder, jvalue).c_str());
  outfile->cd();
  jout->Write("process_config");
  outfile->Close();

  if(skip_chan.size() != 0){
    cout << "did not process channel(s) ";
    for(auto const& ch : skip_chan)
      cout << ch << " ";
    cout << endl;
  }
  
  times[4] += pwatch.RealTime();
  double total_time = twatch.RealTime();
  cout << "Processed " << nentries << " events in " << fixed << setprecision(2)
       << total_time << " seconds: "
       << (int) (nentries/total_time) << " evts/s" << endl;
  cout << "  Setup:           " << fixed << setprecision(2) << times[0] <<endl;
  cout << "  Data in:         " << fixed << setprecision(2) << times[1] <<endl;
  cout << "  CPU Processing:  " << fixed << setprecision(2) << times[2] <<endl;
  cout << "  GPU Processing:  " << fixed << setprecision(2) << times[6] <<endl;
  cout << "  Tree/Histograms: " << fixed << setprecision(2) << times[3] <<endl;
  cout << "  Data out:        " << fixed << setprecision(2) << times[4] <<endl;
  cout << "  Thread mgmt:     " << fixed << setprecision(2) << times[5] <<endl;
  
  return 0;
}
