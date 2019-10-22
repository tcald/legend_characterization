#include "MultiWaveform.hh"
#include "utils.hh"
#include <MGTWaveform.hh>
#include <MGWFPoleZeroCorrection.hh>
#include <MGWFTrapezoidalFilter.hh>
#include <MGWFAsymTrapezoidalFilter.hh>
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#ifdef __CUDA
#include "gutils.hh"
#include <thrust/host_vector.h>
#endif

using namespace std;

int main(){

  // root options
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetLabelFont(132, "XYZ");
  gStyle->SetTitleFont(132, "XYZ");

  const bool plot = false;
  const int wf_per_block = 100000;
  const int n_blocks = 2;
  const int prescale = (wf_per_block*n_blocks) / 100;
  const int nsamples = 2048;
  const float sampling = 10.0;
  const float pz_decay = 50000.0;
  const float ct_decay = 75000.0;
  const float exp_decay = 100000.0;
  const float exp_ampl  = 0.1;
  const float exp_base  = 0.5; 
  const float slow_rise = 4000.0;
  const float slow_flat = 2000.0;
  const float fast_rise = 40.0;
  const float fast_flat = 100.0;
  const float fast_fall = 2000.0;
  const float avse_rise = 100.0;
  const float avse_flat = 0.0;
  
  vector<double> times(4, 0.0);
  TStopwatch timer;
  vector<TCanvas*> cwf;

  for(int ib=0; ib<n_blocks; ib++){
    timer.Start();
    // construct a set of waveforms
    MultiWaveform* mwf = new MultiWaveform(0);
    mwf->SetParam("sampling", sampling);
    mwf->SetParam("resting_base", exp_base);
    mwf->SetParam("pz_decay", -sampling/pz_decay);
    mwf->SetParam("ct_decay", -sampling/ct_decay);
    mwf->SetParam("slow_nrise", slow_rise / sampling);
    mwf->SetParam("slow_nflat", slow_flat / sampling);
    mwf->SetParam("fast_nrise", fast_rise / sampling);
    mwf->SetParam("fast_nflat", fast_flat / sampling);
    mwf->SetParam("fast_nfall", fast_fall / sampling);
    mwf->SetParam("avse_nrise", avse_rise / sampling);
    mwf->SetParam("avse_nflat", avse_flat / sampling);
    mwf->SetParam("slow_rise", slow_rise);
    mwf->SetParam("slow_flat", slow_flat);
    mwf->SetParam("fast_rise", fast_rise);
    mwf->SetParam("fast_flat", fast_flat);
    mwf->SetParam("fast_fall", fast_fall);
    mwf->SetParam("avse_rise", avse_rise);
    mwf->SetParam("avse_flat", avse_flat);
    mwf->SetParam("t0_thresh", 2.0);
    mwf->SetParam("nbase_samples", 200.0);
    mwf->SetParam("nefit_samples", 500);
    mwf->SetParam("ndcr_samples", 100.0);
    mwf->SetParam("ct_frac", 1.0e-6);
    mwf->SetParam("ct_offset", 2000.0);
    for(int i=0; i<wf_per_block; i++){
      vector<double> v(nsamples);
      for(int j=0; j<nsamples; j++)
	v[j] = exp_base + exp_ampl * exp(-j*sampling/exp_decay);
      for(int j=nsamples/2; j<nsamples; j++)
	v[j] += (i+1) * exp(-(j-nsamples/2)*sampling/pz_decay);
      MGTWaveform* wf = new MGTWaveform();;
      wf->SetID(0);
      wf->SetSamplingFrequency(0.1);
      wf->SetWFType(MGWaveform::kADC);
      wf->SetData(v);
      mwf->AddWaveform(wf);
      mwf->SetWFParam(i, "base", exp_base);
      mwf->SetWFParam(i, "ampl", exp_ampl);
      mwf->SetWFParam(i, "tau",  -1.0/exp_decay);
      delete wf;
    }
    times[0] += timer.RealTime();
    timer.Start();
    
    // process waveforms on the gpu  
    #ifdef __CUDA
    map<string, MultiWaveform*> twf = ProcessMultiWaveformGPU(mwf, plot);
    times[1] += timer.RealTime();
    timer.Start();
    #else
    cout << "No CUDA support found, using CPU only" << endl;
    #endif

    // process waveforms on the cpu
    map<string, MultiWaveform*> rwf = ProcessMultiWaveform(mwf, plot);
    times[2] += timer.RealTime();
    timer.Start();
      
    // plot the waveforms and their transforms
    if(plot){
      for(int i=0; i<mwf->GetNWaveforms(); i++){
	int iwf = ib*wf_per_block + i;
	if(iwf % prescale != 0) continue;
	MGTWaveform* wf = mwf->GetWaveform(i);
	TH1D* hwf = wf->GimmeUniqueHist();
	hwf->SetLineColor(2);
	hwf->SetTitle("");
	hwf->SetXTitle("Time (ns)");
	delete wf;
	#ifdef __CUDA
	MGTWaveform* wfb = twf["base_sub"]->GetWaveform(i);
	MGTWaveform* wfp = twf["pz_cor"]->GetWaveform(i);
	MGTWaveform* wfc = twf["pz_ct"]->GetWaveform(i);
	MGTWaveform* wfs = twf["slow_trap"]->GetWaveform(i);
	MGTWaveform* wff = twf["fast_trap"]->GetWaveform(i);
	MGTWaveform* wfa = twf["avse_trap"]->GetWaveform(i);
	#else
	MGTWaveform* wfb = rwf["base_sub"]->GetWaveform(i);
	MGTWaveform* wfp = rwf["pz_cor"]->GetWaveform(i);
	MGTWaveform* wfc = rwf["pz_ct"]->GetWaveform(i);
	MGTWaveform* wfs = rwf["slow_trap"]->GetWaveform(i);
	MGTWaveform* wff = rwf["fast_trap"]->GetWaveform(i);
	MGTWaveform* wfa = rwf["avse_trap"]->GetWaveform(i);
	#endif
	TH1D* hwfb = wfb->GimmeUniqueHist();
	hwfb->SetLineColor(4);
	delete wfb;
	TH1D* hwfp = wfp->GimmeUniqueHist();
	hwfp->SetLineColor(1);
	delete wfp;
	TH1D* hwfc = wfc->GimmeUniqueHist();
	hwfp->SetLineColor(40);
	delete wfc;
	TH1D* hwfs = wfs->GimmeUniqueHist();
	hwfs->SetLineColor(6);
	delete wfs;
	TH1D* hwff = wff->GimmeUniqueHist();
	hwff->SetLineColor(8);
	delete wff;
	TH1D* hwfa = wfa->GimmeUniqueHist();
	hwfa->SetLineColor(28);
	delete wfa;
	double minval = min(hwfp->GetMinimum(), hwf->GetMinimum());
	double maxval = max(hwfp->GetMaximum(), hwf->GetMaximum());
	hwf->GetYaxis()->SetRangeUser(minval-0.02*(maxval-minval),
				      maxval+0.01*(maxval-minval));
	TCanvas* c = new TCanvas(("c_"+to_string(iwf)).c_str(), "", 1200, 800);
	c->cd(1);
	hwf->Draw();
	hwfb->Draw("same");
	hwfp->Draw("same");
	hwfc->Draw("same");
	hwfs->Draw("same");
	hwff->Draw("same");
	hwfa->Draw("same");
	cwf.push_back(c);
      }
    }
    times[3] += timer.RealTime();
    timer.Start();

    #ifdef __CUDA
    for(auto const& p : twf) delete p.second;
    #endif
    for(auto const& p : rwf) delete p.second;
  }

  // write plots to file
  TFile* outfile = new TFile("gtest.root", "recreate");
  for(auto const& c : cwf) outfile->WriteTObject(c);
  outfile->Close();
  delete outfile;
  times[4] += timer.RealTime();

  int nwf = wf_per_block * n_blocks;
  cout << fixed << setprecision(4);
  cout << "Waveform generation: " << times[0] << " sec" << endl;
  #ifdef __CUDA
  cout << "GPU processing:      " << times[1] << " sec - " << fixed
       << setprecision(2) << nwf/times[1] << " events/sec" << endl;
  cout << fixed << setprecision(4);
  #endif
  cout << "CPU processing:      " << times[2] << " sec - " << fixed
       << setprecision(2) << nwf/times[2] << " events/sec" << endl;
  cout << fixed << setprecision(4);
  if(plot){
    cout << "ROOT plotting:       " << times[3] << " sec" << endl;
    cout << "File writing:        " << times[4] << " sec" << endl;
  }
    
  return 0;
}
