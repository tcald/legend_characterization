#include "utils.hh"
#include <MGWFPoleZeroCorrection.hh>
#include <MGWFTrapezoidalFilter.hh>
#include <MGWFAsymTrapezoidalFilter.hh>
#include <TStopwatch.h>
#include <TROOT.h>
#include <iostream>
#include <utility>
#include <algorithm>
#include <numeric>
#include <functional>

using namespace std;

TH1D* Get1DHistOrSum(TFile* file, string hname, TH1D* hist){
  TDirectory* tdir = gROOT->CurrentDirectory();
  file->cd();
  TH1D* h = (TH1D*) file->Get(hname.c_str());
  if(!h) return hist;
  if(hist){
    hist->Add(h);
    return hist;
  }
  else{
    tdir->cd();
    return (TH1D*) h->Clone((string(h->GetName())+"_tmp").c_str());
  }
}

TH2D* Get2DHistOrSum(TFile* file, string hname, TH2D* hist){
  TDirectory* tdir = gROOT->CurrentDirectory();
  file->cd();
  TH2D* h = (TH2D*) file->Get(hname.c_str());
  if(!h) return hist;
  if(hist){
    hist->Add(h);
    return hist;
  }
  else{
    tdir->cd();
    return (TH2D*) h->Clone((string(h->GetName())+"_tmp").c_str());
  }
}

bool SetJson(const Json::Value value, const string name,
	     double& param, bool verbose){
  if(value.isMember(name)){
    param = value[name].asDouble();
    if(verbose) cout<<"setting "<<name<<" to "<<param<<endl;
    return true;
  }
  else{
    if(verbose)
      cout<<"parameter "<<name<<" not found, defaulting to "<<param<<endl;
    return false;
  }
}

bool SetJson(const Json::Value value, const string name,
	     int& param, bool verbose){
  if(value.isMember(name)){
    param = value[name].asInt();
    if(verbose) cout<<"setting "<<name<<" to "<<param<<endl;
    return true;
  }
  else{
    if(verbose)
      cout<<"parameter "<<name<<" not found, defaulting to "<<param<<endl;
    return false;
  }
}

bool SetJson(const Json::Value value, const string name,
	     string& param, bool verbose){
  if(value.isMember(name)){
    param = value[name].asString();
    if(verbose) cout<<"setting "<<name<<" to "<<param<<endl;
    return true;
  }
  else{
    if(verbose)
      cout<<"parameter "<<name<<" not found, defaulting to "<<param<<endl;
    return false;
  }
}

bool SetJson(const Json::Value value, const string name,
	     vector<double>& param, bool verbose){
  if(value.isMember(name)){
    param.resize(value[name].size());
    for(int i=0; i<(int)value[name].size(); i++)
      param[i] = value[name][i].asDouble();
    return true;
  }
  else{
    if(verbose) cout <<"parameter "<<name<<" not found"<<endl;
    return false;
  }
}

bool SetJson(const Json::Value value, const string name,
	     vector<int>& param, bool verbose){
  if(value.isMember(name)){
    param.resize(value[name].size());
    for(int i=0; i<(int)value[name].size(); i++)
      param[i] = value[name][i].asInt();
    return true;
  }
  else{
    if(verbose) cout <<"parameter "<<name<<" not found"<<endl;
    return false;
  }
}

bool SetJson(const Json::Value value, const string name,
	     vector<string>& param, bool verbose){
  if(value.isMember(name)){
    param.resize(value[name].size());
    for(int i=0; i<(int)value[name].size(); i++)
      param[i] = value[name][i].asString();
    return true;
  }
  else{
    if(verbose) cout <<"parameter "<<name<<" not found"<<endl;
    return false;
  }
}

double GetJsonD(const Json::Value value, const string name, bool verbose){
  double ret = 0.0;
  SetJson(value, name, ret, verbose);
  return ret;
}

int GetJsonI(const Json::Value value, const string name, bool verbose){
  int ret = 0;
  SetJson(value, name, ret, verbose);
  return ret;
}

string GetJsonS(const Json::Value value, const string name, bool verbose){
  string ret = "";
  SetJson(value, name, ret, verbose);
  return ret;
}

map<string, MultiWaveform*> ProcessMultiWaveform(MultiWaveform* wf,
						 bool copy_wf, int index){
  TStopwatch ptime;
  ptime.Start();
  // transformed waveforms
  map<string, MultiWaveform*> mwf;
  if(copy_wf){
    mwf["base_sub"]  = new MultiWaveform(*wf, false);
    mwf["pz_cor"]    = new MultiWaveform(*wf, false);
    mwf["pz_ct"]     = new MultiWaveform(*wf, false);
    mwf["slow_trap"] = new MultiWaveform(*wf, false);
    mwf["fast_trap"] = new MultiWaveform(*wf, false);
    mwf["avse_trap"] = new MultiWaveform(*wf, false);
    mwf["base_sub"]->ClearWaveforms();
    mwf["pz_cor"]->ClearWaveforms();
    mwf["pz_ct"]->ClearWaveforms();
    mwf["slow_trap"]->ClearWaveforms();
    mwf["fast_trap"]->ClearWaveforms();
    mwf["avse_trap"]->ClearWaveforms();
  }
  // mgdo transforms
  MGWFPoleZeroCorrection pz("mgwfpz_"+to_string(index));
  pz.SetDecayConstant(wf->GetParam("pz_decay"));
  MGWFPoleZeroCorrection ct("mgwfpz_"+to_string(index));
  ct.SetDecayConstant(wf->GetParam("ct_decay"));
  MGWFTrapezoidalFilter strap(wf->GetParam("slow_rise"),
			      wf->GetParam("slow_flat"), 0.0, 0.0,
			      "mgwftf_slow_"+to_string(index));
  MGWFAsymTrapezoidalFilter ftrap(wf->GetParam("fast_rise"),
				  wf->GetParam("fast_flat"),
				  wf->GetParam("fast_fall"),
				  "mgwftf_fast_"+to_string(index));
  MGWFTrapezoidalFilter atrap(wf->GetParam("avse_rise"),
			      wf->GetParam("avse_flat"), 0.0, 0.0,
			      "mgwftf_avse_"+to_string(index));
  // process each waveform
  double sampling = wf->GetParam("sampling");
  for(int iwf=0; iwf<wf->GetNWaveforms(); iwf++){
    MGTWaveform* w = wf->GetWaveform(iwf);
    vector<double> v = w->GetVectorData();
    int nbase = (int) wf->GetParam("nbase_samples");
    wf->SetWFParam(iwf, "base", accumulate(v.begin()+4,
					   v.begin()+4+nbase, 0.0)/nbase);
    double base = wf->GetParam("resting_base");
    for_each(v.begin(), v.end(), [&](double& s){ s-=base; });
    w->SetData(v);
    wf->SetWFParam(iwf, "base_rms",
		  sqrt(inner_product(v.begin()+4, v.begin()+4+nbase,
				     v.begin()+4, 0.0)/nbase));
    // fit the exponential tail
    vector<double> yfit(v.end()-wf->GetParam("nefit_samples")-4, v.end()-4);
    vector<double> xfit(yfit.size());
    iota(xfit.begin(), xfit.end(), 0.0);
    for_each(xfit.begin(), xfit.end(), [&](double& s){ s*=sampling; });
    double xmean = xfit[xfit.size()/2];
    if(xfit.size() % 2 == 0) xmean = (xmean+xfit[xfit.size()/2-1])/2.;
    double s0 = 0.0;
    double s1 = 0.0;
    double s2 = 0.0;
    for(unsigned i=0; i<xfit.size(); i++){
      double x = xfit[i] - xmean;
      double y = log(yfit[i]);
      s0 += x * y;
      s1 += x * xfit[i];
      s2 += y;
    }
    double tau = -s1 / s0;
    double a = exp(s2/xfit.size() + tau*xmean);
    wf->SetWFParam(iwf, "decay_const", tau);
    wf->SetWFParam(iwf, "decay_amp", a);
    // mgdo transforms
    MGTWaveform twfp, twfs, twff, twfa, twfc, twfe;
    twff.SetTOffset(wf->GetParam("fast_flat")+wf->GetParam("fast_fall"));
    pz.TransformOutOfPlace(*w, twfp);
    ct.TransformOutOfPlace(*w, twfc);
    strap.TransformOutOfPlace(twfp, twfs);
    ftrap.TransformOutOfPlace(twfp, twff);
    atrap.TransformOutOfPlace(twfp, twfa);
    strap.TransformOutOfPlace(twfc, twfe);
    // trapezoidal filter extrema
    v = twfs.GetVectorData();
    wf->SetWFParam(iwf, "trap_max", *max_element(v.begin(), v.end()));
    v = twfa.GetVectorData();
    wf->SetWFParam(iwf, "imax", *max_element(v.begin(), v.end()));
    // time points
    v = twff.GetVectorData();
    vector<double>::iterator it = max_element(v.begin(), v.end());
    double ftbase = accumulate(v. begin()+4, v.begin()+nbase+4, 0.0)/nbase;
    for_each(v.begin(), v.end(), [&](double& s){ s=(s-ftbase)/(*it); } );
    double ftrms = sqrt(inner_product(v.begin()+4, v.begin()+nbase+4,
				      v.begin()+4, 0.) / nbase);
    ftrms *= wf->GetParam("t0_thresh");
    double ftoffset = wf->GetParam("fast_flat") + wf->GetParam("fast_fall");
    ftoffset /= sampling;
    double t0=0.0, t1=0.0, t10=0.0, t50=0.0, t90=0.0, t99=0.0;
    for(int i=distance(v.begin(), it); i>=0; i--){
      double d = v[i+1] - v[i];
      if(d == 0.0) d = 1.0e12;
      if(v[i] < 0.01 && t1  == 0.0) t1  = i+(v[i+1]-0.01)/d + ftoffset;
      if(v[i] < 0.1  && t10 == 0.0) t10 = i+(v[i+1]-0.10)/d + ftoffset;
      if(v[i] < 0.5  && t50 == 0.0) t50 = i+(v[i+1]-0.50)/d + ftoffset;
      if(v[i] < 0.9  && t90 == 0.0) t90 = i+(v[i+1]-0.90)/d + ftoffset;
      if(v[i] < 0.99 && t99 == 0.0) t99 = i+(v[i+1]-0.99)/d + ftoffset;
      if(v[i] < ftrms){ t0 = i + ftoffset; break; }
    }
    wf->SetWFParam(iwf, "t0",  (float) t0);
    wf->SetWFParam(iwf, "t1",  (float) t1);
    wf->SetWFParam(iwf, "t10", (float) t10);
    wf->SetWFParam(iwf, "t50", (float) t50);
    wf->SetWFParam(iwf, "t90", (float) t90);
    wf->SetWFParam(iwf, "t99", (float) t99);
    // integral used in David's "2nd" charge trapping method
    v = twfp.GetVectorData();
    int coff = (int) (wf->GetParam("ct_offset") / sampling);
    double cti = accumulate(v.begin()+(int)(t0+coff),
			    v.begin()+(int)(t0+coff*2), 0.0);
    cti -= accumulate(v.begin()+(int)(t0), v.begin()+(int)(t0+coff), 0.0);
    wf->SetWFParam(iwf, "ct_integral", cti);
    // dcr
    int ndcr = (int) wf->GetParam("ndcr_samples");
    int d0 = min((int)t99, (int)v.size()-2*ndcr-1);
    double ds = accumulate(v.end()-ndcr-4, v.end()-4, 0.0);
    ds -= accumulate(v.begin()+d0, v.begin()+d0+ndcr, 0.0);
    ds /= v.size() - d0 - ndcr;
    wf->SetWFParam(iwf, "dcrslope", (float) ds);
    // fixed time pickoff
    v = twfs.GetVectorData();
    int pickoff = (int)(t0+(wf->GetParam("slow_rise")+0.9*
			    wf->GetParam("slow_flat"))/sampling);
    wf->SetWFParam(iwf, "pickoff", pickoff);
    if(pickoff < 0 || pickoff >= (int) v.size()){
      wf->SetWFParam(iwf, "trappick",     0.0);
      wf->SetWFParam(iwf, "ct1_trappick", 0.0);
      wf->SetWFParam(iwf, "ct2_trappick", 0.0);
    }
    else{
      wf->SetWFParam(iwf, "trappick", (float) v[pickoff]);
      wf->SetWFParam(iwf, "ct1_trappick",(float)twfe.GetVectorData()[pickoff]);
      wf->SetWFParam(iwf, "ct2_trappick",(float)(v[pickoff]+
						 wf->GetParam("ct_frac")*cti));
    }
    // vary the charge trapping integration time for optimization
    if(wf->GetParam("nct_steps") > 0.0)
      for(int ict=0; ict<(int)wf->GetParam("nct_steps"); ict++){
	MGWFPoleZeroCorrection ctpz;
	ctpz.SetDecayConstant(wf->GetParam("ct_decay_"+to_string(ict)));
	MGTWaveform tmp;
	ctpz.TransformOutOfPlace(*w, tmp);
	wf->SetWFParam(iwf, "ct1_trappick_"+to_string(ict),
		      (float)tmp.GetVectorData()[pickoff]);
      }
    // copy the transformed waveforms to the return map
    if(copy_wf){
      mwf["base_sub"]->AddWaveform(w);
      mwf["pz_cor"]->AddWaveform(twfp);
      mwf["pz_ct"]->AddWaveform(twfe);
      mwf["slow_trap"]->AddWaveform(twfs);
      mwf["fast_trap"]->AddWaveform(twff);
      mwf["avse_trap"]->AddWaveform(twfa);
    }
    delete w;
  }
  wf->SetParam("proc_type", 0.0);
  wf->SetParam("proc_time", (float) ptime.RealTime());
  return mwf;
}

TH1D* GetWFHist(MultiWaveform* w, int i, string n){
  MGTWaveform* wf = w->GetWaveform(i);
  vector<double> v = wf->GetVectorData();
  TH1D* h = new TH1D(n.c_str(), "", v.size(), w->GetWFOffset(i),
		     w->GetWFOffset(i)+w->GetParam("sampling")*v.size());
  h->SetXTitle("Time (ns)");
  h->SetYTitle("ADC");
  for(unsigned i=0; i<v.size(); i++) h->SetBinContent(i+1, v[i]);
  delete wf;
  return h;
}
