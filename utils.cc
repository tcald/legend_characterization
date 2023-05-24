#include "utils.hh"
#include <MGWFPoleZeroCorrection.hh>
#include "MGWFPZCorrector.hh"
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
    mwf["pz2_cor"]   = new MultiWaveform(*wf, false);
    mwf["pz_ct"]     = new MultiWaveform(*wf, false);
    mwf["slow_trap"] = new MultiWaveform(*wf, false);
    mwf["fast_trap"] = new MultiWaveform(*wf, false);
    mwf["avse_trap"] = new MultiWaveform(*wf, false);
    mwf["base_sub"]->ClearWaveforms();
    mwf["pz_cor"]->ClearWaveforms();
    mwf["pz2_cor"]->ClearWaveforms();
    mwf["pz_ct"]->ClearWaveforms();
    mwf["slow_trap"]->ClearWaveforms();
    mwf["fast_trap"]->ClearWaveforms();
    mwf["avse_trap"]->ClearWaveforms();
  }
  // mgdo transforms
  MGWFPoleZeroCorrection pz("mgwfpz_"+to_string(index));
  pz.SetDecayConstant(wf->GetParam("pz_decay"));
  pz.SetDecayConstant2(wf->GetParam("pz2_decay"));
  pz.SetDecayFraction2(0.0);
  MGWFPZCorrector pz2("mgwfpz2_"+to_string(index));
  pz2.SetDecay1(exp(-1.0/(wf->GetParam("pz_decay")/16)));
  pz2.SetDecay2(exp(-1.0/(wf->GetParam("pz2_decay")/16)));
  pz2.SetPercentOvershoot(0.0);//wf->GetParam("pz2_frac"));
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
    double base = wf->GetWFParam(iwf, "base");//wf->GetParam("resting_base");
    for_each(v.begin(), v.end(), [&](double& s){ s-=base; });
    w->SetData(v);
    double base_rms = sqrt(inner_product(v.begin()+4, v.begin()+4+nbase,
					 v.begin()+4, 0.0) / nbase);
    wf->SetWFParam(iwf, "base_rms", base_rms);
    // fit the exponential tail
    vector<double> yfit(v.end()-wf->GetParam("nefit_samples")-4, v.end()-4);
    vector<double> xfit(yfit.size());
    iota(xfit.begin(), xfit.end(), 0.0);
    for_each(xfit.begin(), xfit.end(), [&](double& s){ s*=sampling; });
    double sx2y = 0.0, sylny = 0.0, sxy = 0.0, sxylny = 0.0, sy=0.0;
    int nskip = 0;
    for(size_t i=0; i<xfit.size(); i++){
      double lny = 0.0;
      if(yfit[i] > 0.0) lny = log(yfit[i]);
      else{
	nskip ++;
	continue;
      }
      sx2y   += xfit[i] * xfit[i] * yfit[i];
      sylny  += yfit[i] * lny;
      sxy    += xfit[i] * yfit[i];
      sxylny += xfit[i] * yfit[i] * lny;
      sy     += yfit[i];
    }
    double a = exp((sx2y*sylny-sxy*sxylny) / (sy*sx2y-sxy*sxy));
    double tau = -1 / ((sy*sxylny-sxy*sylny) / (sy*sx2y-sxy*sxy));
    double exp_chi2 = 0.0;
    for(size_t i=0; i<xfit.size(); i++){
      if(yfit[i] <= 0) continue;
      double dy = yfit[i] - a * exp(-xfit[i]/tau);
      exp_chi2 += dy * dy;
    }
    exp_chi2 /= (xfit.size() - nskip - 2) * base_rms * base_rms;
    a *= exp((v.size()-wf->GetParam("nefit_samples")-4) * sampling / tau);
    wf->SetWFParam(iwf, "decay_const", tau);
    wf->SetWFParam(iwf, "decay_amp", a);
    wf->SetWFParam(iwf, "decay_chi2", exp_chi2);
    wf->SetWFParam(iwf, "efit_min", (v.size()-xfit.size()-4)*sampling);
    wf->SetWFParam(iwf, "efit_max", (v.size()-4)*sampling);
    // mgdo transforms
    MGTWaveform twfp, twfp2, twfs, twff, twfa, twfc, twfe;
    twff.SetTOffset(wf->GetParam("fast_flat")+wf->GetParam("fast_fall"));
    pz.TransformOutOfPlace(*w, twfp);
    pz.SetDecayFraction2(wf->GetParam("pz2_frac"));
    pz.TransformOutOfPlace(*w, twfp2);
    ct.TransformOutOfPlace(*w, twfc);
    strap.TransformOutOfPlace(twfc, twfe);
    // linear baseline correction
    if((int) wf->GetParam("base_fit") == 1){
      vector<double> vwfp = twfp.GetVectorData();
      int nbase_fit = (int)min((size_t)wf->GetParam("nbase_fit"),vwfp.size());
      double sx = 0.0, sy = 0.0, sxy = 0.0, sx2 = 0.0;
      for(unsigned int x=0; x<(unsigned int)nbase_fit; x++){
	sx  += x;
	sy  += vwfp[x];
	sxy += x*vwfp[x];
	sx2 += x*x;
      }
      double d = nbase_fit * sx2 - sx * sx;
      double b = (sx2 * sy - sx * sxy) / d;
      double m = (nbase_fit * sxy - sx * sy) / d;
      if((int) wf->GetParam("base_cor") == 1){
	for(size_t x=0; x<vwfp.size(); x++) vwfp[x] -= b + m * x;
	twfp.SetData(vwfp);
      }
      double chi2 = 0.0;
      for(size_t x=0; x<vwfp.size(); x++){
	double dy = vwfp[x] - (b + m * x);
	chi2 += dy * dy;
      }
      chi2 /= (nbase_fit - 2) * base_rms * base_rms;
      wf->SetWFParam(iwf, "baseslope", m / sampling);
      wf->SetWFParam(iwf, "baseslope_chi2", chi2);
	
    }
    strap.TransformOutOfPlace(twfp, twfs);
    ftrap.TransformOutOfPlace(twfp, twff);
    atrap.TransformOutOfPlace(twfp, twfa);
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
    double cti = 0.0;
    if(t0+coff*2 < (int) v.size()){
      cti = accumulate(v.begin()+(int)(t0+coff),
		       v.begin()+(int)(t0+coff*2), 0.0);
      cti -= accumulate(v.begin()+(int)(t0), v.begin()+(int)(t0+coff), 0.0);
    }
    wf->SetWFParam(iwf, "ct_integral", cti);
    // dcr
    if(!isnan(t99)){
      int ndcr = (int) wf->GetParam("ndcr_samples");
      int d0 = min((int)t99, (int)v.size()-2*ndcr-1);
      double ds = accumulate(v.end()-ndcr-4, v.end()-4, 0.0);
      ds -= accumulate(v.begin()+d0, v.begin()+d0+ndcr, 0.0);
      ds /= v.size() - d0 - ndcr;
      wf->SetWFParam(iwf, "dcrslope", (float) ds);
    }
    else wf->SetWFParam(iwf, "dcrslope", 0.0);
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
    pickoff = max(0, min(pickoff, (int) v.size()));
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
    // fit the fast exponential decay
    if(pickoff >= v.size()-wf->GetParam("nefit_samples")-4){
      wf->SetWFParam(iwf, "decay2_const", 0.0);
      wf->SetWFParam(iwf, "decay2_amp", 0.0);
      wf->SetWFParam(iwf, "decay2_offset", 0.0);
      wf->SetWFParam(iwf, "decay2_chi2", 0.0);
      wf->SetWFParam(iwf, "fefit_min", pickoff*sampling);
      wf->SetWFParam(iwf, "fefit_max", pickoff*sampling);
      wf->SetWFParam(iwf, "decay2_frac", 0.0);
    }
    else{
      v = twfp.GetVectorData();
      size_t bstart = (size_t) (wf->GetWFParam(iwf, "efit_min") / sampling);
      double pzlate = accumulate(v.begin()+bstart,
				 v.begin()+bstart+nbase, 0.0) / nbase;
      for_each(v.begin(), v.end(), [&](double& s){ s-=pzlate; });
      size_t end = (wf->GetWFParam(iwf, "efit_min")/sampling+pickoff)/2;
      xfit.assign(end-pickoff, 0.0);
      yfit.assign(xfit.size(), 0.0);
      iota(xfit.begin(), xfit.end(), 0.0);
      for_each(xfit.begin(), xfit.end(), [&](double& s){ s*=sampling; });
      for(size_t i=0; i<yfit.size(); i++) yfit[i] = v[pickoff+i];
      sx2y = 0.0; sylny = 0.0; sxy = 0.0; sxylny = 0.0; sy = 0.0;
      nskip = 0;
      for(size_t i=0; i<xfit.size(); i++){
	double lny = 0.0;
	if(yfit[i] > 0.0) lny = log(yfit[i]);
	else{
	  nskip ++;
	  continue;
	}
	sx2y   += xfit[i] * xfit[i] * yfit[i];
	sylny  += yfit[i] * lny;
	sxy    += xfit[i] * yfit[i];
	sxylny += xfit[i] * yfit[i] * lny;
	sy     += yfit[i];
      }
      a = exp((sx2y*sylny-sxy*sxylny) / (sy*sx2y-sxy*sxy));
      tau = -1 / ((sy*sxylny-sxy*sylny) / (sy*sx2y-sxy*sxy));
      exp_chi2 = 0.0;
      for(size_t i=0; i<xfit.size(); i++){
	if(yfit[i] <= 0.0) continue;
	double dy = yfit[i] - a * exp(-xfit[i]/tau);
	exp_chi2 += dy * dy;
      }
      exp_chi2 /= (xfit.size() - nskip - 2) * base_rms * base_rms;
      a *= exp(pickoff * sampling / tau);
      wf->SetWFParam(iwf, "decay2_const", tau);
      wf->SetWFParam(iwf, "decay2_amp", a);
      wf->SetWFParam(iwf, "decay2_offset", pzlate);
      wf->SetWFParam(iwf, "decay2_chi2", exp_chi2);
      wf->SetWFParam(iwf, "fefit_min", pickoff*sampling);
      wf->SetWFParam(iwf, "fefit_max", (pickoff+xfit.size())*sampling);
      double a0 = wf->GetWFParam(iwf, "decay_const") *
	exp(-t99*sampling / wf->GetWFParam(iwf, "decay_amp"));
      double a1 = wf->GetWFParam(iwf, "decay2_const") *
	exp(-t99*sampling / wf->GetWFParam(iwf, "decay2_amp"));
      wf->SetWFParam(iwf, "decay2_frac", a1/(a0+a1));
    }
    // copy the transformed waveforms to the return map
    if(copy_wf){
      mwf["base_sub"]->AddWaveform(w);
      mwf["pz_cor"]->AddWaveform(twfp);
      mwf["pz2_cor"]->AddWaveform(twfp2);
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

TH1D* GetWFHist(MultiWaveform* w, int j, string n){
  MGTWaveform* wf = w->GetWaveform(j);
  vector<double> v = wf->GetVectorData();
  TH1D* h = new TH1D(n.c_str(), "", v.size(), w->GetWFOffset(j),
		     w->GetWFOffset(j)+w->GetParam("sampling")*v.size());
  h->SetXTitle("Time (ns)");
  h->SetYTitle("ADC");
  for(unsigned i=0; i<v.size(); i++) h->SetBinContent(i+1, v[i]);
  delete wf;
  return h;
}
