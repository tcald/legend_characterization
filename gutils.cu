#include "gutils.hh"
#include <TStopwatch.h>
#include <iostream>
#include <iomanip>
#include <functional>
#include <numeric>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/discard_iterator.h>
#include <thrust/tuple.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/reduce.h>
#include <Gpufit/gpufit.h>

using namespace thrust;

struct exp_prod_functor : public binary_function<float, float, float>{
  __host__ __device__ float operator()(const float& x, const float& y) const{
    return expf(x*y);
  }
};

struct square_sum_functor : public binary_function<float, float, float>{
  __host__ __device__ float operator()(const float& x, const float& y) const{
    return x*x+y*y;
  }
};

struct prod_log_functor : public binary_function<float, float, float>{
  __host__ __device__ float operator()(const float& x, const float& y) const{
    return x*logf(y);
  }
};

struct log_sum_functor : public binary_function<float, float, float>{
  __host__ __device__ float operator()(const float& x, const float& y) const{
    return logf(x)+logf(y);
  }
};

struct saxpy_functor{
  const float a;
  saxpy_functor(float _a) : a(_a) {}
  __host__ __device__ float operator()(const float& x, const float& y){
    return x+a*y;
  }
};

struct sqrt_norm_functor{
  const float a;
  sqrt_norm_functor(float _a) : a(_a) {}
  __host__ __device__ float operator()(const float& x){
    return sqrtf(x/a);
  }
};

// maxima of waveforms in x with keys k
std::vector<float> Maxima(const device_vector<float>& x,
			  const device_vector<int>& k,
			  device_vector<float>& y, device_vector<int>& s){
  pair<device_vector<int>::iterator, device_vector<float>::iterator> p =
    reduce_by_key(k.begin(), k.end(), x.begin(), s.begin(), y.begin(),
		  equal_to<int>(), maximum<float>());
  std::vector<float> v(p.first-s.begin());
  thrust::copy(y.begin(), y.begin()+v.size(), v.begin());
  return v;
}

// maxima in waveforms in x with keys k along with the index of the maxima
// pass the number of waveforms n to avoid having to count keys
std::pair<std::vector<float>,
	  std::vector<int> > MaxWithIndex(const device_vector<float>& x,
					  const device_vector<int>& k,
					  const device_vector<int>& t,
					  device_vector<float>& y,
					  device_vector<int>& s, int n){
  zip_iterator<tuple<device_vector<const float>::iterator,
		     device_vector<const int>::iterator> >
    z0 = make_zip_iterator(make_tuple(x.begin(), t.begin()));
  zip_iterator<tuple<device_vector<float>::iterator,
		     device_vector<int>::iterator> > 
    z1 = make_zip_iterator(make_tuple(y.begin(), s.begin()));
  reduce_by_key(k.begin(), k.end(), z0, make_discard_iterator(), z1,
		equal_to<int>(), maximum<tuple<float, int> >());
  // fixme: the lines below aren't actually accessing the reduced values
  std::vector<float> v(n);
  thrust::copy(y.begin(), y.begin()+n, v.begin());
  std::vector<int>   w(n);
  thrust::copy(s.begin(), s.begin()+n, w.begin());
  return std::make_pair(v, w);
}

// minima of waveforms in x with keys k
std::vector<float> Minima(const device_vector<float>& x,
			  const device_vector<int>& k,
			  device_vector<float>& y, device_vector<int>& s){
  pair<device_vector<int>::iterator, device_vector<float>::iterator> p =
    reduce_by_key(k.begin(), k.end(), x.begin(), s.begin(), y.begin(),
		  equal_to<int>(), minimum<float>());
  std::vector<float> v(p.first-s.begin());
  thrust::copy(y.begin(), y.begin()+v.size(), v.begin());
  return v;
}

// minima in waveforms in x with keys k along with the index of the minima
// pass the number of waveforms n to avoid having to count keys
std::pair<std::vector<float>,
	  std::vector<int> > MinWithIndex(const device_vector<float>& x,
					  const device_vector<int>& k,
					  const device_vector<int>& t,
					  device_vector<float>& y,
					  device_vector<int>& s, int n){
  zip_iterator<tuple<device_vector<const float>::iterator,
		     device_vector<const int>::iterator> >
    z0 = make_zip_iterator(make_tuple(x.begin(), t.begin()));
  zip_iterator<tuple<device_vector<float>::iterator,
		     device_vector<int>::iterator> >
    z1 = make_zip_iterator(make_tuple(y.begin(), s.begin()));
  reduce_by_key(k.begin(), k.end(), z0, make_discard_iterator(), z1,
		equal_to<int>(), minimum<tuple<float, int> >());
  std::vector<float> v(n);
  thrust::copy(y.begin(), y.begin()+n, v.begin());
  std::vector<int>   w(n);
  thrust::copy(s.begin(), s.begin()+n, w.begin());
  return std::make_pair(v, w);
}
  
// integral of values in x for keys k from n to n+m in each wf starting at s
std::vector<float> MeanVals(const device_vector<float>& x,
			    const device_vector<int>& k,
			    device_vector<float>& y, int n, int m,
			    std::vector<int>& s, bool do_sum=true){
  if(do_sum) inclusive_scan_by_key(k.begin(), k.end(), x.begin(), y.begin());
  std::vector<float> v(s.size());
  std::vector<float> tmp(y.size());
  thrust::copy(y.begin(), y.end(), tmp.begin());
  for(unsigned i=0; i<v.size(); i++)
    v[i] = (tmp[s[i]+m+n] - tmp[s[i]+m])/n;
  return v;
}

// integral of values in x for keys k from n to n+m in each wf starting at s
std::vector<float> MeanVals(const device_vector<float>& x,
			    const device_vector<int>& k,
			    device_vector<float>& y,
			    std::vector<int>& n, std::vector<int>& m,
			    std::vector<int>& s, bool do_sum=true){
  if(do_sum) inclusive_scan_by_key(k.begin(), k.end(), x.begin(), y.begin());
  std::vector<float> v(s.size());
  std::vector<float> tmp(y.size());
  thrust::copy(y.begin(), y.end(), tmp.begin());
  for(unsigned i=0; i<v.size(); i++)
    v[i] = (tmp[s[i]+m[i]+n[i]] - tmp[s[i]+m[i]]) / n[i];
  return v;
}

// rms of the values in x for keys k from m to m+n in each wf starting at s
std::vector<float> RMSVals(const device_vector<float>& x,
			   const device_vector<int>& k,
			   device_vector<float>& y, int m, int n,
			   std::vector<int>& s, bool do_sum=true){
  if(do_sum) inclusive_scan_by_key(k.begin(), k.end(), x.begin(), y.begin(),
				   equal_to<int>(), square_sum_functor());
  std::vector<float> v(s.size());
  std::vector<float> tmp(y.size());
  thrust::copy(y.begin(), y.end(), tmp.begin());
  for(unsigned i=0; i<v.size(); i++)
    v[i] = sqrt((tmp[s[i]+m+n] - tmp[s[i]+m])/n);
  return v;
}

std::vector<float> ExpDecay(const device_vector<float>& x,
			    const device_vector<float>& y,
			    const device_vector<int>& k,
			    device_vector<float>& z, std::vector<int>& s){
  const int N = (int) x.size();
  const int n = (int) x.size() / s.size();
  std::vector<float> mean(s.size()), tau(s.size()), tmp(s.size());
  inclusive_scan_by_key(k.begin(), k.begin()+N, x.begin(), z.begin(),
			equal_to<int>(), plus<float>());
  thrust::copy(z.begin(), z.begin()+mean.size(), mean.begin());
  std::for_each(mean.begin(), mean.end(), [&](float& v){ v/=n; });
  transform(x.begin(), x.end(), y.begin(), z.begin(), prod_log_functor());
  inclusive_scan_by_key(k.begin(), k.begin()+N, z.begin(), z.begin(),
			equal_to<int>(), plus<float>());
  thrust::copy(z.begin(), z.begin()+tau.size(), tau.begin());
  inclusive_scan_by_key(k.begin(), k.begin()+N, y.begin(), z.begin(),
			equal_to<int>(), log_sum_functor());
  thrust::copy(z.begin(), z.begin()+tmp.size(), tmp.begin());
 for(unsigned i=0; i<tau.size(); i++) tau[i] = tau[i] -  mean[i]*tmp[i];
  inclusive_scan_by_key(k.begin(), k.begin()+N, x.begin(), z.begin(),
			equal_to<int>(), square_sum_functor());
  for(unsigned i=0; i<tau.size(); i++)
    tau[i] = tau[i] / (n*mean[i]*mean[i]-tmp[i]);
  return tau;
}
  
// pole zero correction on x with keys k, pz=1-exp(-sampling/tau)
void PoleZero(device_vector<float>& x,
	      const device_vector<int>& k,
	      device_vector<float>& y, float pz, bool in_place){
  inclusive_scan_by_key(k.begin(), k.end(), x.begin(), y.begin());
  transform(y.begin(), y.end(),
	    make_constant_iterator(pz), y.begin(), multiplies<float>());
  if(in_place)
    transform(y.begin(), y.end(), x.begin(), x.begin(), plus<float>());
  else
    transform(y.begin(), y.end(), x.begin(), y.begin(), plus<float>());
}

// trapezoidal filter from x->y with keys k
// l samples fall time, m samples flat top, n samples rise time
// assumed to be symmetric if last argument is not specified
void Trap(const device_vector<float>& x,
	  const device_vector<int>& k,
	  device_vector<float>& y, device_vector<int>& s, int l,int m,int n=0){
  if(l <= 0 || m < 0) return;
  if(n <= 0) n = l;
  copy(x.begin(), x.end(), y.begin());
  if(n != l) transform(y.begin(), y.end(),make_constant_iterator<float>(1.f/n),
		       y.begin(), multiplies<float>());
  std::vector<int> N({n, m+n, l+m+n});
  std::vector<float> f({-1.f/n, -1.f/l, 1.f/l});
  for(unsigned i=0; i<N.size(); i++){
    transform(k.begin(), k.end()-N[i], k.begin()+N[i],
	      s.begin(), equal_to<int>());
    if(n == l){
      if(i < 2) transform_if(y.begin()+N[i], y.end(), x.begin(), s.begin(),
			     y.begin()+N[i], minus<float>(),identity<float>());
      else transform_if(y.begin()+N[i], y.end(), x.begin(), s.begin(),
			y.begin()+N[i], plus<float>(), identity<float>());
    }
    else transform_if(y.begin()+N[i], y.end(), x.begin(), s.begin(),
		      y.begin()+N[i], saxpy_functor(f[i]), identity<float>());
  }
  if(n == l) transform(y.begin(), y.end(),make_constant_iterator<float>(1.f/n),
		       y.begin(), multiplies<float>());
  inclusive_scan_by_key(k.begin(), k.end(), y.begin(), y.begin());
}

std::map<std::string, MultiWaveform*>
ProcessMultiWaveformGPU(MultiWaveform* wf,
			bool copy_wf){
  TStopwatch ptime;
  ptime.Start();
  std::vector<std::pair<int, int> > blocks = wf->GetBlocks(50000*2048);
  std::map<std::string, MultiWaveform*> mwf;
  std::vector<int> wfstart = wf->GetWFStart();
  if(blocks.size() == 0) return mwf;
  float sampling = wf->GetParam("sampling");
  // transformed waveforms
  MultiWaveform *wfb=NULL, *wfp=NULL, *wfc=NULL;
  MultiWaveform *wfs=NULL, *wff=NULL, *wfa=NULL;
  if(copy_wf){
    wfb = new MultiWaveform(*wf, copy_wf);
    wfp = new MultiWaveform(*wf, copy_wf);
    //wfc = new MultiWaveform(*wf, copy_wf);
    //wfs = new MultiWaveform(*wf, copy_wf);
    //wff = new MultiWaveform(*wf, copy_wf);
    wfa = new MultiWaveform(*wf, copy_wf);
  }
  wfc = new MultiWaveform(*wf, copy_wf);
  wfs = new MultiWaveform(*wf, copy_wf);
  wff = new MultiWaveform(*wf, copy_wf);
  // device vectors
  device_vector<float> y;
  device_vector<int>   k;
  device_vector<int>   t;
  device_vector<int>   s;
  device_vector<float> z;
  device_vector<float> v;
  // process blocks of waveforms until the multiwaveform is fully processed
  for(auto const& pr : blocks){
    int start = wf->GetWFStart(pr.first);
    int end   = wf->GetWFEnd(pr.second);
    std::vector<int> wfst(wfstart.begin()+pr.first, wfstart.begin()+pr.second);
    std::for_each(wfst.begin(), wfst.end(),
		  [&](int& s){ s -= *(wfstart.begin()+pr.first); });
    if(end-start != (int) y.size()){
      y.resize(end-start);
      k.resize(y.size());
      t.resize(y.size());
      s.resize(y.size());
      z.resize(y.size());
      v.resize(y.size());
    }
    copy(wf->wfy.begin()+start, wf->wfy.begin()+end, y.begin());
    copy(wf->wfi.begin()+start, wf->wfi.begin()+end, k.begin());
    std::vector<int> mwfindex(end-start);
    for(int iwf=pr.first; iwf<pr.second; iwf++){
      std::vector<int> tmp(wf->GetWFLength(iwf));
      std::iota(tmp.begin(), tmp.end(), 0);
      std::copy(tmp.begin(), tmp.end(),
		mwfindex.begin()+wf->GetWFStart(iwf)-start);
    }
    thrust::copy(mwfindex.begin(), mwfindex.end(), t.begin());
    fill(s.begin(), s.end(), 0);
    fill(z.begin(), z.end(), 0.0);
    fill(v.begin(), v.end(), 0.0);
    // baseline subtraction
    transform(y.begin(), y.end(),
	      make_constant_iterator(wf->GetParam("resting_base")),
	      y.begin(), minus<float>());
    if(wfb) copy(y.begin(), y.end(), wfb->wfy.begin()+start);
    int nbase = wf->GetParam("nbase_samples");
    wf->SetWFParam("base",    MeanVals(y, k, z, nbase, 4, wfst));
    wf->SetWFParam("base_rms", RMSVals(y, k, z, nbase, 4, wfst));
    // fit the pz decay constants
    std::vector<float> xtmp(wf->GetParam("nefit_samples"));
    std::iota(xtmp.begin(), xtmp.end(), 0.0);
    std::for_each(xtmp.begin(), xtmp.end(), [&](float& s){ s*=sampling; });
    std::vector<float> xfit((pr.second-pr.first)*xtmp.size());
    std::vector<float> yfit(xfit.size());
    std::vector<int> kfit(xfit.size());
    for(int iwf=pr.first; iwf<pr.second; iwf++){
      int i = xtmp.size() * (iwf-pr.first);
      int j = xtmp.size() + i - 1;
      std::copy(xtmp.begin(), xtmp.end(), xfit.begin()+i);
      std::copy(wf->wfy.begin()+wf->GetWFEnd(iwf)-xtmp.size()-4,
		wf->wfy.begin()+wf->GetWFEnd(iwf)-4, yfit.begin()+i);
      std::fill(kfit.begin()+i, kfit.begin()+j+1, iwf-pr.first);
    }
    device_vector<float> dxfit(xfit.begin(), xfit.end());
    device_vector<float> dyfit(yfit.begin(), yfit.end());
    thrust::copy(kfit.begin(), kfit.end(), s.begin());
    std::vector<float> decay = ExpDecay(dxfit, dyfit, s, z, wfst);
    for(int iwf=pr.first; iwf<pr.second; iwf++)
      wf->SetWFParam(iwf, "decay_const", 1/decay[iwf-pr.first]);
    Trap(y, k, z, s,
         (int) (wf->GetParam("fast_fall") / sampling),
	 (int) (wf->GetParam("fast_flat") / sampling),
	 (int) (wf->GetParam("fast_rise") / sampling));
    thrust::copy(z.begin(), z.end(), wff->wfy.begin()+start);
    // get time points, requires copy of fast trap to cpu above
    exclusive_scan_by_key(k.begin(), k.end(),
			  make_constant_iterator(1), s.begin());
    std::pair<std::vector<float>, std::vector<int> >
      pm = MaxWithIndex(z, k, t, v, s, wf->GetNWaveforms());
    std::vector<float> vb = MeanVals(z, k, v, nbase, 4, wfst);
    std::vector<float> vr =  RMSVals(z, k, v, nbase, 4, wfst);
    float t0_thresh = wf->GetParam("t0_thresh");
    std::vector<int> t0(vb.size());
    std::for_each(vr.begin(), vr.end(), [&](float& s){ s *= t0_thresh; });
    for(int iwf=pr.first; iwf<pr.second; iwf++){
      int jwf = iwf-pr.first;
      int si = wff->GetWFStart(iwf);
      float t1=0.0, t10=0.0, t50=0.0, t90=0.0, t99=0.0;
      for(int i=pm.second[jwf]+wfst[0]-si; i>=si; i--){
	float d = (wff->wfy[i+1] - wff->wfy[i]) / pm.first[jwf];
	if(d == 0.0) d = 1.0e12;
	float w = wff->wfy[i] / pm.first[jwf];
	if(w < 0.01 && t1  == 0.0) t1  = i-si+(w-0.01)/d;
	if(w < 0.1  && t10 == 0.0) t10 = i-si+(w-0.10)/d;
	if(w < 0.5  && t50 == 0.0) t50 = i-si+(w-0.50)/d;
	if(w < 0.9  && t90 == 0.0) t90 = i-si+(w-0.90)/d;
	if(w < 0.99 && t99 == 0.0) t99 = i-si+(w-0.99)/d;
	if(w < vr[jwf]){ t0[jwf] = i-si; break; }
      }
      wf->SetWFParam(iwf, "t0",  (float) t0[jwf]);
      wf->SetWFParam(iwf, "t1",  t1);
      wf->SetWFParam(iwf, "t10", t10);
      wf->SetWFParam(iwf, "t50", t50);
      wf->SetWFParam(iwf, "t90", t90);
      wf->SetWFParam(iwf, "t99", t99);
    }
    // energy estimation with charge trapping correction - David's method 1
    PoleZero(y, k, z, 1-exp(-sampling/wf->GetParam("ct_decay")), false);
    Trap(z, k, v, s,
	 (int) wf->GetParam("slow_nrise"), (int) wf->GetParam("slow_nflat"));
    copy(v.begin(), v.end(), wfc->wfy.begin()+start);
    // charge trapping corrected fixed time pickoff
    int poff = (int)((wf->GetParam("slow_rise") +
		      0.9*wf->GetParam("slow_flat")) / sampling);
    for(int iwf=pr.first; iwf<pr.second; iwf++){
      wf->SetWFParam(iwf, "pickoff", t0[iwf-pr.first] + poff);
      wf->SetWFParam(iwf, "ct1_trappick",
		    wfc->wfy[wfst[iwf-pr.first] + t0[iwf-pr.first] + poff]);
    }
    // vary the charge trapping integration time for optimization
    if(wf->GetParam("nct_steps") > 0.0){
      std::vector<float> tmp(wfc->wfy.size());
      for(int ict=0; ict<(int)wf->GetParam("nct_steps"); ict++){
	PoleZero(y, k, z,
		 1-exp(-sampling/
		       wf->GetParam("ct_decay_"+std::to_string(ict))), false);
	Trap(z, k, v, s,
	     (int)wf->GetParam("slow_nrise"), (int)wf->GetParam("slow_nflat"));
	thrust::copy(v.begin(), v.end(), tmp.begin());
	for(int iwf=pr.first; iwf<pr.second; iwf++)
	  wf->SetWFParam(iwf, "ct1_trappick_"+std::to_string(ict),
			tmp[wfst[iwf-pr.first] + t0[iwf-pr.first] + poff]);
      }
    }
    // standard pole zero correction
    PoleZero(y, k, z, 1-exp(-sampling/wf->GetParam("pz_decay")), true);
    if(wfp) copy(y.begin(), y.end(), wfp->wfy.begin()+start);
    // get the integral for David's 2nd charge trapping method
    int coff = (int) (wf->GetParam("ct_offset") / sampling);
    std::vector<int> ct0=t0, ctn(t0.size(), coff);
    std::vector<float> cti0 = MeanVals(v, k, z, ctn, ct0, wfst);
    std::for_each(ct0.begin(), ct0.end(), [&](int& s){ s+= coff; });
    std::vector<float> cti1 = MeanVals(v, k, z, ctn, ct0, wfst);
    for(int iwf=pr.first; iwf<pr.second; iwf++)
      wf->SetWFParam(iwf, "ct_integral",
		    (cti1[iwf-pr.first]-cti0[iwf-pr.second])*coff);
    // calculate dcr
    int ndcr = (int) wf->GetParam("ndcr_samples");
    std::vector<int> d0(t0.size()), d1(t0.size()), d2(t0.size(), ndcr);
    for(int iwf=pr.first; iwf<pr.second; iwf++){
      d0[iwf-pr.first] = std::max(std::min((int) wf->GetWFParam(iwf, "t99"),
					   wf->GetWFLength(iwf)-2*ndcr-1), 0);
      d1[iwf-pr.first] = wf->GetWFLength(iwf)-ndcr-4;
    }
    std::vector<float> ds0 = MeanVals(y, k, z, d2, d0, wfst, false);
    std::vector<float> ds1 = MeanVals(y, k, z, d2, d1, wfst, false);
    for(int iwf=pr.first; iwf<pr.second; iwf++)
      wf->SetWFParam(iwf, "dcrslope",
		    (ds1[iwf-pr.first]-ds0[iwf-pr.first]) /
		    (wf->GetWFLength(iwf) - d0[iwf-pr.first] - ndcr));
    // trap filter for energy estimation
    Trap(y, k, z, s,
	 (int) wf->GetParam("slow_nrise"), (int) wf->GetParam("slow_nflat"));
    copy(z.begin(), z.end(), wfs->wfy.begin()+start);
    wf->SetWFParam("trap_max", Maxima(z, k, v, s));
    // fixed time energy pickoff
    for(int iwf=pr.first; iwf<pr.second; iwf++){
      double val = wfs->wfy[wfst[iwf-pr.first] + t0[iwf-pr.first] + poff];
      wf->SetWFParam(iwf, "trappick", val);
      val += wf->GetWFParam(iwf, "ct_integral") * wf->GetParam("ct_frac");
      wf->SetWFParam(iwf, "ct2_trappick", val);
    }
    // triangle filter for maximum current estimation
    Trap(y, k, z, s,
	 (int) wf->GetParam("avse_nrise"), (int) wf->GetParam("avse_nflat"));
    if(wfa) copy(z.begin(), z.end(), wfa->wfy.begin()+start);
    wf->SetWFParam("imax", Maxima(z, k, v, s));
  }
  // if copyting wf back to host, populate return map
  if(copy_wf){
    mwf["base_sub"]  = wfb;
    mwf["pz_cor"]    = wfp;
    mwf["pz_ct"]     = wfc;
    mwf["slow_trap"] = wfs;
    mwf["fast_trap"] = wff;
    mwf["avse_trap"] = wfa;
  }
  else{
    delete wfc;
    delete wfs;
    delete wff;
  }
  wf->SetParam("proc_type", 1.0);
  wf->SetParam("proc_time", (float) ptime.RealTime());
  return mwf;
}

double FitExponentials(int np,
		       std::vector<REAL>& x, std::vector<REAL>& y,
		       std::vector<REAL>& p,
		       std::vector<REAL>& param, std::vector<int>& state,
		       std::vector<REAL>& chi2, std::vector<int>& iterations,
		       bool fdecay){
  assert(x.size() == y.size());
  if(x.size() == 0) return 0.0;
  assert(x.size() % np == 0);
  const size_t nfit = x.size() / np;
  assert(p.size() == 2*nfit);
  TStopwatch timer;
  timer.Start();
  param.resize(nfit*2);
  chi2.resize(nfit);
  state.resize(nfit);
  iterations.resize(nfit);
  std::vector<int> fparam(2, 1);
  if(!fdecay) fparam[1] = 0;
  const int status = gpufit(nfit, np, y.data(), 0, EXP_1D, p.data(), 0.001,
			    1000, fparam.data(), LSE, x.size()*sizeof(REAL),
			    reinterpret_cast<char*>(x.data()),
			    param.data(), state.data(),
			    chi2.data(), iterations.data());
  double t = timer.RealTime();
  if(status != ReturnState::OK)
    std::cout << "gpu fitting error, aborting fits" << std::endl;
  std::vector<int> states(5, 0);
  
  std::cout << "Fit " << nfit << " event(s) in "
	    << std::fixed << std::setprecision(2)
	    << t << " seconds: " << nfit/t << " fits/s" << std::endl;
  
  for(std::vector<int>::iterator it=state.begin(); it!=state.end(); it++)
    states[*it] ++;
  
  if(states[0] / (float) nfit > 0.95)
    std::cout << "  " << std::fixed << std::setprecision(2)
	      << 100.*states[0]/nfit << "% of fits converged" << std::endl;
  else{
    std::cout << "  ratio converged     " << (float)states[0]/nfit <<std::endl;
    std::cout << "  ratio max iteration " << (float)states[1]/nfit <<std::endl;
    std::cout << "  ratio singular hess " << (float)states[2]/nfit <<std::endl;
    std::cout << "  ratio neg curvature " << (float)states[3]/nfit <<std::endl;
    std::cout << "  ratio gpu not read  " << (float)states[4]/nfit <<std::endl;
  }
  
  return t;
}

size_t GPUMemory(){
  size_t fmem, tmem;
  cudaError_t e = cudaMemGetInfo(&fmem, &tmem);
  if(e == cudaSuccess) return fmem;
  else return 0;
}
