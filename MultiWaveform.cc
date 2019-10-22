#include "MultiWaveform.hh"
#ifdef __CUDA
#include <thrust/copy.h>
#include <thrust/sequence.h>
#endif

MultiWaveform::MultiWaveform(int i){
  SetID(i);
}

MultiWaveform::MultiWaveform(MultiWaveform& wf, bool copy){
  id       = wf.id;
  param    = wf.param;
  sampling = wf.sampling;
  offset   = wf.offset;
  wfparam  = wf.wfparam;
  wfdparam = wf.wfdparam;
  length = wf.GetWFLength();
  start = wf.GetWFStart();
  if(copy){
    wfy = wf.wfy;
    wfi = wf.wfi;
  }
  else{
    wfy.resize(wf.wfy.size());
    wfi.resize(wf.wfi.size());
  }
}

MultiWaveform::~MultiWaveform(){
  ClearWaveforms();
  wfy.shrink_to_fit();
  wfi.shrink_to_fit();
}

int MultiWaveform::GetID() const{
  return id;
}

void MultiWaveform::SetID(int i){
  id = i;
}

int MultiWaveform::GetNWaveforms() const{
  return (int) length.size();
}

MGTWaveform* MultiWaveform::GetWaveform(int i) const{
  if(i < 0 || i >= GetNWaveforms()) return NULL;
  MGTWaveform* wf = new MGTWaveform();
  wf->SetID(id);
  wf->SetSamplingFrequency(1 / sampling[i]);
  wf->SetWFType(MGWaveform::kADC);
  std::vector<double> v(length[i]);
  for(unsigned j=0; j<v.size(); j++) v[j] = (double) wfy[start[i]+j];
  wf->SetData(v);
  return wf;
}

void MultiWaveform::AddWaveform(MGTWaveform& wf){
  std::vector<double> v = wf.GetVectorData();
  int s = 0;
  if(start.size() > 0) s = start.back() + length.back();
  start.push_back(s);
  length.push_back((int) v.size());
  sampling.push_back((float) wf.GetSamplingPeriod());
  offset.push_back((float) wf.GetTOffset());
  wfparam.push_back(std::map<std::string, float>());
  wfdparam.push_back(std::map<std::string, double>());
  wfy.resize(wfy.size() + v.size(), 0.0);
  wfi.resize(wfi.size() + v.size(), 0.0);
  for(unsigned i=0; i<v.size(); i++){
    int j = start.back() + i;
    wfy[j] = (float) v[i];
    #ifndef __CUDA
    wfi[j] = GetNWaveforms()-1;
    #endif
  }
  #ifdef __CUDA
  thrust::fill(wfi.begin()+start.back(), wfi.end(), GetNWaveforms()-1);
  #endif
}

void MultiWaveform::AddWaveform(MGTWaveform* wf){
  AddWaveform(*wf);
}

void MultiWaveform::AddMultiWaveform(MultiWaveform wf){
  int n = GetNWaveforms();
  int s = start.back() + length.back();
  length.insert(length.end(), wf.length.begin(), wf.length.end());
  start.resize(length.size(), 0);
  for(unsigned i=0; i<wf.start.size(); i++){
    start[n+i] = s;
    s += wf.start[i];
  }
  sampling.insert(sampling.end(), wf.sampling.begin(), wf.sampling.end());
  offset.insert(offset.end(), wf.offset.begin(), wf.offset.end());
  wfparam.insert(wfparam.end(), wf.wfparam.begin(), wf.wfparam.end());
  wfdparam.insert(wfdparam.end(), wf.wfdparam.begin(), wf.wfdparam.end());
  wfy.insert(wfy.end(), wf.wfy.begin(), wf.wfy.end());
  wfi.insert(wfi.end(), wf.wfi.begin(), wf.wfi.end());
}

void MultiWaveform::ClearWaveforms(){
  length.resize(0);
  start.resize(0);
  sampling.resize(0);
  offset.resize(0);
  wfparam.resize(0);
  wfdparam.resize(0);
  wfy.resize(0);
  wfi.resize(0);
}

std::vector<int> MultiWaveform::GetWFLength() const{
  return length;
}

std::vector<int> MultiWaveform::GetWFStart() const{
  return start;
}

std::vector<float> MultiWaveform::GetWFSampling() const{
  return sampling;
}

std::vector<float> MultiWaveform::GetWFOffset() const{
  return offset;
}
 
int MultiWaveform::GetWFLength(int i) const{
  if(i < 0 || i >= GetNWaveforms()) return -1;
  return length[i];
}

int MultiWaveform::GetWFStart(int i) const{
  if(i < 0 || i > GetNWaveforms()) return -1;
  return start[i];
}

int MultiWaveform::GetWFEnd(int i) const{
  return GetWFLength(i) + GetWFStart(i);
}

float MultiWaveform::GetWFSampling(int i) const{
  if(i < 0 || i >= GetNWaveforms()) return 0.0;
  return sampling[i];
}

float MultiWaveform::GetWFOffset(int i) const{
  if(i < 0 || i >= GetNWaveforms()) return 0.0;
  return offset[i];
}

std::vector<std::pair<int, int> > MultiWaveform::GetBlocks(int n) const{
  if(length.size() == 0) return std::vector<std::pair<int, int> >(0);
  std::pair<int, int> cblock(0, 0);
  std::vector<std::pair<int, int> > blocks;
  for(int i=1; i<(int)length.size(); i++){
    if((GetWFEnd(cblock.second)-GetWFStart(cblock.first)) + length[i] > n){
      blocks.push_back(cblock);
      cblock.first = i;
      cblock.second = i;
    }
    cblock.second = i;
  }
  blocks.push_back(cblock);
  return blocks;
}
  
std::map<std::string, float> MultiWaveform::GetParams() const{
  return param;
}

float MultiWaveform::GetParam(std::string p) const{
  if(param.find(p) == param.end()) return 0.0;
  return param.at(p);
}

void MultiWaveform::SetParam(std::string p, float v){
  param[p] = v;
}


float MultiWaveform::GetWFParam(int i, std::string p) const{
  if(i < 0 || i >= GetNWaveforms()) return 0.0;
  if(wfparam[i].find(p) == wfparam[i].end()) return 0.0;
  return wfparam[i].at(p);
}

void MultiWaveform::SetWFParam(int i, std::string p, float v){
  if(i < 0 || i >= GetNWaveforms()) return;
  wfparam[i][p] = v;
}

void MultiWaveform::SetWFParam(std::string p, std::vector<float> v){
  for(int i=0; i<(int) std::min(v.size(), length.size()); i++)
    SetWFParam(i, p, v[i]);
}

double MultiWaveform::GetWFDParam(int i, std::string p) const{
  if(i < 0 || i >= GetNWaveforms()) return 0.0;
  if(wfdparam[i].find(p) == wfdparam[i].end()) return 0.0;
  return wfdparam[i].at(p);
}

void MultiWaveform::SetWFDParam(int i, std::string p, double v){
  if(i < 0 || i >= GetNWaveforms()) return;
  wfdparam[i][p] = v;
}

void MultiWaveform::SetWFDParam(std::string p, std::vector<double> v){
  for(int i=0; i<(int) std::min(v.size(), length.size()); i++)
    SetWFDParam(i, p, v[i]);
}
