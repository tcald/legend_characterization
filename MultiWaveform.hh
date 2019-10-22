#ifndef __MULTIWAVEFORM__
#define __MULTIWAVEFORM__

#include <MGTWaveform.hh>
#include <vector>
#include <map>
#ifdef __CUDA
#include <thrust/host_vector.h>
#endif

class MultiWaveform{
public:

  MultiWaveform(int i=0);
  MultiWaveform(MultiWaveform& wf, bool copy);
  virtual ~MultiWaveform();

  virtual int GetID() const;
  virtual void SetID(int i);

  virtual int GetNWaveforms() const;
  virtual MGTWaveform* GetWaveform(int i) const;
  
  virtual void AddWaveform(MGTWaveform& wf);
  virtual void AddWaveform(MGTWaveform* wf);
  virtual void AddMultiWaveform(MultiWaveform wf);
  virtual void ClearWaveforms();

  virtual std::vector<int> GetWFLength() const;
  virtual std::vector<int> GetWFStart() const;
  virtual std::vector<float> GetWFSampling() const;
  virtual std::vector<float> GetWFOffset() const;
  
  virtual int GetWFLength(int i) const;
  virtual int GetWFStart(int i) const;
  virtual int GetWFEnd(int i) const;
  virtual float GetWFSampling(int i) const;
  virtual float GetWFOffset(int i) const;

  virtual std::vector<std::pair<int, int> > GetBlocks(int n) const;
  
  virtual std::map<std::string, float> GetParams() const;
  virtual float GetParam(std::string p) const;
  virtual void  SetParam(std::string p, float v);
  virtual float GetWFParam(int i, std::string p) const;
  virtual void  SetWFParam(int i, std::string p, float v);
  virtual void  SetWFParam(std::string p, std::vector<float> v);
  virtual double GetWFDParam(int i, std::string p) const;
  virtual void   SetWFDParam(int i, std::string p, double d);
  virtual void   SetWFDParam(std::string p, std::vector<double> v);

protected:
  
  int id;
  std::vector<int> length;
  std::vector<int> start;
  std::vector<float> sampling;
  std::vector<float> offset;
  std::map<std::string, float> param;
  std::vector<std::map<std::string, float> > wfparam;
  std::vector<std::map<std::string, double> > wfdparam;

public:
  #ifdef __CUDA
  thrust::host_vector<float> wfy;
  thrust::host_vector<int>   wfi;
  #else
  std::vector<float> wfy;
  std::vector<int>   wfi;
  #endif
};

#endif
