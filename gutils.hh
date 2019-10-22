#ifndef __TRANSFORMS_CUDA_KERNEL__
#define __TRANSFORMS_CUDA_KERNEL__

#include "MultiWaveform.hh"
#include <vector>
#include <map>
#include <string>
#include <cuda.h>
#include <thrust/host_vector.h>
#include <Gpufit/gpufit.h>

// process a multiwaveform on the gpu
// if copy_wf, then copy the waveforms back to the host in the return value
std::map<std::string, MultiWaveform*>
ProcessMultiWaveformGPU(MultiWaveform* wf, bool copy_wf);

// fit a series of exponentials on the gpu
// if fdecay, then allow the decay constant to float
double FitExponentials(int np,
		       std::vector<REAL>& x, std::vector<REAL>& y,
		       std::vector<REAL>& p,
		       std::vector<REAL>& param, std::vector<int>& state,
		       std::vector<REAL>& chi2, std::vector<int>& iterations,
		       bool fdecay=true);

size_t GPUMemory();

#endif
