#ifndef __LEGEND_CHAR_UTILS__
#define __LEGEND_CHAR_UTILS__

#include "MultiWaveform.hh"
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <json/value.h>
#include <vector>
#include <map>

using namespace std;

TH1D* Get1DHistOrSum(TFile* file, string hname, TH1D* hist);
TH2D* Get2DHistOrSum(TFile* file, string hname, TH2D* hist);

bool SetJson(const Json::Value value, const string name,
	     double& param, bool verbose=true);
bool SetJson(const Json::Value value, const string name,
	     int& param, bool verbose=true);
bool SetJson(const Json::Value value, const string name,
	     string& param, bool verbose=true);
bool SetJson(const Json::Value value, const string name,
	     vector<double>& param, bool verbose=true);
bool SetJson(const Json::Value value, const string name,
	     vector<int>& param, bool verbose=true);
bool SetJson(const Json::Value value, const string name,
	     vector<string>& param, bool verbose=true);

double GetJsonD(const Json::Value value, const string name, bool verbose=true);
int    GetJsonI(const Json::Value value, const string name, bool verbose=true);
string GetJsonS(const Json::Value value, const string name, bool verbose=true);

map<string, MultiWaveform*> ProcessMultiWaveform(MultiWaveform* wf,
						 bool copy_wf, int index=0);

TH1D* GetWFHist(MultiWaveform* w, int i, string n);

#endif
