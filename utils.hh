#include <MGTWaveform.hh>
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <json/value.h>

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

void PZDiff(const MGTWaveform& input,
	    MGTWaveform& output, double tau, bool init=true){
  double base = 0.0;
  for(int i=0; i<100; i++) base += input.At(input.GetLength()-i-1);
  base /= 100.0;
  if(init) output.MakeSimilarTo(input);
  output.SetLength(input.GetLength());
  double decay = exp(-1.0/tau/input.GetSamplingFrequency());
  output[0] = input.At(0);
  bool start = false;
  for(int i=1; i<(int)input.GetLength(); i++){
    if(input.At(i) > base) start = true;
    if(start)
      output[i] = output[i-1] + (input.At(i)-base) - decay*(input.At(i-1)-base);
    else output[i] = input[i];
  }
}
  
void OvershootCorrection(const MGTWaveform& input,
			 MGTWaveform& output, double amp, double tau,
			 bool init=true){
  if(init) output.MakeSimilarTo(input);
  output.SetLength(input.GetLength());
  double decay = tau * input.GetSamplingFrequency();
  output[0] = input.At(0);
  output[1] = input.At(1);
  for(int i=2; i<(int)input.GetLength(); i++){
    output[i] = input.At(i);
    for(int j=1; j<i; j++)
      output[i] -= amp * (output[j]-output[j-1]) * exp(-(i-j)/decay);
  }
}

bool SetJson(const Json::Value value, const string name, double& param){
  if(value.isMember(name)){
    param = value[name].asDouble();
    cout<<"setting "<<name<<" to "<<param<<endl;
    return true;
  }
  else{
    cout<<"parameter "<<name<<" not found, defaulting to "<<param<<endl;
    return false;
  }
}

bool SetJson(const Json::Value value, const string name, int& param){
  if(value.isMember(name)){
    param = value[name].asInt();
    cout<<"setting "<<name<<" to "<<param<<endl;
    return true;
  }
  else{
    cout<<"parameter "<<name<<" not found, defaulting to "<<param<<endl;
    return false;
  }
}

bool SetJson(const Json::Value value, const string name, string& param){
  if(value.isMember(name)){
    param = value[name].asString();
    cout<<"setting "<<name<<" to "<<param<<endl;
    return true;
  }
  else{
    cout<<"parameter "<<name<<" not found, defaulting to "<<param<<endl;
    return false;
  }
}

bool SetJson(const Json::Value value, const string name,vector<double>& param){
  if(value.isMember(name)){
    param.resize(value[name].size());
    for(int i=0; i<(int)value[name].size(); i++)
      param[i] = value[name][i].asDouble();
    return true;
  }
  else{
    cout <<"parameter "<<name<<" not found"<<endl;
    return false;
  }
}

bool SetJson(const Json::Value value, const string name,vector<int>& param){
  if(value.isMember(name)){
    param.resize(value[name].size());
    for(int i=0; i<(int)value[name].size(); i++)
      param[i] = value[name][i].asInt();
    return true;
  }
  else{
    cout <<"parameter "<<name<<" not found"<<endl;
    return false;
  }
}

bool SetJson(const Json::Value value, const string name,vector<string>& param){
  if(value.isMember(name)){
    param.resize(value[name].size());
    for(int i=0; i<(int)value[name].size(); i++)
      param[i] = value[name][i].asString();
    return true;
  }
  else{
    cout <<"parameter "<<name<<" not found"<<endl;
    return false;
  }
}

