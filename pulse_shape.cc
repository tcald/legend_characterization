#include "utils.hh"
#include <MGTEvent.hh>
#include <MGTWaveform.hh>
#include <MGWFPoleZeroCorrection.hh>
#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
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

using namespace std;

TDirectory* tdir;

int main(int argc, char* argv[]){

  // pulse shape parameters
  double min_pzd = 4000.0;
  double max_pzd = 6000.0;
  int steps_pzd = 10;
  double min_osa = 0.2;
  double max_osa = 0.4;
  int steps_osa = 10;
  double min_osd = 50.0;
  double max_osd = 200.0;
  int steps_osd = 15;

  double baserms = 1.82;
  int nbase_samples = 200;
  int nslope_samples = 5;

  // root setup
  gStyle->SetOptStat(0);
  gStyle->SetLabelFont(132, "XYZ");
  gStyle->SetTitleFont(132, "XYZ");
  tdir = gROOT->CurrentDirectory();

  // option handling
  map<int, int> chan_map;
  string base_dir = "";
  int max_wf = 0;
  string outfname = "";
  string fnamebase = "";
  int run_start = -1;
  int run_end = -1;
  string infname = "";
  //int write_wf = 0;
  double rc_decay = 0.0;
  int update_percentage = 5;
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
      //{"writewf",    required_argument, NULL, 'w'},
    {"decayconst", required_argument, NULL, 'D'},
    {"update",     required_argument, NULL, 'u'}
  };
  int opt = getopt_long(argc, argv, "hd:c:n:o:f:r:R:i:w:D:u:", opts, NULL);
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
      cout << "  -w plot every n'th wf"              << endl;
      cout << "  -D decay constant in ns"            << endl;
      cout << "  -u update interval (percentage)"    << endl;
      return 0;
    case 'd': base_dir  = string(optarg); break;
    case 'c': chan_map[atoi(optarg)] = (int) chan_map.size()-1; break;
    case 'n': max_wf    = atoi(optarg);   break;
    case 'o': outfname  = string(optarg); break;
    case 'f': fnamebase = string(optarg); break;
    case 'r': run_start = atoi(optarg);   break;
    case 'R': run_end   = atoi(optarg);   break;
    case 'i': infname   = string(optarg); break;
      //case 'w': write_wf  = atoi(optarg);   break;
    case 'D': rc_decay  = atof(optarg);   break;
    case 'u': update_percentage = atoi(optarg); break;
    default: return 1;
    }
    opt = getopt_long(argc, argv, "hd:c:n:o:f:r:R:i:w:D:u:", opts, NULL);
  }

  // compute parameters and initialize histograms
  MGWFPoleZeroCorrection* pole_zero = new MGWFPoleZeroCorrection();
  pole_zero->SetDecayConstant(rc_decay);
  vector<double> pzd_params(steps_pzd, 0.0);
  vector<double> osa_params(steps_osa, 0.0);
  vector<double> osd_params(steps_osd, 0.0);
  vector<vector<vector<TH1D*> > >
    hist(steps_pzd, vector<vector<TH1D*> >(steps_osa,
					   vector<TH1D*>(steps_osd, NULL)));
  for(int ipzd=0; ipzd<steps_pzd; ipzd++){
    pzd_params[ipzd] = min_pzd + ipzd*(max_pzd-min_pzd)/steps_pzd;
    for(int iosa=0; iosa<steps_osa; iosa++){
      if(ipzd == 0)
	osa_params[iosa] = min_osa + iosa*(max_osa-min_osa)/steps_osa;
      for(int iosd=0; iosd<steps_osd; iosd++){
	if(ipzd == 0 && iosa == 0)
	 osd_params[iosd] = min_osd + iosd*(max_osd-min_osd)/steps_osd;
	string s="h_"+to_string(ipzd)+"_"+to_string(iosa)+"_"+to_string(iosd);
	hist[ipzd][iosa][iosd] = new TH1D(s.c_str(), "", 200, -100., 100.);
	
      }
    }
  }
  vector<TH2D*> hslope(steps_pzd, NULL);
  double osa_step = (max_osa-min_osa)/steps_osa;
  double osd_step = (max_osd-min_osd)/steps_osd;
  for(int ipzd=0; ipzd<steps_pzd; ipzd++)
    hslope[ipzd] = new TH2D(("hslope_"+to_string(ipzd)).c_str(), "",
			    steps_osa, osa_params.front()-0.5*osa_step,
			    osa_params.back()+0.5*osa_step, steps_osd,
			    osd_params.front()-0.5*osd_step,
			    osd_params.back()+0.5*osd_step);
  
  // input reader
  TChain tree("MGTree", "MGTree");                                  
  string fname = base_dir + "/" + fnamebase;                
  if(infname != "") tree.Add(infname.c_str());
  else
    for(int irun=run_start; irun<=run_end; irun++)
      tree.Add((fname + to_string(irun) + ".root").c_str());
  //int nentries = (int) tree.GetEntries();
  TTreeReader reader(&tree);
  TTreeReaderValue<MGTEvent> event(reader, "event");
  
  // read events and compute tail slope for all parameter combinations
  int iev = -1;
  int cpercent = -1;
  int lpercent = -1;
  int nvalid = 0;
  MGTWaveform wf0;
  MGTWaveform wf1;
  while(reader.Next()){
    iev ++;
    cpercent = (int) (100.0*nvalid/max_wf);
    if(cpercent % update_percentage == 0 && cpercent != lpercent)
      cout << cpercent << "%" << endl;
    lpercent = cpercent;
    if(nvalid >= max_wf) break;
    TClonesArray* wfs = event->GetWaveforms();
    //int nwf = (int) chan_map.size();
    for(int iwf=0; iwf<(int)wfs->GetEntriesFast(); iwf++){
      if(chan_map.find(iwf) == chan_map.end()) continue;
      MGTWaveform* wf = (MGTWaveform*) wfs->At(iwf);
      vector<double> vwf = wf->GetVectorData();
      double base = accumulate(vwf.begin(),
			       vwf.begin()+nbase_samples, 0) / nbase_samples;
      for_each(vwf.begin(), vwf.end(), [&](double& s){s-=base;});
      double rms = sqrt(inner_product(vwf.begin(), vwf.begin()+nbase_samples,
				      vwf.begin(), 0.0) / nbase_samples);
      if(rms > 2*baserms) continue;
      vector<double>::iterator itmax = max_element(vwf.begin(), vwf.end());
      if(distance(vwf.begin(), itmax) < 4*nslope_samples ||
	 distance(itmax, vwf.end())   < 4*nslope_samples) continue;
      if(*itmax < 2000) continue;
      int ms = distance(vwf.begin(), itmax);
      wf->SetData(vwf);
      if(nvalid == 0){
	wf0.MakeSimilarTo(*wf);
	wf1.MakeSimilarTo(*wf);
      }
      nvalid ++;
      pole_zero->TransformInPlace(*wf);
      for(int ipzd=0; ipzd<steps_pzd; ipzd++){
	PZDiff(*wf, wf0, pzd_params[ipzd], false);
	for(int iosa=0; iosa<steps_osa; iosa++)
	  for(int iosd=0; iosd<steps_osd; iosd++){
	    OvershootCorrection(wf0, wf1, osa_params[iosa],
				osd_params[iosd], false);
	    vector<double> v = wf1.GetVectorData();
	    double y0 = accumulate(v.begin() + ms, // + nslope_samples,
				   v.begin() + ms + /*2**/nslope_samples, 0.0);
	    double y1 = accumulate(vwf.end() - nslope_samples, vwf.end(), 0.0);
	    hist[ipzd][iosa][iosd]->Fill((y1-y0) / nslope_samples /
					 (vwf.size()-ms-2*nslope_samples));
	  }
      }
    }
  }

  double minslope = 1.0e9;
  double minrms   = 1.0e9;
  int mpzd=-1, mosa=-1, mosd=-1;
  for(int ipzd=0; ipzd<steps_pzd; ipzd++)
    for(int iosa=0; iosa<steps_osa; iosa++)
      for(int iosd=0; iosd<steps_osd; iosd++){
	TH1D* h = hist[ipzd][iosa][iosd];
	if(h->Integral(1, h->GetNbinsX()) < 0.5*nvalid){
	  delete h;
	  continue;
	}
	double mean = h->GetMean();
	double rms  = h->GetRMS();
	hslope[ipzd]->SetBinContent(iosa+1, iosd+1, mean);
	hslope[ipzd]->SetBinError(iosa+1, iosd+1, rms);
	if(abs(mean) < abs(minslope)){// && rms < minrms){
	  minslope = mean;
	  minrms = rms;
	  mpzd = ipzd;
	  mosa = iosa;
	  mosd = iosd;
	}
	delete h;
      }

  cout << "optimal slope of " << minslope << "+/-" << minrms << " at:" << endl;
  cout << "  pzd\t" << mpzd << " - " << pzd_params[mpzd] << endl;
  cout << "  osa\t" << mosa << " - " << osa_params[mosa] << endl;
  cout << "  osd\t" << mosd << " - " << osd_params[mosd] << endl;

  TFile* outfile = new TFile(outfname.c_str(), "recreate");
  outfile->cd();
  for(auto const& h : hslope) h->Write();
  outfile->Close();
  
  return 0;
}
			
