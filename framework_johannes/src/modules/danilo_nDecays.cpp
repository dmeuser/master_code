//Investigating number of Neutralino decays in GGM sample
#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGraph2D.h>
#include <TTreeReader.h>
#include <cmath>

Config const &cfg=Config::get();

extern "C"
void run()
{     
	//~ TFile file("/user/dmeuser/master/data/v03D/GGM_GravitinoLSP_M1-200to1500_M2-200to1500.root","read");
	TFile file("/user/dmeuser/master/data/v03D/GGM_GravitinoLSP_M1-50to1500_M3-1000to2500.root","read");
	
	TTree *tree=(TTree*)file.Get(cfg.treeName);
	UShort_t signal_nNeutralinoDecays = 0;
	UShort_t signal_nBinos = 0;
	UShort_t signal_m1 = 0;
	UShort_t signal_m2 = 0;
	tree->SetBranchAddress("signal_nNeutralinoDecays", &signal_nNeutralinoDecays);
	tree->SetBranchAddress("signal_nBinos", &signal_nBinos);
	tree->SetBranchAddress("signal_m1", &signal_m1);
	tree->SetBranchAddress("signal_m2", &signal_m2);
	
	//~ TH2F histDecays("",";M1 (GeV);M2 (GeV);number of neutralino decays",30,25,1525,30,25,1525);
	TH2F histDecays("",";M1 (GeV);M3 (GeV);number of neutralino decays",30,25,1525,31,975,2525);
	//~ TH2F histBinos("",";M1 (GeV);M2 (GeV);number of neutralino_2",30,25,1525,30,25,1525);
	TH2F histBinos("",";M1 (GeV);M3 (GeV);number of neutralino_2",30,25,1525,31,975,2525);
	//~ TH2F nEvents("",";M1 (GeV);M2 (GeV);number of events",30,25,1525,30,25,1525);
	TH2F nEvents("",";M1 (GeV);M3 (GeV);number of events",30,25,1525,31,975,2525);
	
	Long64_t iEvents = tree->GetEntries();
	int processEvents=cfg.processFraction*iEvents;
	for (int iEvent=0; iEvent<iEvents; iEvent++){
		if (iEvent>processEvents) break;
		if (iEvent%(iEvents/100)==0) {io::log*"."; io::log.flush(); };
		tree->GetEvent(iEvent);
		
		histDecays.Fill(signal_m1,signal_m2,signal_nNeutralinoDecays);
		histBinos.Fill(signal_m1,signal_m2,signal_nBinos);
		nEvents.Fill(signal_m1,signal_m2);
		//~ std::cout<<signal_nNeutralinoDecays<<std::endl;
		
	}
	
	histDecays.Divide(&nEvents);
	histBinos.Divide(&nEvents);
	
	TFile out("../output/stuff/GGM_nNeutralinoDecays.root","update");
	//~ histDecays.Write(TString::Format("nDecays%.1f",cfg.processFraction*100),TObject::kOverwrite);
	histDecays.Write(TString::Format("nDecays_M1M3_%.1f",cfg.processFraction*100),TObject::kOverwrite);
	//~ histBinos.Write(TString::Format("nBinos%.1f",cfg.processFraction*100),TObject::kOverwrite);
	histBinos.Write(TString::Format("nBinos_M1M3_%.1f",cfg.processFraction*100),TObject::kOverwrite);
	out.Close();
	file.Close();	
}
