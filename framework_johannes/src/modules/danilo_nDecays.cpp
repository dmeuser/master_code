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

void nDecays() {
	//~ TFile file("/net/data_cms1b/user/dmeuser/master/data/v03D/GGM_GravitinoLSP_M1-200to1500_M2-200to1500.root","read");
	TFile file("/net/data_cms1b/user/dmeuser/master/data/v03D/GGM_GravitinoLSP_M1-50to1500_M3-1000to2500.root","read");
	
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

void nWBosons() {
	TFile file("/net/data_cms1b/user/dmeuser/master/data/v03D/GGM_GravitinoLSP_M1-200to1500_M2-200to1500.root","read");
	
	TTree *tree=(TTree*)file.Get(cfg.treeName);
	int signal_nWBosons = 0;
	int signal_nWBosons_genCol = 0;
	int signal_nWPBosons = 0;
	int signal_nWMBosons = 0;
	UShort_t signal_m1 = 0;
	UShort_t signal_m2 = 0;
	std::vector<tree::IntermediateGenParticle> *intermediateGenParticles=0;
	std::vector<tree::GenParticle> *genParticles=0;
	tree->SetBranchAddress("signal_m1", &signal_m1);
	tree->SetBranchAddress("signal_m2", &signal_m2);
	tree->SetBranchAddress("intermediateGenParticles", &intermediateGenParticles);
	tree->SetBranchAddress("genParticles", &genParticles); 
	
	TH2F histBosons("",";M1 (GeV);M2 (GeV);number of W bosons",30,25,1525,30,25,1525);
	TH2F histBosonsgenCol("",";M1 (GeV);M2 (GeV);number of W bosons (genCol)",30,25,1525,30,25,1525);
	TH2F histBosonsP("",";M1 (GeV);M2 (GeV);number of pos. W bosons",30,25,1525,30,25,1525);
	TH2F histBosonsM("",";M1 (GeV);M2 (GeV);number of neg. W bosons",30,25,1525,30,25,1525);
	TH2F nEvents("",";M1 (GeV);M2 (GeV);number of events",30,25,1525,30,25,1525);
	
	Long64_t iEvents = tree->GetEntries();
	int processEvents=cfg.processFraction*iEvents;
	for (int iEvent=0; iEvent<iEvents; iEvent++){
		if (iEvent>processEvents) break;
		if (iEvent%(iEvents/100)==0) {io::log*"."; io::log.flush(); };
		tree->GetEvent(iEvent);
		signal_nWBosons = 0;
		signal_nWBosons_genCol = 0;
		signal_nWPBosons = 0;
		signal_nWMBosons = 0;
		
		for (tree::IntermediateGenParticle p :*intermediateGenParticles) {
			if (abs(p.pdgId)==24 and p.daughters.size()==2) signal_nWBosons++;
			if (p.pdgId==24 and p.daughters.size()==2) signal_nWPBosons++;
			if (p.pdgId==-24 and p.daughters.size()==2) signal_nWMBosons++;
		}
		for (tree::GenParticle p :*genParticles) {
			if (abs(p.pdgId)==24) signal_nWBosons_genCol++;
		}
		
		histBosons.Fill(signal_m1,signal_m2,signal_nWBosons);
		histBosonsgenCol.Fill(signal_m1,signal_m2,signal_nWBosons_genCol);
		histBosonsP.Fill(signal_m1,signal_m2,signal_nWPBosons);
		histBosonsM.Fill(signal_m1,signal_m2,signal_nWMBosons);
		nEvents.Fill(signal_m1,signal_m2);
		
	}
	
	histBosons.Divide(&nEvents);
	histBosonsgenCol.Divide(&nEvents);
	histBosonsP.Divide(&nEvents);
	histBosonsM.Divide(&nEvents);
	
	TFile out("../output/stuff/GGM_nWBosons.root","update");
	histBosons.Write(TString::Format("nWBosons%.1f",cfg.processFraction*100),TObject::kOverwrite);
	histBosonsgenCol.Write(TString::Format("nWBosonsGenCol%.1f",cfg.processFraction*100),TObject::kOverwrite);
	histBosonsP.Write(TString::Format("nWPBosons%.1f",cfg.processFraction*100),TObject::kOverwrite);
	histBosonsM.Write(TString::Format("nWMBosons%.1f",cfg.processFraction*100),TObject::kOverwrite);
	out.Close();
	file.Close();
}

void dominatProcesses() {
	TFile file("/net/data_cms1b/user/dmeuser/master/data_knut/v26/GGM_GravitinoLSP_M1-200to1500_M2-200to1500_nTuple.root","read");
	
	TTree *tree=(TTree*)file.Get(cfg.treeName);
	int signal_nWBosons = 0;
	UShort_t signal_m1 = 0;
	UShort_t signal_m2 = 0;
	std::vector<tree::IntermediateGenParticle> *intermediateGenParticles=0;
	std::vector<tree::GenParticle> *genParticles=0;
	tree->SetBranchAddress("signal_m1", &signal_m1);
	tree->SetBranchAddress("signal_m2", &signal_m2);
	tree->SetBranchAddress("intermediateGenParticles", &intermediateGenParticles);
	tree->SetBranchAddress("genParticles", &genParticles);
	
	
	TH2F histBosons("",";M1 (GeV);M2 (GeV);number of neutralino decays",30,25,1525,30,25,1525);
	TH2F nEvents("",";M1 (GeV);M2 (GeV);number of events",30,25,1525,30,25,1525);
	
	Long64_t iEvents = tree->GetEntries();
	int processEvents=cfg.processFraction*iEvents;
	int CC = 0;
	int total = 0;
	for (int iEvent=0; iEvent<iEvents; iEvent++){
		if (iEvent>processEvents) break;
		if (iEvent%(iEvents/100)==0) {io::log*"."; io::log.flush(); };
		tree->GetEvent(iEvent);
		bool firstChar = false;
		bool secondChar = false;
		if(signal_m1!=500 || signal_m2!=1000) continue;
		//~ if(signal_m1!=700 || signal_m2!=600) continue;
		signal_nWBosons = 0;
		//~ std::cout<<"-------------------"<<std::endl;
		for (tree::IntermediateGenParticle p :*intermediateGenParticles) {
			if (abs(p.pdgId)==24 && p.daughters.size()==2) signal_nWBosons++;
		}
		for (tree::GenParticle p :*genParticles) {
			if (p.pdgId==1000024) firstChar=true;
			if (p.pdgId==-1000024) secondChar=true;
		}
		if (firstChar && secondChar) {
			CC++;
		}
		else {
			if (signal_nWBosons!=1) {
				std::cout<<"----------------------"<<std::endl;
				for (tree::GenParticle p :*genParticles) {
					std::cout<<p.pdgId<<std::endl;
				}
			}
		}
		total++;
		
		histBosons.Fill(signal_m1,signal_m2,signal_nWBosons);
		nEvents.Fill(signal_m1,signal_m2);
		
	}
	std::cout<<CC<<" "<<signal_nWBosons<<" "<<total<<std::endl;
	
	histBosons.Divide(&nEvents);
	
	//~ TFile out("../output/stuff/GGM_nWBosons.root","update");
	//~ histBosons.Write(TString::Format("nWBosons%.1f",cfg.processFraction*100),TObject::kOverwrite);
	//~ out.Close();
	//~ file.Close();
}

extern "C"
void run()
{     
	//~ nDecays();
	//~ nWBosons();
	dominatProcesses();
}
