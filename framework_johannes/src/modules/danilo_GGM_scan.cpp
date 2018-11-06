//Investigating kinematic variables for the two GGM scans
#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TObject.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGraph2D.h>
#include <TTreeReader.h>
#include <cmath>

Config const &cfg=Config::get();

void M1_M2() {
	TFile file("/net/data_cms1b/user/dmeuser/master/data/v03D/GGM_GravitinoLSP_M1-200to1500_M2-200to1500.root","read");
	
	
	TTree *tree=(TTree*)file.Get(cfg.treeName);
	UShort_t signal_nNeutralinoDecays = 0;
	UShort_t signal_nBinos = 0;
	UShort_t signal_m1 = 0;
	UShort_t signal_m2 = 0;
	std::vector<tree::Photon> *photons=0;
	std::vector<tree::Jet> *jets=0;
	tree::MET *MET=0;
	tree->SetBranchAddress("signal_nNeutralinoDecays", &signal_nNeutralinoDecays);
	tree->SetBranchAddress("signal_nBinos", &signal_nBinos);
	tree->SetBranchAddress("signal_m1", &signal_m1);
	tree->SetBranchAddress("signal_m2", &signal_m2);
	tree->SetBranchAddress("photons", &photons);
	tree->SetBranchAddress("met", &MET);
	tree->SetBranchAddress("jets", &jets);
	
	TH2F phoPt("PhotonPt",";M1 (GeV);M2 (GeV);mean #gamma_{1} PT (GeV)",30,25,1525,30,25,1525);
	TH2F metHist("MET",";M1 (GeV);M2 (GeV);mean MET (GeV)",30,25,1525,30,25,1525);
	TH2F stHist("ST",";M1 (GeV);M2 (GeV);mean ST (GeV)",30,25,1525,30,25,1525);
	TH2F mtHist("MT",";M1 (GeV);M2 (GeV);mean MT (GeV)",30,25,1525,30,25,1525);
	TH2F htgHist("HTG",";M1 (GeV);M2 (GeV);mean HTG (GeV)",30,25,1525,30,25,1525);
	TH2F nEvents("",";M1 (GeV);M2 (GeV);number of events",30,25,1525,30,25,1525);
	
	Long64_t iEvents = tree->GetEntries();
	int processEvents=cfg.processFraction*iEvents;
	for (int iEvent=0; iEvent<iEvents; iEvent++){
		if (iEvent>processEvents) break;
		if (iEvent%(iEvents/100)==0) {io::log*"."; io::log.flush(); };
		tree->GetEvent(iEvent);
		
		//Select loose photons
		std::vector<tree::Photon const *> lPho;
		for (tree::Photon const &ph: *photons){
			if (ph.sigmaIetaIeta<0.001 || ph.sigmaIphiIphi<0.001) continue;
			if (fabs(ph.p.Eta())>1.4442) continue;
			if ((ph.seedCrystalE/ph.p.Pt()) < 0.3) continue;
			if (!ph.hasPixelSeed){
				lPho.push_back(&ph);
			}
		}
		if(lPho.size()==0) continue;
		
		//Calculate various variables
		std::vector<tree::Photon const*> const &pho = lPho;	
		float const MT=phys::M_T(*pho[0],*MET);
		float STg=MET->p.Pt();
		for (auto const &ph: pho){
		 STg+=ph->p.Pt();
		}
		
		float emht = pho[0]->p.Pt();
		for (auto const &jet: *jets) {
		 if (jet.p.Pt() > 30 && fabs(jet.p.Eta()) < 3) {
			if (jet.p.DeltaR(pho[0]->p) > 0.3) {
			   emht += jet.p.Pt();
			}
		 }
		}

		
		phoPt.Fill(signal_m1,signal_m2,pho[0]->p.Pt());
		metHist.Fill(signal_m1,signal_m2,MET->p.Pt());
		stHist.Fill(signal_m1,signal_m2,STg);
		mtHist.Fill(signal_m1,signal_m2,MT);
		htgHist.Fill(signal_m1,signal_m2,emht);
		nEvents.Fill(signal_m1,signal_m2);
		
	}
	io::RootFileSaver aux_saver(TString::Format("../output/GGM_kinematics_%.1f.root",cfg.processFraction*100),"");
	for(TH2F hist:{phoPt,metHist,stHist,mtHist,htgHist}){
		hist.Divide(&nEvents);
		aux_saver.save(hist,TString::Format("M1_M2/%s",hist.GetName()));
	}
		
	
	//~ io::RootFileSaver saver("plots.root",TString::Format("danilo_GGM_scan%.1f/%s",cfg.processFraction*100,"Limits"));
	
	//~ for(TH2F hist:{phoPt,metHist,stHist,mtHist,htgHist}){
	
		//~ hist.Divide(&nEvents);
		//~ TCanvas can;
		//~ can.SetLogz();
		//~ gPad->SetRightMargin(0.2);
		//~ gPad->SetLeftMargin(0.13);
		//~ gPad->SetBottomMargin(0.10);
		
		//~ hist.GetYaxis()->SetTitleOffset(1.3);
		//~ hist.GetXaxis()->SetTitleOffset(0.9);
		//~ hist.GetZaxis()->SetTitleOffset(1.3);
		//~ hist.GetYaxis()->SetTitleSize(0.05);
		//~ hist.GetXaxis()->SetTitleSize(0.05);
		//~ hist.GetZaxis()->SetTitleSize(0.05);
		//~ hist.GetYaxis()->SetLabelSize(0.04);
		//~ hist.GetXaxis()->SetLabelSize(0.04);
		//~ hist.GetZaxis()->SetLabelSize(0.04);
		//~ hist.GetZaxis()->SetLabelOffset(0.01);
		//~ hist.SetAxisRange(250,1500,"X");
		//~ hist.SetAxisRange(250,1500,"Y");  
		//~ hist.SetStats(false);
		//~ hist.Draw("colz");
		
		
		//~ saver.save(can,TString::Format("M1_M2/%s",hist.GetName()),true,true);
	//~ }

}

void M1_M3() {
	TFile file("/net/data_cms1b/user/dmeuser/master/data/v03D/GGM_GravitinoLSP_M1-50to1500_M3-1000to2500.root","read");
	
	
	TTree *tree=(TTree*)file.Get(cfg.treeName);
	UShort_t signal_nNeutralinoDecays = 0;
	UShort_t signal_nBinos = 0;
	UShort_t signal_m1 = 0;
	UShort_t signal_m2 = 0;
	std::vector<tree::Photon> *photons=0;
	std::vector<tree::Jet> *jets=0;
	tree::MET *MET=0;
	tree->SetBranchAddress("signal_nNeutralinoDecays", &signal_nNeutralinoDecays);
	tree->SetBranchAddress("signal_nBinos", &signal_nBinos);
	tree->SetBranchAddress("signal_m1", &signal_m1);
	tree->SetBranchAddress("signal_m2", &signal_m2);
	tree->SetBranchAddress("photons", &photons);
	tree->SetBranchAddress("met", &MET);
	tree->SetBranchAddress("jets", &jets);
	
	TH2F phoPt("PhotonPt",";M1 (GeV);M3 (GeV);mean photon pT (GeV)",30,25,1525,31,975,2525);
	TH2F metHist("MET",";M1 (GeV);M3 (GeV);mean MET (GeV)",30,25,1525,31,975,2525);
	TH2F stHist("ST",";M1 (GeV);M3 (GeV);mean ST (GeV)",30,25,1525,31,975,2525);
	TH2F mtHist("MT",";M1 (GeV);M3 (GeV);mean MT (GeV)",30,25,1525,31,975,2525);
	TH2F htgHist("HTG",";M1 (GeV);M3 (GeV);mean HTG (GeV)",30,25,1525,31,975,2525);
	TH2F nEvents("",";M1 (GeV);M3 (GeV);number of events",30,25,1525,31,975,2525);
	
	Long64_t iEvents = tree->GetEntries();
	int processEvents=cfg.processFraction*iEvents;
	for (int iEvent=0; iEvent<iEvents; iEvent++){
		if (iEvent>processEvents) break;
		if (iEvent%(iEvents/100)==0) {io::log*"."; io::log.flush(); };
		tree->GetEvent(iEvent);
		
		//Select loose photons
		std::vector<tree::Photon const *> lPho;
		for (tree::Photon const &ph: *photons){
			if (ph.sigmaIetaIeta<0.001 || ph.sigmaIphiIphi<0.001) continue;
			if (fabs(ph.p.Eta())>1.4442) continue;
			if ((ph.seedCrystalE/ph.p.Pt()) < 0.3) continue;
			if (!ph.hasPixelSeed){
				lPho.push_back(&ph);
			}
		}
		if(lPho.size()==0) continue;
		
		//Calculate various variables
		std::vector<tree::Photon const*> const &pho = lPho;	
		float const MT=phys::M_T(*pho[0],*MET);
		float STg=MET->p.Pt();
		for (auto const &ph: pho){
		 STg+=ph->p.Pt();
		}
		
		float emht = pho[0]->p.Pt();
		for (auto const &jet: *jets) {
		 if (jet.p.Pt() > 30 && fabs(jet.p.Eta()) < 3) {
			if (jet.p.DeltaR(pho[0]->p) > 0.3) {
			   emht += jet.p.Pt();
			}
		 }
		}

		
		phoPt.Fill(signal_m1,signal_m2,pho[0]->p.Pt());
		metHist.Fill(signal_m1,signal_m2,MET->p.Pt());
		stHist.Fill(signal_m1,signal_m2,STg);
		mtHist.Fill(signal_m1,signal_m2,MT);
		htgHist.Fill(signal_m1,signal_m2,emht);
		nEvents.Fill(signal_m1,signal_m2);
		
	}
	
	io::RootFileSaver aux_saver(TString::Format("../output/GGM_kinematics_%.1f.root",cfg.processFraction*100),"");
	for(TH2F hist:{phoPt,metHist,stHist,mtHist,htgHist}){
		hist.Divide(&nEvents);
		aux_saver.save(hist,TString::Format("M1_M3/%s",hist.GetName()));
	}
	
	//~ io::RootFileSaver saver("plots.root",TString::Format("danilo_GGM_scan%.1f/%s",cfg.processFraction*100,"Limits"));
	
	//~ for(TH2F hist:{phoPt,metHist,stHist,mtHist,htgHist}){
	
		//~ hist.Divide(&nEvents);
		//~ TCanvas can;
		//~ can.SetLogz();
		//~ gPad->SetRightMargin(0.2);
		//~ gPad->SetLeftMargin(0.13);
		//~ gPad->SetBottomMargin(0.10);
		
		//~ hist.GetYaxis()->SetTitleOffset(1.3);
		//~ hist.GetXaxis()->SetTitleOffset(0.9);
		//~ hist.GetZaxis()->SetTitleOffset(1.3);
		//~ hist.GetYaxis()->SetTitleSize(0.05);
		//~ hist.GetXaxis()->SetTitleSize(0.05);
		//~ hist.GetZaxis()->SetTitleSize(0.05);
		//~ hist.GetYaxis()->SetLabelSize(0.04);
		//~ hist.GetXaxis()->SetLabelSize(0.04);
		//~ hist.GetZaxis()->SetLabelSize(0.04);
		//~ hist.GetZaxis()->SetLabelOffset(0.01);
		//~ hist.SetAxisRange(50,1500,"X");
		//~ hist.SetAxisRange(1000,2500,"Y");  
		//~ hist.SetStats(false);
		//~ hist.Draw("colz");
		
		
		//~ saver.save(can,TString::Format("M1_M3/%s",hist.GetName()),true,true);
	//~ }

}

void T5Wg() {
	TFile file("/net/data_cms1b/user/dmeuser/master/data/v03D/SMS-T5Wg.root","read");
	
	
	TTree *tree=(TTree*)file.Get(cfg.treeName);
	UShort_t signal_nNeutralinoDecays = 0;
	UShort_t signal_nBinos = 0;
	UShort_t signal_m1 = 0;
	UShort_t signal_m2 = 0;
	std::vector<tree::Photon> *photons=0;
	std::vector<tree::Jet> *jets=0;
	tree::MET *MET=0;
	tree->SetBranchAddress("signal_nNeutralinoDecays", &signal_nNeutralinoDecays);
	tree->SetBranchAddress("signal_nBinos", &signal_nBinos);
	tree->SetBranchAddress("signal_m1", &signal_m1);
	tree->SetBranchAddress("signal_m2", &signal_m2);
	tree->SetBranchAddress("photons", &photons);
	tree->SetBranchAddress("met", &MET);
	tree->SetBranchAddress("jets", &jets);
	
	TH2F phoPt("PhotonPt",";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV);mean #gamma_{1} PT (GeV)",18,0,2500,21,0,2150);
	TH2F metHist("MET",";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV);mean MET (GeV)",18,0,2500,21,0,2150);
	TH2F stHist("ST",";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV);mean ST (GeV)",18,0,2500,21,0,2150);
	TH2F mtHist("MT",";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV);mean MT (GeV)",18,0,2500,21,0,2150);
	TH2F htgHist("HTG",";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV);mean HTG (GeV)",18,0,2500,21,0,2150);
	TH2F nEvents("",";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV);number of events",18,0,2500,21,0,2150);
	
	Long64_t iEvents = tree->GetEntries();
	int processEvents=cfg.processFraction*iEvents;
	for (int iEvent=0; iEvent<iEvents; iEvent++){
		if (iEvent>processEvents) break;
		if (iEvent%(iEvents/100)==0) {io::log*"."; io::log.flush(); };
		tree->GetEvent(iEvent);
		
		//Select loose photons
		std::vector<tree::Photon const *> lPho;
		for (tree::Photon const &ph: *photons){
			if (ph.sigmaIetaIeta<0.001 || ph.sigmaIphiIphi<0.001) continue;
			if (fabs(ph.p.Eta())>1.4442) continue;
			if ((ph.seedCrystalE/ph.p.Pt()) < 0.3) continue;
			if (!ph.hasPixelSeed){
				lPho.push_back(&ph);
			}
		}
		if(lPho.size()==0) continue;
		
		//Calculate various variables
		std::vector<tree::Photon const*> const &pho = lPho;	
		float const MT=phys::M_T(*pho[0],*MET);
		float STg=MET->p.Pt();
		for (auto const &ph: pho){
		 STg+=ph->p.Pt();
		}
		
		float emht = pho[0]->p.Pt();
		for (auto const &jet: *jets) {
		 if (jet.p.Pt() > 30 && fabs(jet.p.Eta()) < 3) {
			if (jet.p.DeltaR(pho[0]->p) > 0.3) {
			   emht += jet.p.Pt();
			}
		 }
		}

		
		phoPt.Fill(signal_m1,signal_m2,pho[0]->p.Pt());
		metHist.Fill(signal_m1,signal_m2,MET->p.Pt());
		stHist.Fill(signal_m1,signal_m2,STg);
		mtHist.Fill(signal_m1,signal_m2,MT);
		htgHist.Fill(signal_m1,signal_m2,emht);
		nEvents.Fill(signal_m1,signal_m2);
		
	}
	io::RootFileSaver aux_saver(TString::Format("../output/GGM_kinematics_%.1f.root",cfg.processFraction*100),"");
	for(TH2F hist:{phoPt,metHist,stHist,mtHist,htgHist}){
		hist.Divide(&nEvents);
		aux_saver.save(hist,TString::Format("T5Wg/%s",hist.GetName()));
	}
}


void TChiNg() {
	TFile file("/net/data_cms1b/user/dmeuser/master/data/v03D/SMS-TChiNG_BF50N50G.root","read");
	
	
	TTree *tree=(TTree*)file.Get(cfg.treeName);
	UShort_t signal_nNeutralinoDecays = 0;
	UShort_t signal_nBinos = 0;
	UShort_t signal_m1 = 0;
	UShort_t signal_m2 = 0;
	std::vector<tree::Photon> *photons=0;
	std::vector<tree::Jet> *jets=0;
	tree::MET *MET=0;
	tree->SetBranchAddress("signal_nNeutralinoDecays", &signal_nNeutralinoDecays);
	tree->SetBranchAddress("signal_nBinos", &signal_nBinos);
	tree->SetBranchAddress("signal_m1", &signal_m1);
	tree->SetBranchAddress("signal_m2", &signal_m2);
	tree->SetBranchAddress("photons", &photons);
	tree->SetBranchAddress("met", &MET);
	tree->SetBranchAddress("jets", &jets);
	
	TH1F phoPt("PhotonPt",";NLSPmass;x",41,300-12.5,1300+12.5);	
	TH1F metHist("MET",";NLSPmass;x",41,300-12.5,1300+12.5);	
	TH1F stHist("ST",";NLSPmass;x",41,300-12.5,1300+12.5);	
	TH1F mtHist("MT",";NLSPmass;x",41,300-12.5,1300+12.5);	
	TH1F htgHist("HTG",";NLSPmass;x",41,300-12.5,1300+12.5);	
	TH1F nEvents("",";NLSPmass;x",41,300-12.5,1300+12.5);	
	
	Long64_t iEvents = tree->GetEntries();
	int processEvents=cfg.processFraction*iEvents;
	for (int iEvent=0; iEvent<iEvents; iEvent++){
		if (iEvent>processEvents) break;
		if (iEvent%(iEvents/100)==0) {io::log*"."; io::log.flush(); };
		tree->GetEvent(iEvent);
		
		//Select loose photons
		std::vector<tree::Photon const *> lPho;
		for (tree::Photon const &ph: *photons){
			if (ph.sigmaIetaIeta<0.001 || ph.sigmaIphiIphi<0.001) continue;
			if (fabs(ph.p.Eta())>1.4442) continue;
			if ((ph.seedCrystalE/ph.p.Pt()) < 0.3) continue;
			if (!ph.hasPixelSeed){
				lPho.push_back(&ph);
			}
		}
		if(lPho.size()==0) continue;
		
		//Calculate various variables
		std::vector<tree::Photon const*> const &pho = lPho;	
		float const MT=phys::M_T(*pho[0],*MET);
		float STg=MET->p.Pt();
		for (auto const &ph: pho){
		 STg+=ph->p.Pt();
		}
		
		float emht = pho[0]->p.Pt();
		for (auto const &jet: *jets) {
		 if (jet.p.Pt() > 30 && fabs(jet.p.Eta()) < 3) {
			if (jet.p.DeltaR(pho[0]->p) > 0.3) {
			   emht += jet.p.Pt();
			}
		 }
		}

		//~ std::cout<<signal_m1<<"    "<<MET->p.Pt()<<std::endl;
		phoPt.Fill(signal_m1,pho[0]->p.Pt());
		metHist.Fill(signal_m1,MET->p.Pt());
		stHist.Fill(signal_m1,STg);
		mtHist.Fill(signal_m1,MT);
		htgHist.Fill(signal_m1,emht);
		nEvents.Fill(signal_m1);
		
	}
	io::RootFileSaver aux_saver(TString::Format("../output/GGM_kinematics_%.1f.root",cfg.processFraction*100),"");
	for(TH1F hist:{phoPt,metHist,stHist,mtHist,htgHist}){
		hist.Divide(&nEvents);
		aux_saver.save(hist,TString::Format("TChiNg/%s",hist.GetName()));
	}
}

extern "C"
void run()
{     
	//~ M1_M2();
	//~ M1_M3();
	//~ T5Wg();
	TChiNg();
}
