#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/util.hpp"

/*
 * Plotting multiple limts in one canvas
 */
 
#include <TChain.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TLine.h>
#include <TTreeReader.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TColor.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>

using namespace std;

static Config const &cfg=Config::get();

int plot(){
	
	//~ io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_postfitBKG%.1f/%s",cfg.processFraction*100,"Limits"));
	io::RootFileSaver saver("test.root",TString::Format("danilo_plot_postfitBKG%.1f/%s",cfg.processFraction*100,"Limits"));
	
	TFile file("../input/hists_prefit_bkg.root","read");
	
	for (TString bkg :{"GJ","TTcomb","Vg","diboson","efake"}){
		TCanvas can;
		TH1F *hist_noCorr = (TH1F*) file.Get("ST_"+bkg+"_noCorr");
		TH1F *hist_VgGjet = (TH1F*) file.Get("ST_"+bkg+"_VgGjet");
		TH1F *hist_allCorr = (TH1F*) file.Get("ST_"+bkg+"_allCorr");
		
		hist_noCorr->Draw();
		hist_VgGjet->Draw("same");
		hist_allCorr->Draw("same");
		saver.save(can,"ST_"+bkg,true,false);
		
	}
	
	for (TString bkg :{"VGamma","elefakepho","jetfakepho","qcdfakelep","rare"}){
		TCanvas can;
		TH1F *hist_noCorr = (TH1F*) file.Get("Lepton_"+bkg+"_noCorr");
		TH1F *hist_VgGjet = (TH1F*) file.Get("Lepton_"+bkg+"_VgGjet");
		TH1F *hist_allCorr = (TH1F*) file.Get("Lepton_"+bkg+"_allCorr");
		
		hist_noCorr->Draw();
		hist_VgGjet->Draw("same");
		hist_allCorr->Draw("same");
		saver.save(can,"Lepton_"+bkg,true,false);
		
	}
	return 0;
}

extern "C"

void run(){
	plot();
}
