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
	
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_postfitBKG%.1f/%s",cfg.processFraction*100,"Limits"));
	//~ io::RootFileSaver saver("test.root",TString::Format("danilo_plot_postfitBKG%.1f/%s",cfg.processFraction*100,"Limits"));
	
	TFile file("../input/hists_prefit_bkg.root","read");
	
	for (TString bkg :{"GJ","TTcomb","Vg","diboson","efake"}){
		TCanvas can;
		can.SetGrid();
		gStyle->SetHistMinimumZero();
		TH1F *hist_noCorr = (TH1F*) file.Get("ST_"+bkg+"_noCorr");
		TH1F *hist_VgGjet = (TH1F*) file.Get("ST_"+bkg+"_VgGjet");
		TH1F *hist_allCorr = (TH1F*) file.Get("ST_"+bkg+"_allCorr");
		
		hist_noCorr->SetTitle(bkg);
		hist_noCorr->SetMinimum(-1);
		hist_noCorr->SetMaximum(1);
		hist_noCorr->SetFillColor(4);
		hist_noCorr->SetBarWidth(0.3);
		hist_noCorr->SetBarOffset(0.05);
		hist_noCorr->SetStats(0);
		for (int i=1; i<=4; i++) hist_noCorr->GetXaxis()->SetBinLabel(i,(std::to_string(i)).c_str());
		hist_noCorr->GetXaxis()->SetLabelSize(0.1);
		hist_noCorr->Draw("b");
		
		hist_VgGjet->SetFillColor(38);
		hist_VgGjet->SetBarWidth(0.3);
		hist_VgGjet->SetBarOffset(0.35);
		hist_VgGjet->SetStats(0);
		hist_VgGjet->Draw("b same");
		
		hist_allCorr->SetFillColor(9);
		hist_allCorr->SetBarWidth(0.3);
		hist_allCorr->SetBarOffset(0.65);
		hist_allCorr->SetStats(0);
		hist_allCorr->Draw("b same");
		
		gfx::LegendEntries legE;
		legE.append(*hist_noCorr,"no corr","f");
		legE.append(*hist_VgGjet,"VgGjet corr","f");
		legE.append(*hist_allCorr,"VgGjet/efake/rare corr","f");
		double y1=0.7;
		double y2=0.9;
		if (bkg=="Vg") {
			y1=0.2;
			y2=0.4;
			}
		TLegend leg=legE.buildLegend(.2,y1,0.92,y2,1);
		leg.SetTextSize(0.03);
		leg.Draw();
		saver.save(can,"ST_"+bkg,true,false);
		
	}
	
	for (TString bkg :{"VGamma","elefakepho","jetfakepho","qcdfakelep","rare"}){
		TCanvas can;
		can.SetGrid();
		gStyle->SetHistMinimumZero();
		TH1F *hist_noCorr = (TH1F*) file.Get("Lepton_"+bkg+"_noCorr");
		TH1F *hist_VgGjet = (TH1F*) file.Get("Lepton_"+bkg+"_VgGjet");
		TH1F *hist_allCorr = (TH1F*) file.Get("Lepton_"+bkg+"_allCorr");
		
		hist_noCorr->SetTitle(bkg);
		hist_noCorr->SetMinimum(-1.65);
		hist_noCorr->SetMaximum(1.65);
		hist_noCorr->SetFillColor(4);
		hist_noCorr->SetBarWidth(0.3);
		hist_noCorr->SetBarOffset(0.05);
		hist_noCorr->SetStats(0);
		for (int i=1; i<=36; i++) hist_noCorr->GetXaxis()->SetBinLabel(i,(std::to_string(i)).c_str());
		hist_noCorr->Draw("b");
		
		hist_VgGjet->SetFillColor(38);
		hist_VgGjet->SetBarWidth(0.3);
		hist_VgGjet->SetBarOffset(0.35);
		hist_VgGjet->SetStats(0);
		hist_VgGjet->Draw("b same");
		
		hist_allCorr->SetFillColor(9);
		hist_allCorr->SetBarWidth(0.3);
		hist_allCorr->SetBarOffset(0.65);
		hist_allCorr->SetStats(0);
		hist_allCorr->Draw("b same");
		
		gfx::LegendEntries legE;
		legE.append(*hist_noCorr,"no corr","f");
		legE.append(*hist_VgGjet,"VgGjet corr","f");
		legE.append(*hist_allCorr,"VgGjet/efake/rare corr","f");
		double y1=0.7;
		double y2=0.9;
		if (bkg=="VGamma") {
			y1=0.2;
			y2=0.4;
			}
		TLegend leg=legE.buildLegend(.2,y1,0.92,y2,1);
		leg.SetTextSize(0.03);
		leg.Draw();
		
		saver.save(can,"Lepton_"+bkg,true,false);
		
	}
	return 0;
}

extern "C"

void run(){
	plot();
}
