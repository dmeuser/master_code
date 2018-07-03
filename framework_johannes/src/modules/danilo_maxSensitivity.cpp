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
#include <TImage.h>
#include <TLatex.h>


using namespace std;

static Config const &cfg=Config::get();

bool valid(float a) {
	return a!=9999;
}

int plot_M1M2(){
	
	io::RootFileSaver saver("plots.root",TString::Format("danilo_maxSensitivity%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	for (TString scan: {"GGM_M1_M2","GGM_M1_M3","T5Wg"}) {
		TFile file_1("../input/limits/limits_"+scan+"_inclusiv.root","read");
		TFile file_2("../input/limits/limits_"+scan+"_htg.root","read");
		TFile file_3("../input/limits/limits_"+scan+"_lepton.root","read");
		TFile file_4("../input/limits/limits_"+scan+"_diphoton.root","read");
		
		TH2F *st = (TH2F*) file_1.Get("h_exp");
		TH2F *htg = (TH2F*) file_2.Get("h_exp");
		TH2F *lepton = (TH2F*) file_3.Get("h_exp");
		TH2F *diphoton = (TH2F*) file_4.Get("h_exp");
		TH2I maxSen;
		
		if (scan=="GGM_M1_M2") maxSen = TH2I("",";M_{1} (GeV);M_{2} (GeV)",27,175,1525,27,175,1525);
		else if (scan=="GGM_M1_M3") maxSen = TH2I("",";M_{1} (GeV);M_{3} (GeV)",30,25,1525,31,975,2525);
		else if (scan=="T5Wg") maxSen = TH2I("",";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)",18,0,2500,21,0,2150);
		
		float st_temp;
		float htg_temp;
		float lepton_temp;
		float diphoton_temp;
		
		for (int i = 1; i<=st->GetNbinsX(); i++){
			for (int j = 1; j<=st->GetNbinsY(); j++){
				st_temp = st->GetBinContent(i,j);
				htg_temp = htg->GetBinContent(i,j);
				lepton_temp = lepton->GetBinContent(i,j);
				diphoton_temp = diphoton->GetBinContent(i,j);
				
				std::list<float> temps = {st_temp,htg_temp,lepton_temp,diphoton_temp};
				std::replace(temps.begin(),temps.end(),0,9999);
				if (*(std::find_if(temps.begin(),temps.end(),valid)) > 0.001) {
					maxSen.SetBinContent(i,j,std::distance(temps.begin(),std::min_element(temps.begin(),temps.end()))+1);
				}
			}
		}
		
		maxSen.SetMaximum(4);
		maxSen.SetMinimum(1);
		maxSen.SetStats(0);
		maxSen.Draw("col");
		
		TLatex text1(400,1000,"Diphoton");
		TLatex text2(750,700,"Lepton");
		TLatex text3(1000,600,"ST");
		TLatex text4(1250,750,"HTG");
		text1.SetTextSize(0.03);
		text2.SetTextSize(0.03);
		text3.SetTextSize(0.03);
		text4.SetTextSize(0.03);
		text3.SetTextColor(kWhite);
		text1.Draw("same");
		text2.Draw("same");
		text3.Draw("same");
		text4.Draw("same");
			
		saver.save(can,scan,true,false);
	}
	
	return 0;
}

extern "C"

void run(){
	plot_M1M2();
}
