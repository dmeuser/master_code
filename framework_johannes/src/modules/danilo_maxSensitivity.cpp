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
	
	for (TString scan: {"GGM_M1_M2","GGM_M1_M3","T5Wg","TChiNg_BR","CharginoBR_C1C1"}) {
		TFile file_1("../input/limits/limits_"+scan+"_inclusivFinal.root","read");
		TFile file_2("../input/limits/limits_"+scan+"_htgFinal.root","read");
		TFile file_3("../input/limits/limits_"+scan+"_lepton_final.root","read");
		TFile file_4("../input/limits/limits_"+scan+"_diphoton_final.root","read");
		
		TH2F *st = (TH2F*) file_1.Get("h_exp");
		TH2F *htg = (TH2F*) file_2.Get("h_exp");
		TH2F *lepton = (TH2F*) file_3.Get("h_exp");
		TH2F *diphoton = (TH2F*) file_4.Get("h_exp");
		TH2I maxSen;
		
		if (scan=="GGM_M1_M2") maxSen = TH2I("",";#it{M}_{1} (GeV);#it{M}_{2} (GeV)",27,175,1525,27,175,1525);
		else if (scan=="GGM_M1_M3") maxSen = TH2I("",";#it{M}_{1} (GeV);#it{M}_{3} (GeV)",30,25,1525,31,975,2525);
		else if (scan=="TChiNg_BR") maxSen = TH2I("",";BF(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #gamma + #tilde{G}) (%);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)",51,-1,101,41,287.5,1312.5);
		else if (scan=="CharginoBR_C1C1") maxSen = TH2I("",";BF(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #tilde{#chi}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}(#gamma#tilde{G}) + soft) (%);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)",51,-1,101,41,287.5,1312.5);
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
		
        if (scan=="GGM_M1_M2"){
            maxSen.SetAxisRange(250,1500,"X");
            maxSen.SetAxisRange(250,1500,"Y");
        }
        else if (scan=="T5Wg"){
            maxSen.SetAxisRange(1400,2100,"X");
            maxSen.SetAxisRange(0,2500,"Y");
        }
        else if (scan=="TChiNg_BR" or scan=="CharginoBR_C1C1"){
            maxSen.SetAxisRange(1,100,"X");
            maxSen.SetAxisRange(300,1280,"Y");
        }
		maxSen.SetMaximum(4);
		maxSen.SetMinimum(1);
		maxSen.SetStats(0);
		maxSen.Draw("col");
		
		//~ TLatex text1(0.4*maxSen.GetXaxis()->GetXmax(),0.2*maxSen.GetYaxis()->GetXmax(),"Diphoton");
		//~ TLatex text2(0.4*maxSen.GetXaxis()->GetXmax(),0.4*maxSen.GetYaxis()->GetXmax(),"Lepton");
		//~ TLatex text3(0.4*maxSen.GetXaxis()->GetXmax(),0.6*maxSen.GetYaxis()->GetXmax(),"S_{#scale[.8]{T}}^{#scale[.8]{#gamma}}");
		//~ TLatex text4(0.4*maxSen.GetXaxis()->GetXmax(),0.8*maxSen.GetYaxis()->GetXmax(),"H_{#scale[.8]{T}}^{#scale[.8]{#gamma}}");
		//~ text1.SetTextSize(0.03);
		//~ text2.SetTextSize(0.03);
		//~ text3.SetTextSize(0.03);
		//~ text4.SetTextSize(0.03);
		//~ text3.SetTextColor(kWhite);
		//~ text1.Draw("same");
		//~ text2.Draw("same");
		//~ text3.Draw("same");
		//~ text4.Draw("same");
		
		TH2F *diphotonLeg = new TH2F();
		diphotonLeg->SetFillColor(TColor::GetColor(249,249,15));
		diphotonLeg->SetFillStyle(1001);
		diphotonLeg->SetLineWidth(0);
		
		TH2F *leptonLeg = new TH2F();
		leptonLeg->SetFillColor(TColor::GetColor(160,190,109));
		leptonLeg->SetFillStyle(1001);
		leptonLeg->SetLineWidth(0);
		
		TH2F *stLeg = new TH2F();
		stLeg->SetFillColor(TColor::GetColor(52,42,135));
		stLeg->SetFillStyle(1001);
		stLeg->SetLineWidth(0);
		
		TH2F *htgLeg = new TH2F();
		htgLeg->SetFillColor(TColor::GetColor(10,153,206));
		htgLeg->SetFillStyle(1001);
		htgLeg->SetLineWidth(0);
		
		gfx::LegendEntries legE;
		legE.append(*stLeg,"Photon+S_{#scale[.8]{T}}^{#scale[.8]{#gamma}}","f");
		legE.append(*htgLeg,"Photon+H_{#scale[.8]{T}}^{#scale[.8]{#gamma}}","f");
		legE.append(*leptonLeg,"Photon+Lepton","f");
		legE.append(*diphotonLeg,"Diphoton","f");
		TLegend leg=legE.buildLegend(.2,.75,0.65,0.9,2);
		leg.SetFillColor(10);
		leg.SetFillStyle(1001);
		leg.SetTextSize(0.03);
		leg.Draw("f");
			
		saver.save(can,scan,true,false);
	}
	
	return 0;
}

int plot_physmass(){
	
	io::RootFileSaver saver("plots.root",TString::Format("danilo_maxSensitivity%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	//~ for (TString scan: {"GGM_M1_M2_NN","GGM_M1_M3_NN"}) {
	for (TString scan: {"GGM_M1_M2_final","GGM_M1_M3_final"}) {
		TFile file_1("../input/limits/physmass_"+scan+".root","read");
		
		TH2F *st = (TH2F*) file_1.Get("inclusivFinal/hist_Inter");
		TH2F *htg = (TH2F*) file_1.Get("htgFinal/hist_Inter");
		TH2F *lepton = (TH2F*) file_1.Get("lepton_final/hist_Inter");
		TH2F *diphoton = (TH2F*) file_1.Get("diphoton_final/hist_Inter");
		TH2I maxSen;
		
		if (scan=="GGM_M1_M2_final") maxSen = TH2I("",";m_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV);m_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV)",50,50,725,100,135,1240);
		else if (scan=="GGM_M1_M3_final") maxSen = TH2I("",";m_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV)",30,30,750,100,2350,5300);
				
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
					//~ std::cout<<temps.front()<<std::endl;
					maxSen.SetBinContent(i,j,std::distance(temps.begin(),std::min_element(temps.begin(),temps.end()))+1);
				}
			}
		}
		
		maxSen.SetMaximum(4);
		maxSen.SetMinimum(1);
		maxSen.SetStats(0);
		maxSen.Draw("col");
		
		TH2F *diphotonLeg = new TH2F();
		diphotonLeg->SetFillColor(TColor::GetColor(249,249,15));
		diphotonLeg->SetFillStyle(1001);
		diphotonLeg->SetLineWidth(0);
		
		TH2F *leptonLeg = new TH2F();
		leptonLeg->SetFillColor(TColor::GetColor(160,190,109));
		leptonLeg->SetFillStyle(1001);
		leptonLeg->SetLineWidth(0);
		
		TH2F *stLeg = new TH2F();
		stLeg->SetFillColor(TColor::GetColor(52,42,135));
		stLeg->SetFillStyle(1001);
		stLeg->SetLineWidth(0);
		
		TH2F *htgLeg = new TH2F();
		htgLeg->SetFillColor(TColor::GetColor(10,153,206));
		htgLeg->SetFillStyle(1001);
		htgLeg->SetLineWidth(0);
		
		gfx::LegendEntries legE;
		legE.append(*stLeg,"Photon+S_{#scale[.8]{T}}^{#scale[.8]{#gamma}}","f");
		legE.append(*htgLeg,"Photon+H_{#scale[.8]{T}}^{#scale[.8]{#gamma}}","f");
		legE.append(*leptonLeg,"Photon+Lepton","f");
		legE.append(*diphotonLeg,"Diphoton","f");
		TLegend leg=legE.buildLegend(.2,.75,0.65,0.9,2);
		leg.SetFillColor(10);
		leg.SetFillStyle(1001);
		leg.SetTextSize(0.03);
		leg.Draw("f");
			
		saver.save(can,scan+"_physmass",true,false);
	}
	
	return 0;
}

extern "C"

void run(){
	plot_M1M2();
	//~ plot_physmass();
}
