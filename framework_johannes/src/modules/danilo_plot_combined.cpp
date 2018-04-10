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
	
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/T5Wg_limits_mock/limits_T5Wg_johannes.root","read");
	TFile file_2("../input/T5Wg_limits_mock/limits_T5Wg_knut_2.root","read");
	TFile file_3("../input/T5Wg_limits_mock/limits_T5Wg.root","read");
	
	TGraph *johannes_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *knut_exp = (TGraph*) file_2.Get("exp");
	TGraph *combined_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *combined_obs = (TGraph*) file_3.Get("gr_obsC_sm");
	
	johannes_exp->SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	johannes_exp->SetLineColor(kGreen);
	johannes_exp->SetLineWidth(3);
	johannes_exp->SetLineStyle(2);
	johannes_exp->Draw();
	
	knut_exp->SetLineColor(kBlue);
	knut_exp->SetLineWidth(3);
	knut_exp->SetLineStyle(2);
	knut_exp->Draw("same");
	
	combined_exp->SetLineColor(kRed);
	combined_exp->SetLineWidth(3);
	combined_exp->SetLineStyle(2);
	combined_exp->Draw("same");
	
	combined_obs->SetLineColor(kBlack);
	combined_obs->SetLineWidth(3);
	combined_obs->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*johannes_exp,"Photon+ST (exp)","l");
	legE.append(*knut_exp,"Photon+HTG (exp)","l");
	legE.append(*combined_exp,"Combined (exp)","l");
	legE.append(*combined_obs,"Combined (obs)","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_compare_mock",true,false);
	
	return 0;
}

int plot_exclusive(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/T5Wg_limits_mock/limits_T5Wg_johannes.root","read");
	TFile file_2("../input/limits/limits_T5Wg_exclusiv.root","read");
	TFile file_3("../input/limits/limits_T5Wg_leptonVeto.root","read");
	TFile file_4("../input/limits/limits_T5Wg_htgVeto.root","read");
	TFile file_5("../input/limits/limits_T5Wg_diphotonVeto.root","read");
	
	TGraph *johannes_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *exclusiv_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *leptonVeto_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *htgVeto_exp = (TGraph*) file_4.Get("gr_expC_sm");
	TGraph *diphotonVeto_exp = (TGraph*) file_5.Get("gr_expC_sm");
	
	TGraph *leptonAna_exp = new TGraph();
	TGraph *diag = new TGraph();
	
	leptonAna_exp->SetPoint(0,1500,100);
	leptonAna_exp->SetPoint(1,1600,300);
	leptonAna_exp->SetPoint(2,1640,400);
	leptonAna_exp->SetPoint(3,1700,525);
	leptonAna_exp->SetPoint(4,1800,1100);
	leptonAna_exp->SetPoint(5,1840,1600);
	leptonAna_exp->SetPoint(6,1830,1775);
	
	diag->SetPoint(1,500,500);
	diag->SetPoint(2,2200,2200);
	
	johannes_exp->SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	johannes_exp->SetLineColor(kGreen);
	johannes_exp->SetLineWidth(3);
	johannes_exp->SetLineStyle(2);
	johannes_exp->Draw();
	
	exclusiv_exp->SetLineColor(kRed);
	exclusiv_exp->SetLineWidth(3);
	exclusiv_exp->SetLineStyle(2);
	exclusiv_exp->Draw("same");
	
	leptonVeto_exp->SetLineColor(kBlue);
	leptonVeto_exp->SetLineWidth(3);
	leptonVeto_exp->SetLineStyle(2);
	leptonVeto_exp->Draw("same");
	
	htgVeto_exp->SetLineColor(kGray);
	htgVeto_exp->SetLineWidth(3);
	htgVeto_exp->SetLineStyle(2);
	htgVeto_exp->Draw("same");
	
	diphotonVeto_exp->SetLineColor(kMagenta);
	diphotonVeto_exp->SetLineWidth(3);
	diphotonVeto_exp->SetLineStyle(2);
	diphotonVeto_exp->Draw("same");
	
	leptonAna_exp->SetLineColor(kBlack);
	leptonAna_exp->SetLineWidth(3);
	leptonAna_exp->SetLineStyle(2);
	leptonAna_exp->Draw("same");
	
	diag->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*johannes_exp,"Inclusiv (exp)","l");
	legE.append(*exclusiv_exp,"Exclusiv (exp)","l");
	legE.append(*leptonVeto_exp,"LeptonVeto (exp)","l");
	legE.append(*htgVeto_exp,"HTgVeto (exp)","l");
	legE.append(*diphotonVeto_exp,"DiphotonVeto (exp)","l");
	legE.append(*leptonAna_exp,"LeptonAnalysis approx. (exp)","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_compare_exclusive",true,false);
	
	return 0;
}

int plot_TChiNg_BR(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_TChiNg_BR_inclusiv.root","read");
	TFile file_2("../input/limits/limits_TChiNg_BR_exclusiv.root","read");
	TFile file_3("../input/limits/limits_TChiNg_BR_leptonVeto.root","read");
	TFile file_4("../input/limits/limits_TChiNg_BR_htgVeto.root","read");
	TFile file_5("../input/limits/limits_TChiNg_BR_diphotonVeto.root","read");
	
	TGraph *inclusiv_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *exclusiv_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *leptonVeto_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *htgVeto_exp = (TGraph*) file_4.Get("gr_expC_sm");
	TGraph *diphotonVeto_exp = (TGraph*) file_5.Get("gr_expC_sm");
	
	inclusiv_exp->SetTitle(";BR(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #gamma) (%);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	inclusiv_exp->SetLineColor(kGreen);
	inclusiv_exp->SetLineWidth(3);
	inclusiv_exp->SetLineStyle(2);
	inclusiv_exp->Draw();
	
	exclusiv_exp->SetLineColor(kRed);
	exclusiv_exp->SetLineWidth(3);
	exclusiv_exp->SetLineStyle(2);
	exclusiv_exp->Draw("same");
	
	leptonVeto_exp->SetLineColor(kBlue);
	leptonVeto_exp->SetLineWidth(3);
	leptonVeto_exp->SetLineStyle(2);
	leptonVeto_exp->Draw("same");
	
	htgVeto_exp->SetLineColor(kGray);
	htgVeto_exp->SetLineWidth(3);
	htgVeto_exp->SetLineStyle(2);
	htgVeto_exp->Draw("same");
	
	diphotonVeto_exp->SetLineColor(kMagenta);
	diphotonVeto_exp->SetLineWidth(3);
	diphotonVeto_exp->SetLineStyle(2);
	diphotonVeto_exp->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*inclusiv_exp,"Inclusiv (exp)","l");
	legE.append(*exclusiv_exp,"Exclusiv (exp)","l");
	legE.append(*leptonVeto_exp,"LeptonVeto (exp)","l");
	legE.append(*htgVeto_exp,"HTgVeto (exp)","l");
	legE.append(*diphotonVeto_exp,"DiphotonVeto (exp)","l");
	TLegend leg=legE.buildLegend(.7,.2,0.9,.5,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_compare_TChiNg_BR",true,false);
	
	return 0;
}

int plot_htgVeto_FullKnut(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_T5Wg_inclusiv.root","read");
	TFile file_2("../input/T5Wg_limits_mock/limits_T5Wg_knut_2.root","read");
	TFile file_3("../input/limits/limits_T5Wg_htgVeto_KnutFull.root","read");
	TFile file_4("../input/limits/limits_T5Wg_htgHighVeto_KnutHighHtg.root","read");
	
	TGraph *st_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_2.Get("exp");
	TGraph *combined_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *combined_split_exp = (TGraph*) file_4.Get("gr_expC_sm");
	TGraph *diag = new TGraph();
	
	diag->SetPoint(1,500,500);
	diag->SetPoint(2,2200,2200);
	
	st_exp->SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	st_exp->SetLineColor(kGreen);
	st_exp->SetLineWidth(3);
	st_exp->SetLineStyle(2);
	st_exp->Draw();
	
	htg_exp->SetLineColor(kBlue);
	htg_exp->SetLineWidth(3);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	
	combined_exp->SetLineColor(kRed);
	combined_exp->SetLineWidth(3);
	combined_exp->SetLineStyle(2);
	combined_exp->Draw("same");
	
	combined_split_exp->SetLineColor(kMagenta);
	combined_split_exp->SetLineWidth(3);
	combined_split_exp->SetLineStyle(2);
	combined_split_exp->Draw("same");
	
	diag->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*st_exp,"Photon+ST (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*combined_exp,"Combined ST:htgVeto HTG:Full (exp)","l");
	legE.append(*combined_split_exp,"Combined ST:htgHighVeto HTG:High HTG (exp)","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_T5Wg_combine_ST_HTG",true,false);
	
	return 0;
}

int plot_leptonVeto_FullLepton(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_T5Wg_inclusiv.root","read");
	TFile file_2("../input/limits/limits_T5Wg_lepton.root","read");
	TFile file_3("../input/limits/limits_T5Wg_leptonVeto_leptonFull.root","read");
	
	TGraph *st_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *lepton_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *combined_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *diag = new TGraph();
	
	diag->SetPoint(1,500,500);
	diag->SetPoint(2,2200,2200);
	
	st_exp->SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	st_exp->SetLineColor(kGreen);
	st_exp->SetLineWidth(3);
	st_exp->SetLineStyle(2);
	st_exp->Draw();
	
	lepton_exp->SetLineColor(kBlue);
	lepton_exp->SetLineWidth(3);
	lepton_exp->SetLineStyle(2);
	lepton_exp->Draw("same");
	
	combined_exp->SetLineColor(kRed);
	combined_exp->SetLineWidth(3);
	combined_exp->SetLineStyle(2);
	combined_exp->Draw("same");
	
	diag->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*st_exp,"Photon+ST (exp)","l");
	legE.append(*lepton_exp,"Photon+Lepton (exp)","l");
	legE.append(*combined_exp,"Combined ST:leptonVeto Lepton:Full (exp)","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_T5Wg_combine_ST_Lepton",true,false);
	
	return 0;
}
extern "C"

void run(){
	//~ plot();
	//~ plot_exclusive();
	//~ plot_TChiNg_BR();
	//~ plot_htgVeto_FullKnut();
	plot_leptonVeto_FullLepton();
}
