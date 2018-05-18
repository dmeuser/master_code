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
	legE.append(*combined_exp,"Combined (exp)","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_T5Wg_combine_ST_Lepton",true,false);
	
	return 0;
}

int plot_CharginoBR(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_CharginoBR_inclusiv.root","read");
	TFile file_2("../input/limits/limits_CharginoBR_lepton.root","read");
	TFile file_3("../input/limits/limits_CharginoBR_leptonVeto_leptonFull_htgSTLeptonCleaned.root","read");
	TFile file_4("../input/limits/limits_CharginoBR_htg.root","read");
	
	TGraph *inclusiv_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *lepton_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *st_lepton_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_4.Get("gr_expC_sm");
	
	lepton_exp->SetTitle(";BR(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow W^{#pm}#tilde{G}) (%);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	lepton_exp->SetLineColor(kBlue);
	lepton_exp->SetLineWidth(3);
	lepton_exp->SetLineStyle(2);
	lepton_exp->Draw();
	lepton_exp->GetXaxis()->SetRangeUser(0,100);
	lepton_exp->GetYaxis()->SetRangeUser(300,1300);
	lepton_exp->Draw("same");
	can.Update();
	
	inclusiv_exp->SetLineColor(kGreen);
	inclusiv_exp->SetLineWidth(3);
	inclusiv_exp->SetLineStyle(2);
	inclusiv_exp->Draw("same");
	
	st_lepton_exp->SetLineColor(kMagenta);
	st_lepton_exp->SetLineWidth(3);
	st_lepton_exp->SetLineStyle(2);
	st_lepton_exp->Draw("same");
	
	htg_exp->SetLineColor(kRed);
	htg_exp->SetLineWidth(3);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*inclusiv_exp,"Photon+ST (exp)","l");
	legE.append(*lepton_exp,"Photon+Lepton (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*st_lepton_exp,"Combined (exp)","l");
	TLegend leg=legE.buildLegend(.3,.3,0.9,.5,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_compare_CharginoBR",true,false);
	
	return 0;
}

int plot_leptonVeto_FullLepton_corrCheck(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_T5Wg_leptonVeto_leptonFull.root","read");
	TFile file_2("../input/limits/limits_T5Wg_leptonVeto_leptonFull_corrRare.root","read");
	TFile file_3("../input/limits/limits_T5Wg_leptonVeto_leptonFull_corrVg.root","read");
	TFile file_4("../input/limits/limits_T5Wg_leptonVeto_leptonFull_corrEfake.root","read");
	TFile file_5("../input/limits/limits_T5Wg_leptonVeto_leptonFull_corrJetFake.root","read");
	
	TGraph *initial = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *Rare = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *Vg = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *Efake = (TGraph*) file_4.Get("gr_expC_sm");
	TGraph *JetFake = (TGraph*) file_5.Get("gr_expC_sm");
	TGraph *diag = new TGraph();
	
	diag->SetPoint(1,500,500);
	diag->SetPoint(2,2200,2200);
	
	initial->SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	initial->SetLineColor(kGreen);
	initial->SetLineWidth(3);
	initial->SetLineStyle(2);
	initial->Draw();
	
	Rare->SetLineColor(kBlue);
	Rare->SetLineWidth(3);
	Rare->SetLineStyle(2);
	Rare->Draw("same");
	
	Vg->SetLineColor(kRed);
	Vg->SetLineWidth(3);
	Vg->SetLineStyle(2);
	Vg->Draw("same");
	
	Efake->SetLineColor(kMagenta);
	Efake->SetLineWidth(3);
	Efake->SetLineStyle(2);
	Efake->Draw("same");
	
	JetFake->SetLineColor(kOrange);
	JetFake->SetLineWidth(3);
	JetFake->SetLineStyle(2);
	JetFake->Draw("same");
	
	diag->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*initial,"Without correlation","l");
	legE.append(*Vg,"Vg correlated","l");
	legE.append(*Efake,"Efake correlated","l");
	legE.append(*Rare,"Rare correlated","l");
	legE.append(*JetFake,"JetFake correlated","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_T5Wg_correlations_ST_Lepton",true,false);
	
	return 0;
}

int plot_T5Wg_htgHighLeptonVeto_leptonFull_htgHigh(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_T5Wg_inclusiv.root","read");
	TFile file_2("../input/limits/limits_T5Wg_lepton.root","read");
	TFile file_3("../input/limits/limits_T5Wg_htg.root","read");
	TFile file_4("../input/limits/limits_T5Wg_htgHighLeptonVeto_leptonFull_htgHigh.root","read");
	
	TGraph *st_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *lepton_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *combined_exp = (TGraph*) file_4.Get("gr_expC_sm");
	TGraph *diag = new TGraph();
	
	diag->SetPoint(1,500,500);
	diag->SetPoint(2,2200,2200);
	
	st_exp->SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	st_exp->SetLineColor(kGreen);
	st_exp->SetLineWidth(3);
	st_exp->SetLineStyle(2);
	st_exp->Draw();
	
	lepton_exp->SetLineColor(kMagenta);
	lepton_exp->SetLineWidth(3);
	lepton_exp->SetLineStyle(2);
	lepton_exp->Draw("same");
	
	htg_exp->SetLineColor(kRed);
	htg_exp->SetLineWidth(3);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	
	combined_exp->SetLineColor(kBlue);
	combined_exp->SetLineWidth(3);
	combined_exp->SetLineStyle(2);
	combined_exp->Draw("same");
	
	diag->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*st_exp,"Photon+ST (exp)","l");
	legE.append(*lepton_exp,"Photon+Lepton (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*combined_exp,"Combined (exp)","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_T5Wg_combine_ST_Lepton_HTG",true,false);
	
	return 0;
}

int plot_T6gg_htgHighVeto_htgHigh(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_T6gg_inclusiv.root","read");
	TFile file_2("../input/limits/limits_T6gg_htg.root","read");
	TFile file_3("../input/limits/limits_T6gg_htgHighVeto_htgHigh.root","read");
	
	TGraph *st_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *combined_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *diag = new TGraph();
	
	diag->SetPoint(1,500,500);
	diag->SetPoint(2,2200,2200);
	
	st_exp->SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
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
	
	
	diag->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*st_exp,"Photon+ST (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*combined_exp,"Combined (exp)","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_T6gg_combine_ST_HTG",true,false);
	
	return 0;
}

int plot_T6Wg_htgHighVeto_htgHigh(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_T6Wg_inclusiv.root","read");
	TFile file_2("../input/limits/limits_T6Wg_htg.root","read");
	TFile file_3("../input/limits/limits_T6Wg_htgHighVeto_htgHigh.root","read");
	
	TGraph *st_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *combined_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *diag = new TGraph();
	
	diag->SetPoint(1,500,500);
	diag->SetPoint(2,2200,2200);
	
	st_exp->SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
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
	
	
	diag->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*st_exp,"Photon+ST (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*combined_exp,"Combined (exp)","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_T6Wg_combine_ST_HTG",true,false);
	
	return 0;
}

int plot_T5gg_htgHighVeto_htgHigh(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_T5gg_inclusiv.root","read");
	TFile file_2("../input/limits/limits_T5gg_htg.root","read");
	TFile file_3("../input/limits/limits_T5gg_htgHighVeto_htgHigh.root","read");
	
	TGraph *st_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *combined_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *diag = new TGraph();
	
	diag->SetPoint(1,500,500);
	diag->SetPoint(2,2200,2200);
	
	st_exp->SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
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
	
	
	diag->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*st_exp,"Photon+ST (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*combined_exp,"Combined (exp)","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_T5gg_combine_ST_HTG",true,false);
	
	return 0;
}

int plot_T5Wg_htgHighVeto_htgHigh(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_T5Wg_inclusiv.root","read");
	TFile file_2("../input/limits/limits_T5Wg_htg.root","read");
	TFile file_4("../input/limits/limits_T5Wg_htgHighVeto_highHtg.root","read");
	
	TGraph *st_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_2.Get("gr_expC_sm");
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
	
	combined_split_exp->SetLineColor(kRed);
	combined_split_exp->SetLineWidth(3);
	combined_split_exp->SetLineStyle(2);
	combined_split_exp->Draw("same");
	
	diag->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*st_exp,"Photon+ST (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*combined_split_exp,"Combined (exp)","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_T5Wg_combine_ST_HTG_2",true,false);
	
	return 0;
}

int plot_GGM_htgHighVeto_htgHigh(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_GGM_M1_M2_inclusiv.root","read");
	TFile file_2("../input/limits/limits_GGM_M1_M2_htg.root","read");
	TFile file_4("../input/limits/limits_GGM_M1_M2_htgHighVeto_htgHigh.root","read");
	
	
	TGraph *st_exp = (TGraph*) file_1.Get("gr_expC");
	TGraph *st_exp_extra = (TGraph*) file_1.Get("multicont_gr_8");
	TGraph *htg_exp = (TGraph*) file_2.Get("gr_expC");
	TGraph *combined_exp = (TGraph*) file_4.Get("gr_expC");
	TGraph *combined_exp_extra = (TGraph*) file_4.Get("multicont_gr_8");

	st_exp->SetTitle(";M1 (GeV);M2 (GeV)");
	
	st_exp->SetLineColor(kGreen);
	st_exp->SetLineWidth(3);
	st_exp->SetLineStyle(2);
	st_exp->Draw();
	st_exp_extra->SetLineColor(kGreen);
	st_exp_extra->SetLineWidth(3);
	st_exp_extra->SetLineStyle(2);
	st_exp_extra->Draw("same");
	
	htg_exp->SetLineColor(kBlue);
	htg_exp->SetLineWidth(3);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	
	combined_exp->SetLineColor(kRed);
	combined_exp->SetLineWidth(3);
	combined_exp->SetLineStyle(2);
	combined_exp->Draw("same");
	combined_exp_extra->SetLineColor(kRed);
	combined_exp_extra->SetLineWidth(3);
	combined_exp_extra->SetLineStyle(2);
	combined_exp_extra->Draw("same");
 
	gfx::LegendEntries legE;
	legE.append(*st_exp,"Photon+ST (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*combined_exp,"Combined (exp)","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_GGM_combine_ST_HTG",true,false);
	
	return 0;
}

int plot_T6gg_diffCombi(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_T6gg_inclusiv_stVeto.root","read");
	TFile file_2("../input/limits/limits_T6gg_htgVeto_htgFull.root","read");
	TFile file_3("../input/limits/limits_T6gg_htgHighVeto_htgHigh.root","read");
	TFile file_4("../input/limits/limits_T6gg_inclusiv.root","read");
	TFile file_5("../input/limits/limits_T6gg_htg.root","read");
	
	TGraph *stVeto_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *htgVeto_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *htgHigh_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *st_exp = (TGraph*) file_4.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_5.Get("gr_expC_sm");
	TGraph *diag = new TGraph();
	
	diag->SetPoint(1,500,500);
	diag->SetPoint(2,2200,2200);
	
	stVeto_exp->SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	stVeto_exp->SetLineColor(kMagenta);
	stVeto_exp->SetLineWidth(3);
	stVeto_exp->SetLineStyle(2);
	stVeto_exp->Draw();
	
	htgVeto_exp->SetLineColor(kOrange);
	htgVeto_exp->SetLineWidth(3);
	htgVeto_exp->SetLineStyle(2);
	htgVeto_exp->Draw("same");
	
	htgHigh_exp->SetLineColor(kBlack);
	htgHigh_exp->SetLineWidth(3);
	htgHigh_exp->SetLineStyle(2);
	htgHigh_exp->Draw("same");
	
	st_exp->SetLineColor(kGreen);
	st_exp->SetLineWidth(3);
	st_exp->SetLineStyle(2);
	st_exp->Draw("same");
	
	htg_exp->SetLineColor(kRed);
	htg_exp->SetLineWidth(3);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	
	
	diag->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*htgHigh_exp,"ST:highHTG veto HTG:highHTG bins (exp)","l");
	legE.append(*stVeto_exp,"ST:full HTG:ST veto (exp)","l");
	legE.append(*htgVeto_exp,"ST:HTG veto HTG:full (exp)","l");
	legE.append(*st_exp,"Photon+ST (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_T6gg_diffCombi",true,false);
	
	return 0;
}

int plot_TChiNg_BR_combi(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_TChiNg_BR_inclusiv.root","read");
	TFile file_2("../input/limits/limits_TChiNg_BR_htg.root","read");
	TFile file_3("../input/limits/limits_TChiNg_BR_stFull_htgSTcleaned.root","read");
	
	TGraph *inclusiv_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *combi_exp = (TGraph*) file_3.Get("gr_expC_sm");
	
	inclusiv_exp->SetTitle(";BR(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #gamma) (%);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	inclusiv_exp->SetLineColor(kGreen);
	inclusiv_exp->SetLineWidth(3);
	inclusiv_exp->SetLineStyle(2);
	inclusiv_exp->Draw();
	
	htg_exp->SetLineColor(kRed);
	htg_exp->SetLineWidth(3);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	
	combi_exp->SetLineColor(kBlue);
	combi_exp->SetLineWidth(3);
	combi_exp->SetLineStyle(2);
	combi_exp->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*inclusiv_exp,"Photon+ST (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*combi_exp,"ST:full HTG:ST veto (exp)","l");
	TLegend leg=legE.buildLegend(.7,.2,0.9,.5,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_TChiNg_BR_combi",true,false);
	
	return 0;
}

int plot_CharginoBRstrongN1700(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_CharginoBRstrongN1700_inclusiv.root","read");
	TFile file_2("../input/limits/limits_CharginoBRstrongN1700_htg.root","read");
	TFile file_3("../input/limits/limits_CharginoBRstrongN1700_htgHighVeto_htgHigh.root","read");
	
	TGraph *inclusiv_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *combi_exp = (TGraph*) file_3.Get("gr_expC_sm");
	
	combi_exp->SetTitle(";BR(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow W^{#pm}#tilde{G}) (%);m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV)");
	
	combi_exp->SetLineColor(kBlue);
	combi_exp->SetLineWidth(3);
	combi_exp->SetLineStyle(2);
	combi_exp->Draw();
	
	inclusiv_exp->SetLineColor(kGreen);
	inclusiv_exp->SetLineWidth(3);
	inclusiv_exp->SetLineStyle(2);
	inclusiv_exp->Draw("same");
	
	htg_exp->SetLineColor(kRed);
	htg_exp->SetLineWidth(3);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*inclusiv_exp,"Photon+ST (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*combi_exp,"ST:highHTGVeto HTG:highHTG (exp)","l");
	TLegend leg=legE.buildLegend(.7,.2,0.9,.5,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_CharginoBRstrongN1700",true,false);
	
	return 0;
}

int plot_CharginoBRstrongG1950(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_CharginoBRstrongG1950_inclusiv.root","read");
	TFile file_2("../input/limits/limits_CharginoBRstrongG1950_htg.root","read");
	TFile file_3("../input/limits/limits_CharginoBRstrongG1950_htgHighVeto_htgHigh.root","read");
	
	TGraph *inclusiv_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *combi_exp = (TGraph*) file_3.Get("gr_expC_sm");
	
	combi_exp->SetTitle(";BR(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow W^{#pm}#tilde{G}) (%);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	combi_exp->SetLineColor(kBlue);
	combi_exp->SetLineWidth(3);
	combi_exp->SetLineStyle(2);
	combi_exp->Draw();
	
	inclusiv_exp->SetLineColor(kGreen);
	inclusiv_exp->SetLineWidth(3);
	inclusiv_exp->SetLineStyle(2);
	inclusiv_exp->Draw("same");
	
	htg_exp->SetLineColor(kRed);
	htg_exp->SetLineWidth(3);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*inclusiv_exp,"Photon+ST (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*combi_exp,"ST:highHTGVeto HTG:highHTG (exp)","l");
	TLegend leg=legE.buildLegend(.7,.2,0.9,.5,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_CharginoBRstrongG1950",true,false);
	
	return 0;
}
extern "C"

void run(){
	//~ plot();
	//~ plot_exclusive();
	//~ plot_TChiNg_BR();
	//~ plot_htgVeto_FullKnut();
	//~ plot_leptonVeto_FullLepton();
	//~ plot_leptonVeto_FullLepton_corrCheck();
	//~ plot_CharginoBR();
	//~ plot_T5Wg_htgHighLeptonVeto_leptonFull_htgHigh();
	//~ plot_T6gg_htgHighVeto_htgHigh();
	//~ plot_T6Wg_htgHighVeto_htgHigh();
	//~ plot_T5gg_htgHighVeto_htgHigh();
	//~ plot_T5Wg_htgHighVeto_htgHigh();
	//~ plot_GGM_htgHighVeto_htgHigh();
	//~ plot_T6gg_diffCombi();
	//~ plot_TChiNg_BR_combi();
	//~ plot_CharginoBRstrongN1700();
	plot_CharginoBRstrongG1950();
}
