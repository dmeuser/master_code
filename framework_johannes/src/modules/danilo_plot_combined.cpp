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
	
	//~ leptonAna_exp->SetLineColor(kBlack);
	//~ leptonAna_exp->SetLineWidth(3);
	//~ leptonAna_exp->SetLineStyle(2);
	//~ leptonAna_exp->Draw("same");
	
	diag->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*johannes_exp,"Inclusiv (exp)","l");
	legE.append(*exclusiv_exp,"Exclusiv (exp)","l");
	legE.append(*leptonVeto_exp,"LeptonVeto (exp)","l");
	legE.append(*htgVeto_exp,"HTgVeto (exp)","l");
	legE.append(*diphotonVeto_exp,"DiphotonVeto (exp)","l");
	//~ legE.append(*leptonAna_exp,"LeptonAnalysis approx. (exp)","l");
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
	TH2F axis("","",26,0,100,26,300,1300);
	axis.SetTitle(";BR(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #gamma + #tile{G}) (%);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	axis.Draw("axis");
	axis.SetStats(0);
		
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
	
	axis.Draw("axis same");
	
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
	
	TFile file_1("../input/limits/limits_CharginoBR_inclusivNN.root","read");
	TFile file_2("../input/limits/limits_CharginoBR_lepton.root","read");
	TFile file_3("../input/limits/limits_CharginoBR_allCombined_htgHighNN.root","read");
	TFile file_4("../input/limits/limits_CharginoBR_htgNN.root","read");
	TFile file_5("../input/limits/limits_CharginoBR_diphoton.root","read");
	
	TGraph *inclusiv_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *lepton_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *st_lepton_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_4.Get("gr_expC_sm");
	TGraph *diphoton_exp = (TGraph*) file_5.Get("gr_expC_sm");
	
	lepton_exp->SetTitle(";BR(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow W^{#pm}#tilde{G}) (%);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	lepton_exp->SetLineColor(kMagenta);
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
	
	st_lepton_exp->SetLineColor(kBlue);
	st_lepton_exp->SetLineWidth(3);
	st_lepton_exp->SetLineStyle(2);
	st_lepton_exp->Draw("same");
	
	htg_exp->SetLineColor(kRed+1);
	htg_exp->SetLineWidth(3);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	
	diphoton_exp->SetLineColor(kCyan);
	diphoton_exp->SetLineWidth(3);
	diphoton_exp->SetLineStyle(2);
	diphoton_exp->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*inclusiv_exp,"Photon+ST (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*lepton_exp,"Photon+Lepton (exp)","l");
	legE.append(*diphoton_exp,"Diphoton (exp)","l");
	legE.append(*st_lepton_exp,"Combined (exp)","l");
	TLegend leg=legE.buildLegend(.3,.3,0.9,.5,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_combi_CharginoBR",true,false);
	
	return 0;
}

int plot_CharginoBR_C1C1(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	//~ TFile file_1("../input/limits/limits_CharginoBR_C1C1_inclusivNN.root","read");
	//~ TFile file_2("../input/limits/limits_CharginoBR_C1C1_lepton.root","read");
	TFile file_1("../input/limits/limits_CharginoBR_C1C1_inclusivFinal.root","read");
	TFile file_2("../input/limits/limits_CharginoBR_C1C1_lepton_final.root","read");
	//~ TFile file_3("../input/limits/limits_CharginoBR_C1C1_allCombined_htgHighNN.root","read");
	TFile file_3("../input/limits/limits_CharginoBR_C1C1_allCombined_final.root","read");
	//~ TFile file_4("../input/limits/limits_CharginoBR_C1C1_htgNN.root","read");
	TFile file_4("../input/limits/limits_CharginoBR_C1C1_htgFinal.root","read");
	//~ TFile file_5("../input/limits/limits_CharginoBR_C1C1_diphoton.root","read");
	TFile file_5("../input/limits/limits_CharginoBR_C1C1_diphoton_final.root","read");
	
	TGraph *inclusiv_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *lepton_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *st_lepton_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_4.Get("multicont_gr_2_sm");
	TGraph *diphoton_exp = (TGraph*) file_5.Get("gr_expC_sm");
	TGraph *inclusiv_obs = (TGraph*) file_1.Get("gr_obsC_sm");
	TGraph *lepton_obs = (TGraph*) file_2.Get("gr_obsC_sm");
	TGraph *st_lepton_obs = (TGraph*) file_3.Get("gr_obsC_sm");
	TGraph *htg_obs = (TGraph*) file_4.Get("gr_obsC_sm");
	TGraph *diphoton_obs = (TGraph*) file_5.Get("gr_obsC_sm");
	
	TH2F axis("","",26,0,100,26,300,1300);
	axis.SetTitle(";BF(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #tilde{#chi}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}(#gamma#tilde{G}) + soft) (%);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	axis.Draw("axis");
	axis.SetStats(0);
		
	lepton_exp->SetLineColor(kMagenta);
	lepton_exp->SetLineWidth(2);
	lepton_exp->SetLineStyle(2);
	lepton_exp->Draw("same");
	lepton_obs->SetLineColor(kMagenta);
	lepton_obs->SetLineWidth(2);
	lepton_obs->SetLineStyle(1);
	lepton_obs->Draw("same");
	
	inclusiv_exp->SetLineColor(kGreen+1);
	inclusiv_exp->SetLineWidth(2);
	inclusiv_exp->SetLineStyle(2);
	inclusiv_exp->Draw("same");
	inclusiv_obs->SetLineColor(kGreen+1);
	inclusiv_obs->SetLineWidth(2);
	inclusiv_obs->SetLineStyle(1);
	inclusiv_obs->Draw("same");
	
	st_lepton_exp->SetLineColor(kBlue);
	st_lepton_exp->SetLineWidth(2);
	st_lepton_exp->SetLineStyle(2);
	st_lepton_exp->Draw("same");
	st_lepton_obs->SetLineColor(kBlue);
	st_lepton_obs->SetLineWidth(2);
	st_lepton_obs->SetLineStyle(1);
	st_lepton_obs->Draw("same");
	
	//~ htg_exp->SetLineColor(kRed+1);
	//~ htg_exp->SetLineWidth(2);
	//~ htg_exp->SetLineStyle(2);
	//~ htg_exp->Draw("same");
	//~ htg_obs->SetLineColor(kRed+1);
	//~ htg_obs->SetLineWidth(2);
	//~ htg_obs->SetLineStyle(1);
	//~ htg_obs->Draw("same");
	
	diphoton_exp->SetLineColor(kOrange-3);
	diphoton_exp->SetLineWidth(2);
	diphoton_exp->SetLineStyle(2);
	diphoton_exp->Draw("same");
	diphoton_obs->SetLineColor(kOrange-3);
	diphoton_obs->SetLineWidth(2);
	diphoton_obs->SetLineStyle(1);
	diphoton_obs->Draw("same");
	
	TGraph *exp = new TGraph();
	TGraph *obs = new TGraph();
	exp->SetLineWidth(2);
	obs->SetLineWidth(2);
	exp->SetLineStyle(2);
	
	axis.Draw("axis same");
	
	gfx::LegendEntries legE;
	legE.append(*inclusiv_obs,"Photon+S_{#scale[.8]{T}}^{#scale[.8]{#gamma}}","l");
	//~ legE.append(*htg_obs,"Photon+H_{#scale[.8]{T}}^{#scale[.8]{#gamma}}","l");
	legE.append(*lepton_obs,"Photon+Lepton","l");
	legE.append(*diphoton_obs,"Diphoton","l");
	legE.append(*st_lepton_obs,"Combination","l");
	legE.prepend(*exp,"Expected","l");
	legE.prepend(*obs,"Observed","l");
	TLegend leg=legE.buildLegend(.2,.7,0.7,.9,2);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_combi_CharginoBR_C1C1",true,false);
	
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
	
	//~ TFile file_1("../input/limits/limits_T5Wg_inclusivNN.root","read");
	//~ TFile file_2("../input/limits/limits_T5Wg_lepton.root","read");
	//~ TFile file_3("../input/limits/limits_T5Wg_htgNN.root","read");
	//~ TFile file_4("../input/limits/limits_T5Wg_htgHighLeptonVetoNN_leptonFull_htgHighNN.root","read");
	TFile file_1("../input/limits/limits_T5Wg_inclusivFinal.root","read");
	TFile file_2("../input/limits/limits_T5Wg_lepton_final.root","read");
	TFile file_3("../input/limits/limits_T5Wg_htgFinal.root","read");
	TFile file_4("../input/limits/limits_T5Wg_ST_HTG_lepton_final.root","read");
	
	TGraph *st_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *lepton_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *combined_exp = (TGraph*) file_4.Get("gr_expC_sm");
	TGraph *st_obs = (TGraph*) file_1.Get("gr_obsC_sm");
	TGraph *lepton_obs = (TGraph*) file_2.Get("gr_obsC_sm");
	TGraph *htg_obs = (TGraph*) file_3.Get("gr_obsC_sm");
	TGraph *combined_obs = (TGraph*) file_4.Get("gr_obsC_sm");
	TGraph *diag = new TGraph();
	
	diag->SetPoint(1,500,490);
	diag->SetPoint(2,2200,2190);
	
	TH2F axis("","",27,1400,2100,27,0,2400);
	axis.SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	axis.Draw("axis");
	axis.SetStats(0);
	
	st_exp->SetLineColor(kGreen+1);
	st_exp->SetLineWidth(2);
	st_exp->SetLineStyle(2);
	st_exp->Draw("same");
	st_obs->SetLineColor(kGreen+1);
	st_obs->SetLineWidth(2);
	st_obs->SetLineStyle(1);
	st_obs->Draw("same");
	
	int points = lepton_exp->GetN();
	int diff = points-lepton_obs->GetN();
	for (int i=points; i>=points-5; i--) {
		lepton_exp->RemovePoint(i);
		lepton_obs->RemovePoint(i-diff);
	}
	
	lepton_exp->SetLineColor(kMagenta);
	lepton_exp->SetLineWidth(2);
	lepton_exp->SetLineStyle(2);
	lepton_exp->Draw("same");
	lepton_obs->SetLineColor(kMagenta);
	lepton_obs->SetLineWidth(2);
	lepton_obs->SetLineStyle(1);
	lepton_obs->Draw("same");
	
	points = htg_exp->GetN();
	diff = points-htg_obs->GetN();
	for (int i=points; i>=points-5; i--) {
		htg_exp->RemovePoint(i);
		htg_obs->RemovePoint(i-diff);
	}
	
	htg_exp->SetLineColor(kRed+1);
	htg_exp->SetLineWidth(2);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	htg_obs->SetLineColor(kRed+1);
	htg_obs->SetLineWidth(2);
	htg_obs->SetLineStyle(1);
	htg_obs->Draw("same");
	
	combined_exp->SetLineColor(kBlue);
	combined_exp->SetLineWidth(2);
	combined_exp->SetLineStyle(2);
	combined_exp->Draw("same");
	combined_obs->SetLineColor(kBlue);
	combined_obs->SetLineWidth(2);
	combined_obs->SetLineStyle(1);
	combined_obs->Draw("same");
	
	diag->Draw("same F");
	diag->Draw("same");
	
	axis.Draw("axis same");
	
	TGraph *exp = new TGraph();
	TGraph *obs = new TGraph();
	exp->SetLineWidth(2);
	obs->SetLineWidth(2);
	exp->SetLineStyle(2);
	
	gfx::LegendEntries legE;
	legE.append(*st_obs,"Photon+S_{#scale[.8]{T}}^{#scale[.8]{#gamma}}","l");
	legE.append(*lepton_obs,"Photon+Lepton","l");
	legE.append(*htg_obs,"Photon+H_{#scale[.8]{T}}^{#scale[.8]{#gamma}}","l");
	legE.append(*combined_obs,"Combination","l");
	legE.prepend(*exp,"Expected","l");
	legE.prepend(*obs,"Observed","l");
	TLegend leg=legE.buildLegend(.2,.7,0.7,.9,2);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_T5Wg_combine_ST_Lepton_HTG",true,false);
	
	return 0;
}

int plot_T5Wg_diphotonCase(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_T5Wg_ST_HTG_lepton_final.root","read");
	TFile file_2("../input/limits/limits_T5Wg_allCombined_withDiFinal.root","read");
	
	TGraph *without_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *with_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *diag = new TGraph();
	
	diag->SetPoint(1,500,490);
	diag->SetPoint(2,2200,2190);
	
	TH2F axis("","",27,1400,2100,27,0,2400);
	axis.SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	axis.Draw("axis");
	axis.SetStats(0);
	
	without_exp->SetLineColor(kBlue);
	without_exp->SetLineWidth(3);
	without_exp->SetLineStyle(2);
	without_exp->Draw("same");
	
	int points = with_exp->GetN();
	for (int i=points; i>=points-5; i--) {
		with_exp->RemovePoint(i);
	}
	
	with_exp->SetLineColor(kGray+1);
	with_exp->SetLineWidth(3);
	with_exp->SetLineStyle(2);
	with_exp->Draw("same");
	
	diag->Draw("same F");
	diag->Draw("same");
	
	axis.Draw("axis same");
	
	gfx::LegendEntries legE;
	legE.append(*without_exp,"Without Diphoton","l");
	legE.append(*with_exp,"With Diphoton","l");
	TLegend leg=legE.buildLegend(.2,.7,0.7,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_T5Wg_diphotonCase",true,false);
	
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
	
	htgHigh_exp->SetLineColor(kBlue);
	htgHigh_exp->SetLineWidth(3);
	htgHigh_exp->SetLineStyle(2);
	htgHigh_exp->Draw("same");
	
	st_exp->SetLineColor(kGreen);
	st_exp->SetLineWidth(3);
	st_exp->SetLineStyle(2);
	st_exp->Draw("same");
	
	htg_exp->SetLineColor(kRed+1);
	htg_exp->SetLineWidth(3);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	
	
	diag->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*htgHigh_exp,"Strategy 3 (exp)","l");
	legE.append(*htgVeto_exp,"Strategy 2 (exp)","l");
	legE.append(*stVeto_exp,"Strategy 1 (exp)","l");
	legE.append(*st_exp,"Photon+ST (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_T6gg_diffCombi",true,false);
	
	return 0;
}


int plot_T5Wg_diffCombi(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_T5Wg_inclusive_STcleaned.root","read");
	TFile file_2("../input/limits/limits_T5Wg_htgVeto_KnutFull.root","read");
	TFile file_3("../input/limits/limits_T5Wg_htgHighVeto_highHtg.root","read");

	
	TGraph *stVeto_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *htgVeto_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *htgHigh_exp = (TGraph*) file_3.Get("gr_expC_sm");

	TGraph *diag = new TGraph();
	
	diag->SetPoint(1,500,500);
	diag->SetPoint(2,2200,2200);
	
	stVeto_exp->SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	int points = stVeto_exp->GetN();
	
	for (int i=points; i>=points-6; i--) {
		stVeto_exp->RemovePoint(i);
	}
	stVeto_exp->SetLineColor(kGreen+1);
	stVeto_exp->SetLineWidth(3);
	stVeto_exp->SetLineStyle(2);
	stVeto_exp->Draw();
	
	htgVeto_exp->SetLineColor(kRed+1);
	htgVeto_exp->SetLineWidth(3);
	htgVeto_exp->SetLineStyle(2);
	htgVeto_exp->Draw("same");
	
	htgHigh_exp->SetLineColor(kBlue);
	htgHigh_exp->SetLineWidth(3);
	htgHigh_exp->SetLineStyle(2);
	htgHigh_exp->Draw("same");
	
	
	diag->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*htgHigh_exp,"Strategy 3 (exp)","l");
	legE.append(*htgVeto_exp,"Strategy 2 (exp)","l");
	legE.append(*stVeto_exp,"Strategy 1 (exp)","l");
	TLegend leg=legE.buildLegend(.6,.75,0.92,.9,1);
	leg.SetTextSize(0.033);
	leg.Draw();
	
	saver.save(can,"Limits_T5Wg_diffCombi",true,false);
	
	return 0;
}


int plot_TChiNg_BR_diffcombi(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_TChiNg_BR_inclusiv.root","read");
	TFile file_2("../input/limits/limits_TChiNg_BR_htg.root","read");
	TFile file_3("../input/limits/limits_TChiNg_BR_stFull_htgSTcleaned.root","read");
	TFile file_4("../input/limits/limits_TChiNg_BR_htgHighVeto_htgHigh.root","read");
	TFile file_5("../input/limits/limits_TChiNg_BR_htgVeto_htg.root","read");
	
	TGraph *inclusiv_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *combi_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *combi2_exp = (TGraph*) file_4.Get("gr_expC_sm");
	TGraph *combi3_exp = (TGraph*) file_5.Get("gr_expC_sm");
	
	inclusiv_exp->SetTitle(";BR(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #gamma) (%);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	inclusiv_exp->SetLineColor(kGreen);
	inclusiv_exp->SetLineWidth(3);
	inclusiv_exp->SetLineStyle(2);
	inclusiv_exp->Draw();
	
	htg_exp->SetLineColor(kRed+1);
	htg_exp->SetLineWidth(3);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	
	combi_exp->SetLineColor(kBlue);
	combi_exp->SetLineWidth(3);
	combi_exp->SetLineStyle(2);
	combi_exp->Draw("same");
	
	combi2_exp->SetLineColor(kMagenta);
	combi2_exp->SetLineWidth(3);
	combi2_exp->SetLineStyle(2);
	combi2_exp->Draw("same");
	
	combi3_exp->SetLineColor(kBlack);
	combi3_exp->SetLineWidth(3);
	combi3_exp->SetLineStyle(2);
	combi3_exp->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*inclusiv_exp,"Photon+ST (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*combi_exp,"ST:full HTG:ST veto (exp)","l");
	legE.append(*combi2_exp,"ST:highHTG veto HTG:highHTG bins (exp)","l");
	legE.append(*combi3_exp,"ST:HTG veto HTG:full (exp)","l");
	TLegend leg=legE.buildLegend(.7,.2,0.9,.5,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_TChiNg_BR_diffcombi",true,false);
	
	return 0;
}

int plot_TChiNg_BR_combi(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	//~ TFile file_1("../input/limits/limits_TChiNg_BR_inclusivNN.root","read");
	//~ TFile file_2("../input/limits/limits_TChiNg_BR_htgNN.root","read");
	//~ TFile file_3("../input/limits/limits_TChiNg_BR_lepton.root","read");
	TFile file_1("../input/limits/limits_TChiNg_BR_inclusivFinal.root","read");
	TFile file_2("../input/limits/limits_TChiNg_BR_htgFinal.root","read");
	TFile file_3("../input/limits/limits_TChiNg_BR_lepton_final.root","read");
	//~ TFile file_4("../input/limits/limits_TChiNg_BR_diphoton.root","read");
	TFile file_4("../input/limits/limits_TChiNg_BR_diphoton_final.root","read");
	//~ TFile file_5("../input/limits/limits_TChiNg_BR_allCombined_highHtgNN.root","read");
	TFile file_5("../input/limits/limits_TChiNg_BR_allCombined_final.root","read");
	
	TGraph *inclusiv_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *lepton_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *diphoton_exp = (TGraph*) file_4.Get("gr_expC_sm");
	TGraph *combi_exp = (TGraph*) file_5.Get("gr_expC_sm");
	TGraph *inclusiv_obs = (TGraph*) file_1.Get("gr_obsC_sm");
	TGraph *htg_obs = (TGraph*) file_2.Get("gr_obsC_sm");
	TGraph *lepton_obs = (TGraph*) file_3.Get("gr_obsC_sm");
	TGraph *diphoton_obs = (TGraph*) file_4.Get("gr_obsC_sm");
	TGraph *combi_obs = (TGraph*) file_5.Get("gr_obsC_sm");
	
	TH2F axis("","",26,0,100,26,300,1300);
	axis.SetTitle(";BF(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #gamma + #tilde{G}) (%);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	axis.Draw("axis");
	axis.SetStats(0);
	
	inclusiv_exp->SetLineColor(kGreen+1);
	inclusiv_exp->SetLineWidth(2);
	inclusiv_exp->SetLineStyle(2);
	inclusiv_exp->Draw("same");
	inclusiv_obs->SetLineColor(kGreen+1);
	inclusiv_obs->SetLineWidth(2);
	inclusiv_obs->SetLineStyle(1);
	inclusiv_obs->Draw("same");
	
	htg_exp->SetLineColor(kRed+1);
	htg_exp->SetLineWidth(2);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	htg_obs->SetLineColor(kRed+1);
	htg_obs->SetLineWidth(2);
	htg_obs->SetLineStyle(1);
	htg_obs->Draw("same");
	
	lepton_exp->SetLineColor(kMagenta);
	lepton_exp->SetLineWidth(2);
	lepton_exp->SetLineStyle(2);
	lepton_exp->Draw("same");
	lepton_obs->SetLineColor(kMagenta);
	lepton_obs->SetLineWidth(2);
	lepton_obs->SetLineStyle(1);
	lepton_obs->Draw("same");
	
	diphoton_exp->SetLineColor(kOrange-3);
	diphoton_exp->SetLineWidth(2);
	diphoton_exp->SetLineStyle(2);
	diphoton_exp->Draw("same");
	diphoton_obs->SetLineColor(kOrange-3);
	diphoton_obs->SetLineWidth(2);
	diphoton_obs->SetLineStyle(1);
	diphoton_obs->Draw("same");
	
	combi_exp->SetLineColor(kBlue);
	combi_exp->SetLineWidth(2);
	combi_exp->SetLineStyle(2);
	combi_exp->Draw("same");
	combi_obs->SetLineColor(kBlue);
	combi_obs->SetLineWidth(2);
	combi_obs->SetLineStyle(1);
	combi_obs->Draw("same");
	
	TGraph *exp = new TGraph();
	TGraph *obs = new TGraph();
	exp->SetLineWidth(2);
	obs->SetLineWidth(2);
	exp->SetLineStyle(2);
	
	axis.Draw("axis same");
	
	gfx::LegendEntries legE;
	legE.append(*inclusiv_obs,"Photon+S_{#scale[.8]{T}}^{#scale[.8]{#gamma}}","l");
	legE.append(*htg_obs,"Photon+H_{#scale[.8]{T}}^{#scale[.8]{#gamma}}","l");
	legE.append(*lepton_obs,"Photon+Lepton","l");
	legE.append(*diphoton_obs,"Diphoton","l");
	legE.append(*combi_obs,"Combination","l");
	legE.prepend(*exp,"Expected","l");
	legE.prepend(*obs,"Observed","l");
	TLegend leg=legE.buildLegend(.2,.7,0.7,.9,2);
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
	
	htg_exp->SetLineColor(kRed+1);
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
	
	TFile file_1("../input/limits/limits_CharginoBRstrongG1950_inclusivNN.root","read");
	TFile file_2("../input/limits/limits_CharginoBRstrongG1950_htgNN.root","read");
	TFile file_3("../input/limits/limits_CharginoBRstrongG1950_leptonHighHtgVetoNN_LEPcleanedHighHtgNN_leptonFull.root","read");
	
	TGraph *inclusiv_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *combi_exp = (TGraph*) file_3.Get("gr_expC_sm");
	
	combi_exp->SetTitle(";BR(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow W^{#pm}#tilde{G}) (%);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	combi_exp->SetLineColor(kBlue);
	combi_exp->SetLineWidth(3);
	combi_exp->SetLineStyle(2);
	combi_exp->Draw();
	
	inclusiv_exp->SetLineColor(kGreen+1);
	inclusiv_exp->SetLineWidth(3);
	inclusiv_exp->SetLineStyle(2);
	inclusiv_exp->Draw("same");
	
	htg_exp->SetLineColor(kRed+1);
	htg_exp->SetLineWidth(3);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*inclusiv_exp,"Photon+ST (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*combi_exp,"Combined (exp)","l");
	TLegend leg=legE.buildLegend(.7,.2,0.9,.5,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_CharginoBRstrongG1950",true,false);
	
	return 0;
}

int plot_CharginoBRstrongG1800(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_CharginoBRstrongG1800_inclusivNN.root","read");
	TFile file_2("../input/limits/limits_CharginoBRstrongG1800_htgNN.root","read");
	TFile file_3("../input/limits/limits_CharginoBRstrongG1800_leptonHighHtgVetoNN_LEPcleanedHighHtgNN_fullLepton.root","read");
	
	TGraph *inclusiv_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *combi_exp = (TGraph*) file_3.Get("gr_expC_sm");
	
	combi_exp->SetTitle(";BR(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow W^{#pm}#tilde{G}) (%);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	combi_exp->SetLineColor(kBlue);
	combi_exp->SetLineWidth(3);
	combi_exp->SetLineStyle(2);
	combi_exp->Draw();
	
	inclusiv_exp->SetLineColor(kGreen);
	inclusiv_exp->SetLineWidth(3);
	inclusiv_exp->SetLineStyle(2);
	inclusiv_exp->Draw("same");
	
	htg_exp->SetLineColor(kRed+1);
	htg_exp->SetLineWidth(3);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*inclusiv_exp,"Photon+ST (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*combi_exp,"Combined (exp)","l");
	TLegend leg=legE.buildLegend(.7,.2,0.9,.5,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_CharginoBRstrongG1800",true,false);
	
	return 0;
}

int plot_CharginoBRstrongG1700(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_CharginoBRstrongG1700_inclusivNN.root","read");
	TFile file_2("../input/limits/limits_CharginoBRstrongG1700_htgNN.root","read");
	TFile file_3("../input/limits/limits_CharginoBRstrongG1700_leptonHighHtgVetoNN_LEPcleanedHighHtgNN_leptonFull.root","read");
	
	TGraph *inclusiv_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *combi_exp = (TGraph*) file_3.Get("gr_expC_sm");
	
	combi_exp->SetTitle(";BR(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow W^{#pm}#tilde{G}) (%);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	
	combi_exp->SetLineColor(kBlue);
	combi_exp->SetLineWidth(3);
	combi_exp->SetLineStyle(2);
	combi_exp->Draw();
	
	inclusiv_exp->SetLineColor(kGreen);
	inclusiv_exp->SetLineWidth(3);
	inclusiv_exp->SetLineStyle(2);
	inclusiv_exp->Draw("same");
	
	htg_exp->SetLineColor(kRed+1);
	htg_exp->SetLineWidth(3);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	
	gfx::LegendEntries legE;
	legE.append(*inclusiv_exp,"Photon+ST (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*combi_exp,"Combined (exp)","l");
	TLegend leg=legE.buildLegend(.7,.2,0.9,.5,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_CharginoBRstrongG1700",true,false);
	
	return 0;
}

int plot_GGM_combi(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_GGM_M1_M2_inclusivFinal.root","read");
	TFile file_2("../input/limits/limits_GGM_M1_M2_htgFinal.root","read");
	TFile file_3("../input/limits/limits_GGM_M1_M2_lepton_final.root","read");
	TFile file_4("../input/limits/limits_GGM_M1_M2_diphoton_final.root","read");
	TFile file_5("../input/limits/limits_GGM_M1_M2_allCombined_final.root","read");
	//~ TFile file_6("../input/limits/limits_GGM_M1_M2_allCombined_highHtg.root","read");
	
	TH2F axis("","",26,250,1500,26,250,1500);
	axis.SetTitle(";#it{M}_{1} (GeV);#it{M}_{2} (GeV)");
	axis.Draw("axis");
	axis.SetStats(0);
	
	TGraph diphoton;
	TGraph st;
	TGraph lepton;
	TGraph combi;
	
	//Diphoton Limit
	for(int i=1; i<=2; i++){
		TString temp="multicont_gr_"+std::to_string(i)+"_sm";
		std::cout<<temp<<std::endl;
		TGraph *st_exp = (TGraph*) file_4.Get(temp);
		if (diphoton.GetN()==0) diphoton.SetPoint(diphoton.GetN(),250,250);
		if (i==1) {
			for (int j=0; j<=st_exp->GetN()-1; j++) {
				diphoton.SetPoint(diphoton.GetN(),st_exp->GetX()[j],st_exp->GetY()[j]);
			}
		}
		else {
			for (int j=st_exp->GetN()-1; j>=0; j--) {
				diphoton.SetPoint(diphoton.GetN(),st_exp->GetX()[j],st_exp->GetY()[j]);
			}
			diphoton.SetPoint(diphoton.GetN(),1500,250);
		}
	}
	diphoton.SetFillColorAlpha(kCyan,0.3);
	diphoton.SetLineWidth(3);
	diphoton.SetLineStyle(2);
	diphoton.Draw("same F");
	
	//Photon+ST Limit
	for(int i :{1,3,2}){
		TString temp="multicont_gr_"+std::to_string(i)+"_sm";
		TGraph *st_exp = (TGraph*) file_1.Get(temp);
		std::cout<<temp<<" "<<st_exp->GetN()<<std::endl;
		if (i==1) {
			st.SetPoint(st.GetN(),250,250);
			for (int j=0; j<=st_exp->GetN()-1; j++) {
					st.SetPoint(st.GetN(),st_exp->GetX()[j],st_exp->GetY()[j]);
			}
		}
		else if (i==3) {
			for (int j=0; j<=st_exp->GetN()-1; j++) {
					st.SetPoint(st.GetN(),st_exp->GetX()[j],st_exp->GetY()[j]);
			}
		}
		else {
			for (int j=st_exp->GetN()-1; j>=0; j--) {
				st.SetPoint(st.GetN(),st_exp->GetX()[j],st_exp->GetY()[j]);
			}
			st.SetPoint(st.GetN(),1500,250);
		}
	}
	st.SetFillColorAlpha(kGreen,0.2);
	st.SetLineWidth(3);
	st.SetLineStyle(2);
	st.Draw("same F");
	
	//Lepton Limit
	for(int i=1; i<=2; i++){
		TString temp="multicont_gr_"+std::to_string(i)+"_sm";
		std::cout<<temp<<std::endl;
		TGraph *st_exp = (TGraph*) file_3.Get(temp);
		if (lepton.GetN()==0) lepton.SetPoint(lepton.GetN(),250,250);
		if (i==1) {
			for (int j=0; j<=st_exp->GetN()-1; j++) {
				lepton.SetPoint(lepton.GetN(),st_exp->GetX()[j],st_exp->GetY()[j]);
			}
		}
		else {
			for (int j=st_exp->GetN()-1; j>=0; j--) {
				lepton.SetPoint(lepton.GetN(),st_exp->GetX()[j],st_exp->GetY()[j]);
			}
			lepton.SetPoint(lepton.GetN(),1500,250);
		}
	}
	lepton.SetFillColorAlpha(kOrange,0.2);
	lepton.SetLineWidth(3);
	lepton.SetLineStyle(2);
	lepton.Draw("same F");
	
	//Combi Limit
	for(int i=1; i<=2; i++){
		TString temp="multicont_gr_"+std::to_string(i)+"_sm";
		std::cout<<temp<<std::endl;
		TGraph *st_exp = (TGraph*) file_5.Get(temp);
		if (combi.GetN()==0) combi.SetPoint(combi.GetN(),250,250);
		if (i==1) {
			for (int j=0; j<=st_exp->GetN()-1; j++) {
				combi.SetPoint(combi.GetN(),st_exp->GetX()[j],st_exp->GetY()[j]);
			}
		}
		else {
			for (int j=st_exp->GetN()-1; j>=0; j--) {
				combi.SetPoint(combi.GetN(),st_exp->GetX()[j],st_exp->GetY()[j]);
			}
			combi.SetPoint(combi.GetN(),1500,250);
		}
	}
	combi.SetFillColorAlpha(kBlack,0.1);
	combi.SetLineWidth(3);
	combi.SetLineStyle(2);
	combi.Draw("same F");
	
	//~ //Limit without correlated bkg uncertainties
	//~ TGraph *uncorr_exp = (TGraph*) file_6.Get("multicont_gr_1_sm");
	//~ uncorr_exp->SetLineColor(kBlack);
	//~ uncorr_exp->SetLineWidth(3);
	//~ uncorr_exp->SetLineStyle(2);
	//~ uncorr_exp->Draw("same");
	
	axis.Draw("axis same");
	TLatex text1(400,1000,"Lepton/Diphoton");
	TLatex text2(750,700,"All");
	TLatex text3(1000,600,"SinglePhoton/Lepton");
	TLatex text4(1250,750,"Lepton");
	TLatex text5(575,1300,"Combination");
	TLatex text6(450,1225,"Diphoton");
	text1.SetTextSize(0.03);
	text2.SetTextSize(0.03);
	text3.SetTextSize(0.03);
	text4.SetTextSize(0.03);
	text5.SetTextSize(0.03);
	text6.SetTextSize(0.03);
	text1.Draw("same");
	text2.Draw("same");
	text3.Draw("same");
	text4.Draw("same");
	text5.Draw("same");
	text6.Draw("same");
	
	
	saver.save(can,"Limits_GGM_combine_all",true,false);
	
	return 0;
}

int plot_T5Wg_allCombined(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_T5Wg_inclusiv.root","read");
	TFile file_2("../input/limits/limits_T5Wg_lepton.root","read");
	TFile file_3("../input/limits/limits_T5Wg_htg.root","read");
	TFile file_4("../input/limits/limits_T5Wg_htgHighLeptonVeto_leptonFull_htgHigh.root","read");
	TFile file_5("../input/limits/limits_T5Wg_diphoton.root","read");
	TFile file_6("../input/limits/limits_T5Wg_allExclusiv.root ","read");
	
	TGraph *st_exp = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *lepton_exp = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *htg_exp = (TGraph*) file_3.Get("gr_expC_sm");
	TGraph *combined_exp = (TGraph*) file_4.Get("gr_expC_sm");
	TGraph *combinedAll_exp = (TGraph*) file_6.Get("gr_expC_sm");
	TGraph *diphoton_exp = (TGraph*) file_5.Get("gr_expC_sm");
	TGraph *diag = new TGraph();
	
	diag->SetPoint(1,500,500);
	diag->SetPoint(2,2200,2200);
	
	TH2F axis("","",27,1400,2100,27,0,2200);
	axis.SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	axis.Draw("axis");
	axis.SetStats(0);
	
	st_exp->SetLineColor(kGreen);
	st_exp->SetLineWidth(3);
	st_exp->SetLineStyle(2);
	st_exp->Draw("same");
	
	lepton_exp->SetLineColor(kMagenta);
	lepton_exp->SetLineWidth(3);
	lepton_exp->SetLineStyle(2);
	lepton_exp->Draw("same");
	
	htg_exp->SetLineColor(kRed+1);
	htg_exp->SetLineWidth(3);
	htg_exp->SetLineStyle(2);
	htg_exp->Draw("same");
	
	combined_exp->SetLineColor(kBlue);
	combined_exp->SetLineWidth(3);
	combined_exp->SetLineStyle(2);
	combined_exp->Draw("same");
	
	combinedAll_exp->SetLineColor(kOrange);
	combinedAll_exp->SetLineWidth(3);
	combinedAll_exp->SetLineStyle(2);
	combinedAll_exp->Draw("same");
	
	diphoton_exp->SetLineColor(kCyan);
	diphoton_exp->SetLineWidth(3);
	diphoton_exp->SetLineStyle(2);
	diphoton_exp->Draw("same");
	
	diag->Draw("same");
	
	axis.Draw("axis same");
	
	gfx::LegendEntries legE;
	legE.append(*st_exp,"Photon+ST (exp)","l");
	legE.append(*lepton_exp,"Photon+Lepton (exp)","l");
	legE.append(*htg_exp,"Photon+HTG (exp)","l");
	legE.append(*diphoton_exp,"Diphoton (exp)","l");
	legE.append(*combined_exp,"Combined without Diphoton (exp)","l");
	legE.append(*combinedAll_exp,"Combined all (exp)","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_T5Wg_combineAll",true,false);
	
	return 0;
}

int plot_GGM2_diffcombi(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_GGM_M1_M3_lepton.root","read");
	TFile file_2("../input/limits/limits_GGM_M1_M3_htgHighVetoLeptonVeto_htgHighLeptonVeto_LeptonFull.root","read");
	TFile file_3("../input/limits/limits_GGM_M1_M3_leptonVeto_STLeptonVeto_leptonFull.root","read");
	
	
	TGraph *lepton = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *combined = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *combined2 = (TGraph*) file_3.Get("gr_expC_sm");
	
	TH2F axis("","",27,150,1500,27,1000,2500);
	axis.SetTitle(";M_{1} (GeV);M_{3} (GeV)");
	axis.Draw("axis");
	axis.SetStats(0);

	
	lepton->SetLineColor(kRed);
	lepton->SetLineWidth(3);
	lepton->SetLineStyle(2);
	lepton->Draw("same");
	
	combined->SetLineColor(kGreen);
	combined->SetLineWidth(3);
	combined->SetLineStyle(2);
	combined->Draw("same");
 
	combined2->SetLineColor(kBlue);
	combined2->SetLineWidth(3);
	combined2->SetLineStyle(2);
	combined2->Draw("same");
 
	gfx::LegendEntries legE;
	legE.append(*lepton,"Photon+Lepton (exp)","l");
	legE.append(*combined,"ST:highHTG veto HTG:highHTG bins (exp)","l");
	legE.append(*combined2,"ST:full HTG:ST veto (exp)","l");
	TLegend leg=legE.buildLegend(.2,.7,0.92,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_GGM2_diffCombis",true,false);
	
	return 0;
}

int plot_GGM2_combi(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_GGM_M1_M3_lepton.root","read");
	TFile file_2("../input/limits/limits_GGM_M1_M3_diphoton.root","read");
	TFile file_3("../input/limits/limits_GGM_M1_M3_allCombined_highHtgNN.root","read");
	
	
	TGraph *lepton = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *diphoton = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *combined = (TGraph*) file_3.Get("gr_expC_sm");
	
	TH2F axis("","",27,150,1500,27,1000,2500);
	axis.SetTitle(";M_{1} (GeV);M_{3} (GeV)");
	axis.Draw("axis");
	axis.SetStats(0);

	
	lepton->SetLineColor(kMagenta);
	lepton->SetLineWidth(3);
	lepton->SetLineStyle(2);
	lepton->Draw("same");
	
	diphoton->SetLineColor(kCyan);
	diphoton->SetLineWidth(3);
	diphoton->SetLineStyle(2);
	diphoton->Draw("same");
 
	combined->SetLineColor(kBlue);
	combined->SetLineWidth(3);
	combined->SetLineStyle(2);
	combined->Draw("same");
 
	gfx::LegendEntries legE;
	legE.append(*lepton,"Photon+Lepton (exp)","l");
	legE.append(*diphoton,"Diphoton (exp)","l");
	legE.append(*combined,"Combined (exp)","l");
	TLegend leg=legE.buildLegend(.7,.75,-1,-1,1);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"Limits_GGM2_allCombined",true,false);
	
	return 0;
}

int plot_GGM_combi_physmass(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	//~ TFile file_1("../input/limits/physmass_GGM_M1_M2_NN.root","read");
	TFile file_1("../input/limits/physmass_GGM_M1_M2_final.root","read");
	
	TH2F axis("","",27,120,720,27,120,1220);
	axis.SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
	axis.Draw("axis");
	axis.SetStats(0);
	
	TGraph *st = (TGraph*) file_1.Get("inclusivFinal/cont_sm");
	TGraph *lepton = (TGraph*) file_1.Get("lepton_final/cont_sm");
	TGraph *htg = (TGraph*) file_1.Get("htgFinal/cont_sm");
	TGraph *combined = (TGraph*) file_1.Get("allCombined_final/cont_sm");
	TGraph *diphoton = (TGraph*) file_1.Get("diphoton_final/cont_sm");
	TGraph *st_obs = (TGraph*) file_1.Get("inclusivFinal/cont_obs_sm");
	TGraph *lepton_obs = (TGraph*) file_1.Get("lepton_final/cont_obs_sm");
	TGraph *htg_obs = (TGraph*) file_1.Get("htgFinal/cont_obs_sm");
	TGraph *combined_obs = (TGraph*) file_1.Get("allCombined_final/cont_obs_sm");
	TGraph *diphoton_obs = (TGraph*) file_1.Get("diphoton_final/cont_obs_sm");
	TGraph *diag = new TGraph();
	
	diag->SetPoint(1,100,220);
	diag->SetPoint(2,2200,2320);
	
	lepton->SetLineColor(kMagenta);
	lepton->SetLineWidth(2);
	lepton->SetLineStyle(2);
	lepton->Draw("same");
	lepton_obs->SetLineColor(kMagenta);
	lepton_obs->SetLineWidth(2);
	lepton_obs->SetLineStyle(1);
	lepton_obs->Draw("same");
	
	st->SetLineColor(kGreen+1);
	st->SetLineWidth(2);
	st->SetLineStyle(2);
	st->Draw("same");
	st_obs->SetLineColor(kGreen+1);
	st_obs->SetLineWidth(2);
	st_obs->SetLineStyle(1);
	st_obs->Draw("same");
	
	htg->SetLineColor(kRed+1);
	htg->SetLineWidth(2);
	htg->SetLineStyle(2);
	htg->Draw("same");
	htg_obs->SetLineColor(kRed+1);
	htg_obs->SetLineWidth(2);
	htg_obs->SetLineStyle(1);
	htg_obs->Draw("same");
	
	diphoton->SetLineColor(kOrange-3);
	diphoton->SetLineWidth(2);
	diphoton->SetLineStyle(2);
	diphoton->Draw("same");
	diphoton_obs->SetLineColor(kOrange-3);
	diphoton_obs->SetLineWidth(2);
	diphoton_obs->SetLineStyle(1);
	diphoton_obs->Draw("same");
	
	combined->SetLineColor(kBlue);
	combined->SetLineWidth(2);
	combined->SetLineStyle(2);
	combined->Draw("same");
	combined_obs->SetLineColor(kBlue);
	combined_obs->SetLineWidth(2);
	combined_obs->SetLineStyle(1);
	combined_obs->Draw("same");
	
	diag->SetFillColor(kWhite);
	diag->SetLineColor(kBlack);
	
	diag->Draw("same F");
	diag->Draw("same");
	axis.Draw("axis same");
	
	TGraph *exp = new TGraph();
	TGraph *obs = new TGraph();
	exp->SetLineWidth(2);
	obs->SetLineWidth(2);
	exp->SetLineStyle(2);
 
	gfx::LegendEntries legE;
	legE.append(*st_obs,"Photon+S_{#scale[.8]{T}}^{#scale[.8]{#gamma}} ","l");
	legE.append(*lepton_obs,"Photon+Lepton ","l");
	legE.append(*htg_obs,"Photon+H_{#scale[.8]{T}}^{#scale[.8]{#gamma}} ","l");
	legE.append(*diphoton_obs,"Diphoton ","l");
	legE.append(*combined_obs,"Combination","l");
	legE.append(*diag,"m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} = m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} + 120 GeV","l");
	legE.prepend(*exp,"Expected","l");
	legE.prepend(*obs,"Observed","l");
	TLegend leg=legE.buildLegend(.2,.7,0.7,.9,2);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"LimitsPhysMass_GGM_combine",true,false);
	
	return 0;
}

int plot_GGM2_combi_physmass(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	//~ TFile file_1("../input/limits/physmass_GGM_M1_M3_NN.root","read");
	TFile file_1("../input/limits/physmass_GGM_M1_M3_final.root","read");
	
	TH2F axis("","",100,75,750,100,2350,5300);
	axis.SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV)");
	axis.Draw("axis");
	axis.SetStats(0);
	
	TGraph *lepton = (TGraph*) file_1.Get("lepton_final/cont_sm");
	TGraph *combined = (TGraph*) file_1.Get("allCombined_final/cont_sm");
	TGraph *diphoton = (TGraph*) file_1.Get("diphoton_final/cont_sm");
	TGraph *leptonObs = (TGraph*) file_1.Get("lepton_final/cont_obs_sm");
	TGraph *combinedObs = (TGraph*) file_1.Get("allCombined_final/cont_obs_sm");
	TGraph *diphotonObs = (TGraph*) file_1.Get("diphoton_final/cont_obs_sm");
	
	lepton->SetLineColor(kMagenta);
	lepton->SetLineWidth(2);
	lepton->SetLineStyle(2);
	lepton->Draw("same");
	leptonObs->SetLineColor(kMagenta);
	leptonObs->SetLineWidth(2);
	leptonObs->SetLineStyle(1);
	leptonObs->Draw("same");
	
	diphoton->SetLineColor(kOrange-3);
	diphoton->SetLineWidth(2);
	diphoton->SetLineStyle(2);
	diphoton->Draw("same");
	diphotonObs->SetLineColor(kOrange-3);
	diphotonObs->SetLineWidth(2);
	diphotonObs->SetLineStyle(1);
	diphotonObs->Draw("same");
	
	combined->SetLineColor(kBlue);
	combined->SetLineWidth(2);
	combined->SetLineStyle(2);
	combined->Draw("same");
	combinedObs->SetLineColor(kBlue);
	combinedObs->SetLineWidth(2);
	combinedObs->SetLineStyle(1);
	combinedObs->Draw("same");
	
	axis.Draw("axis same");
	
	TGraph *exp = new TGraph();
	TGraph *obs = new TGraph();
	exp->SetLineWidth(2);
	obs->SetLineWidth(2);
	exp->SetLineStyle(2);
 
	gfx::LegendEntries legE;
	legE.append(*leptonObs,"Photon+Lepton","l");
	legE.append(*diphotonObs,"Diphoton","l");
	legE.append(*combinedObs,"Combination","l");
	legE.prepend(*exp,"Expected","l");
	legE.prepend(*obs,"Observed","l");
	TLegend leg=legE.buildLegend(.3,.75,0.7,.9,2);
	leg.SetTextSize(0.03);
	leg.Draw();
	
	saver.save(can,"LimitsPhysMass_GGM2_combine",true,false);
	
	return 0;
}

int plot_TChiNg_diffcombi(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	can.SetLogy();
	
	TFile file_1("../output/aux_limits_TChiNg.root","read");
	
	TGraph *xsec = (TGraph*) file_1.Get("inclusiv/xs");
	//~ TGraph *st_exp = (TGraph*) file_1.Get("inclusiv/exp");
	TGraph *combi_exp = (TGraph*) file_1.Get("inclusiv_stVeto/exp");
	TGraph *combi2_exp = (TGraph*) file_1.Get("htgHighVeto_htgHigh/exp");
	TGraph *combi3_exp = (TGraph*) file_1.Get("htgVeto_fullHtg/exp");
	TH2F axis("","",100,950,1100,100,0.00095,0.0035);
	//~ axis.GetYaxis()->SetNdivisions(0);
	axis.GetYaxis()->SetMoreLogLabels();
	axis.GetYaxis()->SetNoExponent(0);
	axis.SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV); 95% CL upper limit / cross section (pb)");
	axis.Draw("axis");
	axis.SetStats(0);
	
	xsec->SetLineColor(kBlack);
	xsec->SetLineWidth(3);
	xsec->Draw("same");
	
	//~ st_exp->SetLineColor(kGreen);
	//~ st_exp->SetLineWidth(3);
	//~ st_exp->SetLineStyle(2);
	//~ st_exp->Draw("same");
	
	combi_exp->SetLineColor(kGreen+1);
	combi_exp->SetLineWidth(3);
	combi_exp->SetLineStyle(2);
	combi_exp->Draw("same");
	
	combi2_exp->SetLineColor(kBlue);
	combi2_exp->SetLineWidth(3);
	combi2_exp->SetLineStyle(2);
	combi2_exp->Draw("same");
	
	combi3_exp->SetLineColor(kRed);
	combi3_exp->SetLineWidth(3);
	combi3_exp->SetLineStyle(2);
	combi3_exp->Draw("same");
	
	axis.Draw("axis same");
	
	gfx::LegendEntries legE;
	legE.append(*xsec,"Theory cross section","l");
	legE.append(*combi2_exp,"Strategy 3 (exp)","l");
	legE.append(*combi3_exp,"Strategy 2 (exp)","l");
	legE.append(*combi_exp,"Strategy 1 (exp)","l");
	//~ legE.append(*st_exp,"Photon+ST (exp)","l");
	TLegend leg=legE.buildLegend(.7,.3,0.9,.5,1);
	leg.SetTextSize(0.033);
	leg.Draw();
	
	saver.save(can,"Limits_TChiNg_diffcombi",true,false);
	
	return 0;
}

int plot_GGM2_combi_filledAreas(){
	io::RootFileSaver saver("plots.root",TString::Format("danilo_plot_combined%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_GGM_M1_M3_lepton_final.root","read");
	TFile file_2("../input/limits/limits_GGM_M1_M3_diphoton_final.root","read");
	TFile file_3("../input/limits/limits_GGM_M1_M3_allCombined_final.root","read");
	
	
	TGraph *lepton = (TGraph*) file_1.Get("gr_expC_sm");
	TGraph *diphoton = (TGraph*) file_2.Get("gr_expC_sm");
	TGraph *combined = (TGraph*) file_3.Get("gr_expC_sm");
	
	TH2F axis("","",27,150,1500,27,1000,2500);
	axis.SetTitle(";#it{M}_{1} (GeV);#it{M}_{3} (GeV)");
	axis.Draw("axis");
	axis.SetStats(0);

	lepton->SetPoint(lepton->GetN(),0,0);
	lepton->SetFillColorAlpha(kOrange,0.2);
	lepton->SetLineWidth(3);
	lepton->SetLineStyle(2);
	lepton->Draw("same F");
	
	diphoton->SetPoint(diphoton->GetN(),0,0);
	diphoton->SetFillColorAlpha(kCyan,0.3);
	diphoton->SetLineWidth(3);
	diphoton->SetLineStyle(2);
	diphoton->Draw("same F");
	
	combined->SetPoint(combined->GetN(),1500,1000);
	combined->SetPoint(combined->GetN(),0,0);
	combined->SetFillColorAlpha(kBlack,0.1);
	combined->SetLineWidth(3);
	combined->SetLineStyle(2);
	combined->Draw("same F");
	
	TLatex text1(1230,1010,"Lepton");
	TLatex text2(960,1320,"Combination");
	TLatex text3(170,2020,"Diphoton");
	TLatex text4(400,1350,"Lepton/Diphoton");
	text1.SetTextSize(0.03);
	text2.SetTextSize(0.03);
	text3.SetTextSize(0.03);
	text4.SetTextSize(0.03);
	text1.Draw("same");
	text2.Draw("same");
	text3.Draw("same");
	text4.Draw("same");
	
	saver.save(can,"Limits_GGM2_allCombined_filledAreas",true,false);
	
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
	//~ plot_T6gg_htgHighVeto_htgHigh();
	//~ plot_T6Wg_htgHighVeto_htgHigh();
	//~ plot_T5gg_htgHighVeto_htgHigh();
	//~ plot_T5Wg_htgHighVeto_htgHigh();
	//~ plot_GGM_htgHighVeto_htgHigh();
	//~ plot_T6gg_diffCombi();
	//~ plot_TChiNg_BR_diffcombi();
	//~ plot_CharginoBRstrongN1700();
	//~ plot_CharginoBRstrongG1950();
	//~ plot_CharginoBRstrongG1800();
	//~ plot_CharginoBRstrongG1700();
	//~ plot_T5Wg_allCombined();
	//~ plot_GGM2_diffcombi();
	//~ plot_GGM_combi();
	//~ plot_GGM2_combi();
	//~ plot_GGM2_combi_filledAreas();
	plot_GGM_combi_physmass();
	plot_GGM2_combi_physmass();
	//~ plot_TChiNg_diffcombi();
	plot_CharginoBR_C1C1();
	plot_T5Wg_htgHighLeptonVeto_leptonFull_htgHigh();
	//~ plot_T5Wg_diffCombi();
	plot_TChiNg_BR_combi();
	//~ plot_T5Wg_diphotonCase();
}
