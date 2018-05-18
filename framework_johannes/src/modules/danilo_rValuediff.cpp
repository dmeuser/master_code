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
	
	io::RootFileSaver saver("plots.root",TString::Format("danilo_rValuediff%.1f/%s",cfg.processFraction*100,"Limits"));
	TCanvas can;
	
	TFile file_1("../input/limits/limits_T6gg_inclusiv.root","read");
	TFile file_2("../input/limits/limits_T6gg_htg.root","read");
	
	TH2F *johannes_exp = (TH2F*) file_1.Get("h_exp");
	TH2F *knut_exp = (TH2F*) file_2.Get("h_exp");
	
	johannes_exp->Add(knut_exp,-1);
	johannes_exp->SetStats(false);
	
	johannes_exp->SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV);r value diff");
	
	johannes_exp->Draw("colz");
	
	saver.save(can,"T6gg_ST_HTG",true,false);
	
	return 0;
}
extern "C"

void run(){
	plot();
}
