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


extern "C"
void run()
{
      for (TString scan :{"M1_M2","M1_M3","T5Wg"}){
            io::RootFileSaver saver("plots.root","danilo_plotSig");
            TCanvas can;
            
            TString fileLoc="../input/limits/Significance/significance_GGM_"+scan+"_allCombined_final.root";
            if (scan=="T5Wg"){
				fileLoc="../input/limits/Significance/significance_"+scan+"_ST_HTG_lepton_final.root";
			}

            TFile file3(fileLoc,"read");

			TH2F *hist = (TH2F*) file3.Get("h_sig");
			can.cd();

			gPad->SetRightMargin(0.2);
			gPad->SetLeftMargin(0.13);
			gPad->SetBottomMargin(0.10);
			hist->GetXaxis()->SetTitle("#it{M}_{1} (GeV)");
			hist->GetYaxis()->SetTitle("#it{M}_{2} (GeV)");
			if (scan=="M1_M3") hist->GetYaxis()->SetTitle("#it{M}_{3} (GeV)");
			else if (scan=="T5Wg"){
				hist->SetTitle(";m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV);m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
			}

			hist->GetZaxis()->SetTitle("Observed significance");

			hist->GetYaxis()->SetTitleOffset(1.3);
			hist->GetXaxis()->SetTitleOffset(0.9);
			hist->GetZaxis()->SetTitleOffset(1.3);
			hist->GetYaxis()->SetTitleSize(0.05);
			hist->GetXaxis()->SetTitleSize(0.05);
			hist->GetZaxis()->SetTitleSize(0.05);
			hist->GetYaxis()->SetLabelSize(0.04);
			hist->GetXaxis()->SetLabelSize(0.04);
			hist->GetZaxis()->SetLabelSize(0.04);
			hist->SetAxisRange(250,1500,"X");
			hist->SetAxisRange(250,1500,"Y");
			if (scan=="M1_M3") {
				hist->SetAxisRange(50,1500,"X");
				hist->SetAxisRange(1000,2500,"Y");
			}
			else if (scan=="T5Wg") {
				hist->SetAxisRange(1400,2100,"X");
				hist->SetAxisRange(0,2400,"Y");
			}
			hist->SetStats(false);
			hist->Draw("colz");
			can.RedrawAxis();
			saver.save(can,scan,true,false);
			can.Clear();
      }
}
