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
      for (TString scan :{"M1_M2","M1_M3"}){
            io::RootFileSaver saver("plots.root","danilo_plotGGMparameters");
            TCanvas can;
            
            TString fileLoc="../input/GGM_scan/hist2d.root";
            if (scan=="M1_M3") fileLoc="../input/GGM_scan/hist2d_scan2.root";

            TFile file3(fileLoc,"read");

            for (TString svar :{"massNLSP","massChargino","massGluino", "massNeutralino2","massDiff","BRtoPhoton","BRtoZ"}) {
                  
                  if (scan=="M1_M3" && svar=="massDiff") continue;
                  
                  TH2F *hist = (TH2F*) file3.Get(svar);
                  can.cd();
                  
                  gPad->SetRightMargin(0.2);
                  gPad->SetLeftMargin(0.13);
                  gPad->SetBottomMargin(0.10);
                  hist->GetXaxis()->SetTitle("#it{M}_{1} (GeV)");
                  hist->GetYaxis()->SetTitle("#it{M}_{2} (GeV)");
                  if (scan=="M1_M3") hist->GetYaxis()->SetTitle("#it{M}_{3} (GeV)");
                  
                  if (svar=="massChargino") hist->GetZaxis()->SetTitle("m_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV)");
                  else if (svar=="massNLSP") hist->GetZaxis()->SetTitle("m_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV)");
                  else if (svar=="massNeutralino2") hist->GetZaxis()->SetTitle("m_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{2}}}} (GeV)");
                  else if (svar=="massGluino") hist->GetZaxis()->SetTitle("m_{#lower[-0.12]{#tilde{g}}} (GeV)");
                  else if (svar=="massDiff") hist->GetZaxis()->SetTitle("|m_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}-m_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}| (GeV)");
                  else if (svar=="BRtoPhoton") hist->GetZaxis()->SetTitle("BF(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #gamma + #tilde{G}) (%)");
                  else if (svar=="BRtoZ") hist->GetZaxis()->SetTitle("BF(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow Z + #tilde{G}) (%)");
                  
                  if (svar=="BRtoZ" or svar=="BRtoPhoton") {
					  for(int i=1; i<=hist->GetNbinsX(); i++){
						  for(int j=1; j<=hist->GetNbinsY(); j++){
							  if(hist->GetBinContent(i,j)==0) hist->SetBinContent(i,j,0.0001);
						  }
					  }
				  }
                  
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
                  hist->SetStats(false);
                  hist->Draw("colz");
                  can.RedrawAxis();
                  TString plotLoc=scan+"/"+svar;
                  saver.save(can,plotLoc,true,true);
                  can.Clear();
            }
      }
}
