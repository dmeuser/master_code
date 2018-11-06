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

void M1_M2(){
      io::RootFileSaver saver("plots.root","danilo_GGM_scan100.0");
      TCanvas can;
   

      TFile file1("../output/GGM_kinematics_100.0.root","read");
   
      for (TString svar :{"PhotonPt","MET","ST", "MT", "HTG"}) {
            
            TH2F *hist = (TH2F*) file1.Get("M1_M2/"+svar);
            can.cd();
            gPad->SetRightMargin(0.2);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.10);
            hist->GetXaxis()->SetTitle("M_{1} (GeV)");
            hist->GetYaxis()->SetTitle("M_{2} (GeV)");
            
            if (svar=="PhotonPt") {
                  hist->GetZaxis()->SetTitle("mean p_{#scale[.8]{T}}^{#gamma_{1}} (GeV)");
            }
            else if (svar=="MET") {
                  hist->GetZaxis()->SetTitle("mean p_{#scale[.8]{T}}^{#scale[.8]{miss}} (GeV)");
            }
            else if (svar=="ST") {
                  hist->GetZaxis()->SetTitle("mean S_{#scale[.8]{T}}^{#scale[.8]{#gamma}} (GeV)");
            }
            else if (svar=="MT") {
                  hist->GetZaxis()->SetTitle("mean MT(#gamma_{1},p_{T}^{miss}) (GeV)");
            }
            else if (svar=="HTG") {
                  hist->GetZaxis()->SetTitle("mean H_{#scale[.8]{T}}^{#scale[.8]{#gamma}} (GeV)");
            }
            
            hist->GetYaxis()->SetTitleOffset(1.3);
            hist->GetXaxis()->SetTitleOffset(0.9);
            hist->GetZaxis()->SetTitleOffset(1.2);
            hist->GetYaxis()->SetTitleSize(0.05);
            hist->GetXaxis()->SetTitleSize(0.05);
            hist->GetZaxis()->SetTitleSize(0.05);
            hist->GetYaxis()->SetLabelSize(0.04);
            hist->GetXaxis()->SetLabelSize(0.04);
            hist->GetZaxis()->SetLabelSize(0.04);
            hist->GetZaxis()->SetLabelOffset(0.01);
            hist->SetAxisRange(250,1500,"X");
            hist->SetAxisRange(250,1500,"Y");
            hist->SetStats(false);
            hist->Draw("colz");
            can.RedrawAxis();
            saver.save(can,"M1_M2/"+svar,true,true);
            can.Clear();
      }
}

void M1_M3(){
     io::RootFileSaver saver("plots.root","danilo_GGM_scan100.0");
      TCanvas can;
   

      TFile file1("../output/GGM_kinematics_100.0.root","read");
   
      for (TString svar :{"PhotonPt","MET","ST", "MT", "HTG"}) {
            
            TH2F *hist = (TH2F*) file1.Get("M1_M3/"+svar);
            can.cd();
            gPad->SetRightMargin(0.2);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.10);
            hist->GetXaxis()->SetTitle("M_{1} (GeV)");
            hist->GetYaxis()->SetTitle("M_{3} (GeV)");
            
            if (svar=="PhotonPt") {
                  hist->GetZaxis()->SetTitle("mean p_{#scale[.8]{T}}^{#gamma_{1}} (GeV)");
            }
            else if (svar=="MET") {
                  hist->GetZaxis()->SetTitle("mean p_{#scale[.8]{T}}^{#scale[.8]{miss}} (GeV)");
            }
            else if (svar=="ST") {
                  hist->GetZaxis()->SetTitle("mean S_{#scale[.8]{T}}^{#scale[.8]{#gamma}} (GeV)");
            }
            else if (svar=="MT") {
                  hist->GetZaxis()->SetTitle("mean MT(#gamma_{1},p_{T}^{miss}) (GeV)");
            }
            else if (svar=="HTG") {
                  hist->GetZaxis()->SetTitle("mean H_{#scale[.8]{T}}^{#scale[.8]{#gamma}} (GeV)");
            }
            
            hist->GetYaxis()->SetTitleOffset(1.3);
            hist->GetXaxis()->SetTitleOffset(0.9);
            hist->GetZaxis()->SetTitleOffset(1.2);
            hist->GetYaxis()->SetTitleSize(0.05);
            hist->GetXaxis()->SetTitleSize(0.05);
            hist->GetZaxis()->SetTitleSize(0.05);
            hist->GetYaxis()->SetLabelSize(0.04);
            hist->GetXaxis()->SetLabelSize(0.04);
            hist->GetZaxis()->SetLabelSize(0.04);
            hist->GetZaxis()->SetLabelOffset(0.01);
            hist->SetAxisRange(150,1500,"X");
            hist->SetAxisRange(1000,2500,"Y");
            hist->SetStats(false);
            hist->Draw("colz");
            can.RedrawAxis();
            saver.save(can,"M1_M3/"+svar,true,true);
            can.Clear();
      }
}

void T5Wg(){
      io::RootFileSaver saver("plots.root","danilo_GGM_scan100.0");
      TCanvas can;
   

      TFile file1("../output/GGM_kinematics_100.0.root","read");
   
      for (TString svar :{"PhotonPt","MET","ST", "MT", "HTG"}) {
            
            TH2F *hist = (TH2F*) file1.Get("T5Wg/"+svar);
            can.cd();
            gPad->SetRightMargin(0.2);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.10);
            hist->GetXaxis()->SetTitle("m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV)");
            hist->GetYaxis()->SetTitle("m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
            
            if (svar=="PhotonPt") {
                  hist->GetZaxis()->SetTitle("mean p_{#scale[.8]{T}}^{#gamma_{1}} (GeV)");
            }
            else if (svar=="MET") {
                  hist->GetZaxis()->SetTitle("mean p_{#scale[.8]{T}}^{#scale[.8]{miss}} (GeV)");
            }
            else if (svar=="ST") {
                  hist->GetZaxis()->SetTitle("mean S_{#scale[.8]{T}}^{#scale[.8]{#gamma}} (GeV)");
            }
            else if (svar=="MT") {
                  hist->GetZaxis()->SetTitle("mean MT(#gamma_{1},p_{T}^{miss}) (GeV)");
            }
            else if (svar=="HTG") {
                  hist->GetZaxis()->SetTitle("mean H_{#scale[.8]{T}}^{#scale[.8]{#gamma}} (GeV)");
            }
            
            hist->GetYaxis()->SetTitleOffset(1.3);
            hist->GetXaxis()->SetTitleOffset(0.9);
            hist->GetZaxis()->SetTitleOffset(1.2);
            hist->GetYaxis()->SetTitleSize(0.05);
            hist->GetXaxis()->SetTitleSize(0.05);
            hist->GetZaxis()->SetTitleSize(0.05);
            hist->GetYaxis()->SetLabelSize(0.04);
            hist->GetXaxis()->SetLabelSize(0.04);
            hist->GetZaxis()->SetLabelSize(0.04);
            hist->GetZaxis()->SetLabelOffset(0.01);
            hist->SetAxisRange(1400,2100,"X");
            hist->SetAxisRange(0,2200,"Y");
            hist->SetStats(false);
            hist->Draw("colz");
            can.RedrawAxis();
            saver.save(can,"T5Wg/"+svar,true,true);
            can.Clear();
      }
}


extern "C"

void run()
{
   M1_M2();
   M1_M3();
   T5Wg();
}
