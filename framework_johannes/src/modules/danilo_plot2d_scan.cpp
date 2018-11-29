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
io::RootFileSaver saver("plots.root","danilo_plot2d_scan");
TCanvas can;

//~ for (TString order :{"NLO","LO"}){
for (TString order :{"NLO"}){
         //With LHAPDF modified prospino
         TString fileLoc="../input/xsecs_M1_M2.root";
         if (order=="LO") fileLoc="../input/xsecs_M1_M2_LO.root";

         TFile file3(fileLoc,"read");
         
         for (TString svar :{"xsecs","pdf_uncertainty","alphaS_uncertainty", "scaleUp_uncertainty", "scaleDown_uncertainty","total_uncertainty","alphaUp","alphaDown","alphaDiff","N1C1_MassPlane_hist","hist_KN2C1","hist_KC1C1","hist_KN1C1"}) {
            TH2F *hist = (TH2F*) file3.Get(svar);
            can.cd();
            if (svar == "xsecs" || svar == "N1C1_MassPlane_hist") {
               can.SetLogz();
               hist->SetTitle(";;");
            }
            else can.SetLogz(false);
            gPad->SetRightMargin(0.2);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.10);
            hist->GetXaxis()->SetTitle("#it{M}_{1} (GeV)");
            hist->GetYaxis()->SetTitle("#it{M}_{2} (GeV)");
            
            if (svar == "N1C1_MassPlane_hist") {
                  hist->GetXaxis()->SetTitle("M_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV)");
                  hist->GetYaxis()->SetTitle("M_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV)");
            }
            
            if (svar == "xsecs" || svar == "N1C1_MassPlane_hist") {
                  hist->GetZaxis()->SetTitle("Cross section (pb)");
                  hist->SetMinimum(0.0001);
            }
            else if (svar=="pdf_uncertainty") hist->GetZaxis()->SetTitle("PDF uncertainty (%)");
            else if (svar=="alphaS_uncertainty") hist->GetZaxis()->SetTitle("#alpha_{s} uncertainty (%)");
            else if (svar=="scaleUp_uncertainty") hist->GetZaxis()->SetTitle("Scale-up uncertainty (%)");
            else if (svar=="scaleDown_uncertainty") hist->GetZaxis()->SetTitle("Scale-down uncertainty (%)");
            else if (svar=="total_uncertainty") hist->GetZaxis()->SetTitle("Total uncertainty (%)");
            
            if (svar == "scale_uncertainty") {
                  hist->SetMaximum(1.5);
            }
            
            if (svar == "alphaS_uncertainty") {
                  for (int i = 1; i<=hist->GetNbinsX(); i++){
                        for (int j = 1; j<=hist->GetNbinsY(); j++){
                              if (hist->GetBinContent(i,j)==0) hist->SetBinContent(i,j,0.00000001);
                        }
                  }
            }
            
            if (svar == "hist_KN2C1" || svar == "hist_KN1C1" || svar == "hist_KC1C1") {
				hist->SetMinimum(1.04);
				hist->SetMaximum(1.40);
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
            hist->SetStats(false);
            hist->Draw("colz");
            can.RedrawAxis();
            TString plotLoc="M1_M2_modProspino/"+svar;
            if (order=="LO") plotLoc="M1_M2_modProspino_LO/"+svar;
            saver.save(can,plotLoc,true,true);
            can.Clear();
            
         }
         
         //With LHAPDF modified prospino second scan
         fileLoc="../input/xsecs_M1_M3.root";
         if (order=="LO") fileLoc="../input/xsecs_M1_M3_LO.root";

         TFile file4(fileLoc,"read");
         
         for (TString svar :{"xsecs","pdf_uncertainty","alphaS_uncertainty", "scaleUp_uncertainty", "scaleDown_uncertainty","total_uncertainty","alphaUp","alphaDown","alphaDiff","N1C1_MassPlane_hist","hist_KN2C1","hist_KC1C1"}) {
            
            TH2F *hist = (TH2F*) file4.Get(svar);
            can.cd();
            if (svar == "xsecs" || svar == "N1C1_MassPlane_hist") {
               can.SetLogz();
               hist->SetTitle(";;");
            }
            else can.SetLogz(false);
            gPad->SetRightMargin(0.2);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.10);
            hist->GetXaxis()->SetTitle("#it{M}_{1} (GeV)");
            hist->GetYaxis()->SetTitle("#it{M}_{3} (GeV)");
            
            if (svar == "N1C1_MassPlane_hist") {
                  hist->GetXaxis()->SetTitle("M_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV)");
                  hist->GetYaxis()->SetTitle("M_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV)");
            }
            
            if (svar == "xsecs" || svar == "N1C1_MassPlane_hist") {
                  hist->GetZaxis()->SetTitle("Cross section (pb)");
                  hist->SetMinimum(0.00001);
            }
            else if (svar=="pdf_uncertainty") hist->GetZaxis()->SetTitle("PDF uncertainty (%)");
            else if (svar=="alphaS_uncertainty") hist->GetZaxis()->SetTitle("#alpha_{s} uncertainty (%)");
            else if (svar=="scaleUp_uncertainty") hist->GetZaxis()->SetTitle("Scale-up uncertainty (%)");
            else if (svar=="scaleDown_uncertainty") hist->GetZaxis()->SetTitle("Scale-down uncertainty (%)");
            else if (svar=="total_uncertainty") hist->GetZaxis()->SetTitle("Total uncertainty (%)");
            
            if (svar == "scale_uncertainty") {
                  hist->SetMaximum(1.5);
            }
            
            if (svar == "hist_KN2C1" || svar == "hist_KC1C1") {
				hist->SetMinimum(0.95);
				hist->SetMaximum(1.30);
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
            hist->SetAxisRange(50,1500,"X");
            hist->SetAxisRange(1000,2500,"Y");
            hist->SetStats(false);
            hist->Draw("colz");
            can.RedrawAxis();
            TString plotLoc="M1_M3_modProspino/"+svar;
            if (order=="LO") plotLoc="M1_M3_modProspino_LO/"+svar;
            saver.save(can,plotLoc,true,true);
            can.Clear();
            
         }
   }
}
