//Transforming the Tgraph2d into a TH2D for the acceptance
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
#include <TGraph2D.h>
#include <TMath.h>


extern "C"
void run()
{     
      for (TString scan :{"GGM_M1_M2","GGM_M1_M3","T5Wg"}) {
            
            std::vector<TString> sSelections = {"exclusiv","inclusiv","htgVeto","leptonVeto","diphotonVeto"};
            if (scan == "T5Wg") {
                  sSelections.push_back("newTrig");
                  sSelections.push_back("newTrig200");
                  sSelections.push_back("htgHighVeto");
            }
            
            for (TString selection :sSelections) {
                  
                  TFile file("../output/signal_scan_"+selection+"_v03D.root","read");
                  
                  std::vector<TString> sVar = {"acceptance"};
                  if (scan == "GGM_M1_M2" && selection=="inclusiv") {
                        sVar.push_back("nEvents");
                        sVar.push_back("nEventsBin3");
                        sVar.push_back("nEventsBin4");
                        sVar.push_back("sqrt(N)/N");
                        sVar.push_back("sqrt(N)/N Bin3");
                        sVar.push_back("sqrt(N)/N Bin4");
                  }
                  
                  for (TString var :sVar) {
                  
                        io::RootFileSaver saver("plots.root","danilo_acceptanceHist");
                        TCanvas can;
                        TString path;
                        
                        if (var == "sqrt(N)/N") {
                              path = scan+"/pre_ph165/c_MET300/MT300/STg/"+scan+"_nEvents";
                        }
                        else if (var == "sqrt(N)/N Bin3") {
                              path = scan+"/pre_ph165/c_MET300/MT300/STg/"+scan+"_nEventsBin3";
                        }
                        else if (var == "sqrt(N)/N Bin4") {
                              path = scan+"/pre_ph165/c_MET300/MT300/STg/"+scan+"_nEventsBin4";
                        }
                        else {
                              path = scan+"/pre_ph165/c_MET300/MT300/STg/"+scan+"_"+var;
                        }
                        TGraph2D *graph = (TGraph2D*) file.Get(path);
                        can.cd();
                        can.SetLogz();
                        TH2D *hist = graph->GetHistogram();
                        
                        if (var == "sqrt(N)/N" || var == "sqrt(N)/N Bin3" || var == "sqrt(N)/N Bin4") {
                              for (int i=1; i<hist->GetNbinsX()+1; i++) {
                                    for (int j=1; j<hist->GetNbinsX()+1; j++) {
                                          if (hist->GetBinContent(i,j)==0) hist->SetBinContent(i,j,0.);
                                          else hist->SetBinContent(i,j,1.0/TMath::Sqrt(hist->GetBinContent(i,j)));
                                    }
                              }
                        }
                        gPad->SetRightMargin(0.2);
                        gPad->SetLeftMargin(0.13);
                        gPad->SetBottomMargin(0.10);
                        hist->SetBit(TH1::kNoTitle);
                        hist->GetXaxis()->SetTitle("M_{1} (GeV)");
                        if (scan=="GGM_M1_M2") hist->GetYaxis()->SetTitle("M_{2} (GeV)");
                        else if (scan=="GGM_M1_M3") hist->GetYaxis()->SetTitle("M_{3} (GeV)");
                        else if (scan=="T5Wg") {
                              hist->GetXaxis()->SetTitle("m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV)");
                              hist->GetYaxis()->SetTitle("m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
                        }
                        if (var == "acceptance") {
                              hist->GetZaxis()->SetTitle("acceptance x eff");
                        }
                        else {
                              hist->GetZaxis()->SetTitle(var);
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
                        hist->GetZaxis()->SetLabelOffset(0.02);
                        hist->SetStats(false);
                        if (var=="acceptance") hist->SetMaximum(0.4);
                        hist->Draw("colz");
                  
                        can.RedrawAxis();
                        saver.save(can,scan+"_"+selection+"_"+var,true,true);
                        can.Clear();
                  }
            }
      }
      
      //Save acceptance for TChiNg and possible new trigger
      TFile file1("../output/signal_scan_newTrig200_v03D.root","read");
      TFile file2("../output/signal_scan_newTrig_v03D.root","read");
                  
      io::RootFileSaver saver("plots.root","danilo_acceptanceHist");
      TCanvas can;
      
      TString path = "TChiNg/pre_ph165/c_MET300/MT300/STg/TChiNg_acceptance";
      TGraph *graph1 = (TGraph*) file1.Get(path);
      TGraph *graph2 = (TGraph*) file2.Get(path);
      can.cd();
      graph1->GetXaxis()->SetTitle("m_{NLSP} (GeV)");
      graph1->GetYaxis()->SetTitle("acceptance");
      graph2->SetMarkerColor(kRed);
      //~ graph1->SetStats(false);
      graph1->Draw("AP");
      graph2->Draw("P same");
      
      gfx::LegendEntries legE;
	legE.append(*graph1,"Loose Photon p_{T}>200 GeV","p");
	legE.append(*graph2,"Tight Photon p_{T}>100 GeV","p");
      TLegend leg=legE.buildLegend(.2,.7,.7,.9,1);
	leg.SetTextSize(0.03);
	leg.Draw();

      can.RedrawAxis();
      saver.save(can,"TChiNg_newTrig_acceptance",true,true);
      can.Clear();
      
      
}
