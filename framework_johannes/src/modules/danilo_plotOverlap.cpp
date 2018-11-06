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

int plotOverlap_ST() {
      io::RootFileSaver saver("plots.root","danilo_plotOverlap");
      TCanvas can;

      TString fileLoc="../output/overlap.root";
      TFile file(fileLoc,"read");
      TH2F *input = (TH2F*) file.Get("Overlap_ST100.0");
      TH2F hist("Overlap_2d","",4,0.5,4.5,5,0.5,5.5);
      can.cd();
      
      int i=1;
      for(TString bin: {"600 - 800","800 - 1000","1000 - 1300",">1300"}){
            hist.GetXaxis()->SetBinLabel(i,bin);
            i++;
      }
      
      std::vector<int> count;
      count = {281,101,65,24};
      i=1;
      for(TString bin: {"Inclusive","Exclusive","Photon+H_{#scale[.8]{T}}^{#scale[.8]{#gamma}}","Lepton","Diphoton"}){
            for (int j = 1; j <= 4; j++) {
                  if (bin!="Inclusive") hist.SetBinContent(j,i,input->GetBinContent(j,i-1));
                  else hist.SetBinContent(j,i,count[j-1]);
            }
            hist.GetYaxis()->SetBinLabel(i,bin);
            i++;
      }
      
      
      hist.SetMarkerColor(kRed);
      hist.GetYaxis()->SetTitleOffset(2);
      hist.GetXaxis()->SetTickLength(0);
      hist.GetYaxis()->SetTickLength(0);
      hist.SetTitle(";STg;Analysis");
      hist.SetStats(false);
      hist.SetMarkerSize(2.3);
      hist.Draw("col TEXT");
      can.SetLogz(true);
      saver.save(can,"Overlap_2d_ST",true,false);
      can.Clear();
      
      return 0;
}

int plotOverlap_HTG() {
      io::RootFileSaver saver("plots.root","danilo_plotOverlap");
      TCanvas can;

      TString fileLoc="../input/overlap_HTG.root";
      TFile file(fileLoc,"read");
      TH2F hist("Overlap_2d","",6,0.5,6.5,5,0.5,5.5);
      can.cd();
      
      int i = 2;
      for (TString analys: {"HTG","ST","Lepton","Diphoton"}) {
            TH1F *histTemp = (TH1F*) file.Get(analys);
            for (int j = 1; j <= 6; j++) {
                  hist.SetBinContent(j,i,histTemp->GetBinContent(j));
            }
            hist.GetYaxis()->SetBinLabel(i,analys);
            i++;
      }
      
      i=1;
      for (int count: {103,82,21,6,10,4}) {
            hist.SetBinContent(i,1,count);
            i++;
      }
      
      i=1;
      for(TString bin: {"350 - 450","450 - 600",">600","350 - 450","450 - 600",">600"}){
            hist.GetXaxis()->SetBinLabel(i,bin);
            i++;
      }
      
      i=1;
      for(TString bin: {"Inclusive","Exclusive","Photon+S_{#scale[.8]{T}}^{#scale[.8]{#gamma}}","Lepton","Diphoton"}){
            hist.GetYaxis()->SetBinLabel(i,bin);
            i++;
      }
      
      hist.SetMarkerColor(kRed);
      hist.GetYaxis()->SetTitleOffset(2);
      hist.GetXaxis()->SetTickLength(0);
      hist.GetYaxis()->SetTickLength(0);
      hist.SetTitle(";%MET;Analysis");
      hist.SetStats(false);
      hist.SetMarkerSize(2.3);
      hist.Draw("col TEXT");
      
      TLine line(3.5,0.5,3.5,5.5);
      line.SetLineWidth(3);
      line.SetLineStyle(2);
      line.Draw("same");
      
      TLatex text;
      text.SetTextSize(0.7*text.GetTextSize());
      text.DrawLatex(1.5,0.6,"H_{#scale[.8]{T}}^{#scale[.8]{#gamma}} < 2 TeV");
      text.DrawLatex(4.5,0.6,"H_{#scale[.8]{T}}^{#scale[.8]{#gamma}} > 2 TeV");
      
      can.SetLogz(true);
      saver.save(can,"Overlap_2d_HTG",true,false);
      can.Clear();
      
      return 0;
}

extern "C"
void run(){
      plotOverlap_ST();
      plotOverlap_HTG();
}
