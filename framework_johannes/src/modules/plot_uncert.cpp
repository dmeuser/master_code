#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/io.hpp"
#include "tools/util.hpp"

#include <TGraph2D.h>
#include <TGraphSmooth.h>

Config const &cfg=Config::get();

enum Scan_t
{
   TChiWg,
   TChiNg,
   T5gg,
   T5Wg,
   T6gg,
   T6Wg,
   GGM,
};
std::map<int,float> getXsecs(Scan_t scan)
{
   std::string path=CMAKE_SOURCE_DIR;
   if (scan==TChiWg) path += "/xsec_N2C1_wino.csv";
   if (scan==TChiNg) path += "/xsec_comb_wino.csv";  
   unsigned const nCol=3;
   std::ifstream fstream (path,std::ifstream::in);
   std::string line;
   std::map<int,float> m_M_XS;
   while (fstream.good()) {
      std::getline(fstream, line);
      if (line.find("#") == 0) continue;
      line=util::rm_duplicate_spaces(line);
      std::vector<float> values=util::to_vector<float>(line,' ');
      if (values.size() == nCol) {
         if (scan==TChiWg || scan==TChiNg) {
            m_M_XS[(int)values[0]]=values[1]/1000.; // convert to pb
         } else {
            debug<<"unkown scan";
         }
      }
   }
   return m_M_XS;
}

std::map<int,float> getXsec_unc(Scan_t scan)
{
   std::string path=CMAKE_SOURCE_DIR;
   if (scan==TChiWg) path += "/xsec_N2C1_wino.csv";
   if (scan==TChiNg) path += "/xsec_comb_wino.csv";
   unsigned const nCol=3;
   std::ifstream fstream (path,std::ifstream::in);
   std::string line;
   std::map<int,float> m_M_XS;
   while (fstream.good()) {
      std::getline(fstream, line);
      if (line.find("#") == 0) continue;
      line=util::rm_duplicate_spaces(line);
      std::vector<float> values=util::to_vector<float>(line,' ');
      if (values.size() == nCol) {
         if (scan==TChiWg || scan==TChiNg) {
            m_M_XS[(int)values[0]]=values[2]/1000.; // convert to pb
         } else {
            debug<<"unkown scan";
         }
      }
   }
   return m_M_XS;
}

void plot_1d(Scan_t scan)
{
   assert(scan==TChiWg || scan==TChiNg);
   TString sScan;
   if (scan==TChiWg) sScan="TChiWg";
   if (scan==TChiNg) sScan="TChiNg";

   TString title, title2;
   if (scan==TChiWg)  title="pp #rightarrow"
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}, ";
   if (scan==TChiWg)  title2=
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #gamma/W^{#pm}#tilde{G}";
   if (scan==TChiNg)  title="pp #rightarrow"
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}/"
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}, ";
   if (scan==TChiNg)  title2=
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}#rightarrow"
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}+soft, "                        
                         "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #gamma/(H,Z)#tilde{G}";
   io::RootFileReader fileReader(TString::Format("uncert_%s.root",sScan.Data()));
   io::RootFileSaver saver("plots.root","limit_1d");
   std::map<TString,TGraph> gr,grISR,grMet;
   TCanvas can;
   for (TString lvl: {"statB1","statB2","statB3","statB4"}) {
      TGraph2D gr2(*fileReader.read<TGraph2D>("g_"+lvl));    
      gr[lvl]=TGraph(gr2.GetN());
      gr[lvl].SetLineWidth(2);
      double *x=gr2.GetX();
      double *z=gr2.GetZ();
      for (int i=0; i<gr2.GetN(); i++) {
         gr[lvl].SetPoint(i,x[i],z[i]);
      }
      gr[lvl].Sort();
      gr[lvl].GetYaxis()->SetRangeUser(0,20);
   }
   gr["statB1"].GetYaxis()->SetTitleSize(.055);
   gr["statB1"].GetYaxis()->SetTitleOffset(1.5);
   gr["statB1"].SetTitle(";m_{NLSP} [(GeV);#sigma_{stat} in %   ");
   gr["statB1"].Draw("a l");
   gr["statB2"].Draw("same l");
   gr["statB3"].Draw("same l");
   gr["statB4"].Draw("same l");
   gr["statB1"].SetLineColor(kBlue);
   gr["statB2"].SetLineColor(kBlue+2);
   gr["statB3"].SetLineColor(kBlue+4);
   gr["statB4"].SetLineColor(kBlue+8);
   gfx::LegendEntries le;
   le.append(gr["statB1"],"S_{T}^{#gamma}: 600-800 GeV","l");
   le.append(gr["statB2"],"S_{T}^{#gamma}: 800-1000 GeV","l");
   le.append(gr["statB3"],"S_{T}^{#gamma}: 1000-1300 GeV","l");
   le.append(gr["statB4"],"S_{T}^{#gamma}: 1300-#infty GeV","l");
   TLegend l1=le.buildLegend(.42,.6,.92,.92);
   l1.DrawClone();
   TLatex l=gfx::cornerLabel(title,3);
   l.SetY(0.85);
   TLatex l3=gfx::cornerLabel(title2,3);
   l3.SetY(0.8);

   l.DrawClone();
   l3.DrawClone();
   saver.save(can,sScan+"_uncert",true,false);
   
   TCanvas can2;
   can2.cd();
   for (TString lvl: {"ISRB1","ISRB2","ISRB3","ISRB4"}) {
      TGraph2D gr2(*fileReader.read<TGraph2D>("g_"+lvl));
      grISR[lvl]=TGraph(gr2.GetN()); 
      grISR[lvl].SetLineWidth(2);
      double *x=gr2.GetX();
      double *z=gr2.GetZ();
      for (int i=0; i<gr2.GetN(); i++) {
         grISR[lvl].SetPoint(i,x[i],z[i]);
      }      
      grISR[lvl].Sort();
      grISR[lvl].GetYaxis()->SetRangeUser(0,20);
   }
   grISR["ISRB1"].GetYaxis()->SetTitleSize(.055);
   grISR["ISRB1"].GetYaxis()->SetTitleOffset(1.5);
   grISR["ISRB1"].SetTitle(";m_{NLSP} (GeV);#sigma_{ISR} in %   ");
   grISR["ISRB1"].Draw("a l");
   grISR["ISRB2"].Draw("same l");
   grISR["ISRB3"].Draw("same l");
   grISR["ISRB4"].Draw("same l");
   grISR["ISRB1"].SetLineColor(kBlue);
   grISR["ISRB2"].SetLineColor(kBlue+2);
   grISR["ISRB3"].SetLineColor(kBlue+4);
   grISR["ISRB4"].SetLineColor(kBlue+8);
   l1.DrawClone();
   l.DrawClone();
   l3.DrawClone();
   saver.save(can2,sScan+"_uncert_ISR",true,false);

   TCanvas can3;
   can3.cd();
   for (TString lvl: {"MetB1","MetB2","MetB3","MetB4"}) {
      TGraph2D gr2(*fileReader.read<TGraph2D>("g_"+lvl));
      grMet[lvl]=TGraph(gr2.GetN()); 
      grMet[lvl].SetLineWidth(2);
      double *x=gr2.GetX();
      double *z=gr2.GetZ();
      for (int i=0; i<gr2.GetN(); i++) {
         grMet[lvl].SetPoint(i,x[i],z[i]);
      }      
      grMet[lvl].Sort();
      grMet[lvl].GetYaxis()->SetRangeUser(0,20);
   }
   grMet["MetB1"].GetYaxis()->SetTitleSize(.055);
   grMet["MetB1"].GetYaxis()->SetTitleOffset(1.5);
   grMet["MetB1"].SetTitle(";m_{NLSP} (GeV);#sigma_{FastSim p_{T}^{miss}} in %   ");
   grMet["MetB1"].Draw("a l");
   grMet["MetB2"].Draw("same l");
   grMet["MetB3"].Draw("same l");
   grMet["MetB4"].Draw("same l");
   grMet["MetB1"].SetLineColor(kBlue);
   grMet["MetB2"].SetLineColor(kBlue+2);
   grMet["MetB3"].SetLineColor(kBlue+4);
   grMet["MetB4"].SetLineColor(kBlue+8);
   l1.DrawClone();
   l.DrawClone();
   l3.DrawClone();
   saver.save(can3,sScan+"_uncert_Met",true,false);

}

extern "C"
void run()
{
  // plot_1d(TChiWg);
   plot_1d(TChiNg);
   
   
 //  save_aux(T5gg);
 //  save_aux(T5Wg);
 //  save_aux(GGM);
}
