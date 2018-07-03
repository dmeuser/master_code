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

void plot_1d(Scan_t scan, TString selection)
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
   std::map<int,float> mXsecs=getXsecs(scan);
   std::map<int,float> mXsec_unc=getXsec_unc(scan);

   io::RootFileReader fileReader(TString::Format("../input/limits/limits_TChiNg_%s.root",selection.Data()));
   io::RootFileSaver aux_saver(TString::Format("aux_limits_%s.root",sScan.Data()),selection.Data());
   std::map<TString,TGraph> gr;
   for (TString lvl: {"exp","exp+1","exp-1","exp+2","exp-2","obs"}) {
      TGraph2D gr2(*fileReader.read<TGraph2D>("gr_"+lvl));    
      gr[lvl]=TGraph(gr2.GetN());
      gr[lvl].SetLineWidth(2);
      double *x=gr2.GetX();
      double *z=gr2.GetZ();
      for (int i=0; i<gr2.GetN(); i++) {
         gr[lvl].SetPoint(i,x[i],z[i]*mXsecs[int(x[i])]);
      }
      gr[lvl].Sort();
      aux_saver.save(gr[lvl],lvl);
   }
   //~ TGraph2D grSig(*fileReader.read<TGraph2D>("gr_significance"));
   //~ TGraph gr1DSig=TGraph(grSig.GetN());
   //~ double *xSig=grSig.GetX();
   //~ double *zSig=grSig.GetZ();
   //~ for (int i=0; i<grSig.GetN(); i++) {
      //~ gr1DSig.SetPoint(i,xSig[i],zSig[i]);
   //~ }
   //~ gr1DSig.Sort();

   for (TString lvl: {"1","2"}) {
      gr[lvl]=TGraph(gr["exp+"+lvl]);
      double x, y;
      for (int i=gr["exp-"+lvl].GetN()-1; i>=0; i--) {
         gr["exp-"+lvl].GetPoint(i,x,y);
         gr[lvl].SetPoint(gr[lvl].GetN(),x,y);
      }
      gr[lvl].GetPoint(0,x,y);
      gr[lvl].SetPoint(gr[lvl].GetN(),x,y);
   }
   gr["xs"]=TGraph(mXsecs.size());
   gr["xs"].SetLineWidth(2);
   gr["xs+1"]=TGraph(gr["xs"]);
   gr["xs-1"]=TGraph(gr["xs"]);
   int iP=0;
   for (auto m_xs :mXsecs) {
      int const m=m_xs.first;
      gr["xs"].SetPoint(iP,m,m_xs.second);
      gr["xs+1"].SetPoint(iP,m,m_xs.second+mXsec_unc[m]);
      gr["xs-1"].SetPoint(iP,m,m_xs.second-mXsec_unc[m]);
      iP++;
   }

   aux_saver.save(gr["xs"],"xs");
   aux_saver.save(gr["xs+1"],"xs+1");
   aux_saver.save(gr["xs-1"],"xs-1");

   io::RootFileSaver saver("plots.root",TString::Format("limit_1d/%s",selection.Data()));
   TCanvas can;
   can.SetLogy();

   gr["2"].Draw("afl");
   gr["2"].SetFillColor(kYellow-7);
   gr["2"].SetLineColor(kYellow-7);
   gr["1"].Draw("fl");
   gr["1"].SetFillColor(kGreen-7);
   gr["1"].SetLineColor(kGreen-7);
   gr["exp"].Draw("l");
   gr["exp"].SetLineStyle(2);
   gr["exp"].SetLineColor(kRed);
   gr["obs"].Draw("l");
   gr["obs"].SetLineColor(kBlack);
   gr["xs"].Draw("l");
   gr["xs"].SetLineColor(kBlue);
   gr["xs"].SetLineStyle(2);
   gr["xs+1"].Draw("l");
   gr["xs+1"].SetLineColor(kBlue);
   gr["xs-1"].Draw("l");
   gr["xs-1"].SetLineColor(kBlue);

   if (scan == TChiNg) gr["2"].GetXaxis()->SetRangeUser(300,1150);
   if (scan == TChiWg) gr["2"].GetXaxis()->SetRangeUser(300,1000);
   if (scan == TChiNg) gr["2"].GetYaxis()->SetRangeUser(0.0004,0.8);
   if (scan == TChiWg) gr["2"].GetYaxis()->SetRangeUser(0.0008,0.4);
   gr["2"].GetYaxis()->SetTitleSize(.055);
   gr["2"].GetYaxis()->SetTitleOffset(1.5);
   gr["2"].SetTitle(";m_{NLSP} (GeV);95% CL upper limit / cross section (pb)");

   gfx::LegendEntries le;
   le.append(gr["xs"],"Cross section","l");
   le.append(gr["xs+1"],"Theory uncert.","l");
   TLegend l1=le.buildLegend(.55,.82,.97,.92,1);
   l1.SetColumnSeparation(-0.2);
   l1.DrawClone();
   le.clear();
   le.append(gr["obs"],"Observed limit","l");
   le.append(gr["exp"],"Expected limit","l");
   le.buildLegend(.55,.72,.97,.82,1).DrawClone();
   le.clear();
   le.append(gr["1"],"68% expected","f");
   le.append(gr["2"],"95% expected","f");
   TLegend l2=le.buildLegend(.55,.62,.97,.72,1);
   l2.SetColumnSeparation(-0.2);
   l2.DrawClone();


   TLatex l=gfx::cornerLabel(title,3);
   l.SetY(0.24);
   l.SetX(0.2);
   l.SetTextSize(0.048);
   l.DrawClone();
   TLatex l3=gfx::cornerLabel(title2,3);
   l3.SetY(0.17);
   l3.SetX(0.2);
   l3.SetTextSize(0.048);
   l3.DrawClone();

   can.RedrawAxis();
   saver.save(can,sScan,true,false);

   //~ TCanvas can2;
   //~ gr1DSig.Draw("al");
   //~ gr1DSig.SetLineColor(kBlue);
   //~ gr1DSig.SetLineWidth(2);
   //~ gr1DSig.SetTitle(";m_{NLSP} (GeV);Significance (s.d.)");
   //~ gr1DSig.GetYaxis()->SetTitleOffset(0.8);
   //~ if (scan == TChiNg) gr1DSig.GetXaxis()->SetRangeUser(300,1150);
   //~ if (scan == TChiWg) gr1DSig.GetXaxis()->SetRangeUser(300,1000);
   //~ if (scan == TChiNg) gr1DSig.GetYaxis()->SetRangeUser(-3,3);
   //~ if (scan == TChiWg) gr1DSig.GetYaxis()->SetRangeUser(-3,3);
   //~ l.DrawClone();
   //~ l3.DrawClone();
   //~ can2.RedrawAxis();
   //~ saver.save(can2,sScan+"_significance",true,false);
}

extern "C"
void run()
{  
   //~ TString selection="inclusiv";
   TString selection="exclusiv";
   //~ plot_1d(TChiWg);
   plot_1d(TChiNg,selection);
}
