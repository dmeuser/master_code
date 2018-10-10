#include "hist.h"

//#include "tools/io.hpp"
#include "gfx.h"
//#include "tools/util.hpp"

#include <TMath.h>
#include <cassert>
#include <iostream>


std::vector<double> hist::getBinVector(std::vector<float> edges, std::vector<float> widths)
{
   assert(edges.size()==widths.size()+1);
   std::vector<double> xbins={edges[0]};
   for (uint i=0; i<edges.size()-1; i++){
      float x=edges[i];
      while ((edges[i+1]-x)>1e-5){
         x+=widths[i];
         xbins.push_back(x);
      }
   }
   return xbins;
}

TH1F hist::fromWidths(const char *name, const char *title,std::vector<float> edges, std::vector<float> widths)
{
   std::vector<double> xbins=getBinVector(edges, widths);
   return TH1F(name,title,xbins.size()-1,&xbins[0]);
}

TH2F hist::fromWidths_2d(const char *name, const char *title, std::vector<float> edges_x, std::vector<float> widths_x, std::vector<float> edges_y, std::vector<float> widths_y)
{
   std::vector<double> xbins=getBinVector(edges_x, widths_x);
   std::vector<double> ybins=getBinVector(edges_y, widths_y);
   return TH2F(name,title,xbins.size()-1,&xbins[0],ybins.size()-1,&ybins[0]);
}

TH1F hist::rebinned(TH1F const &h, std::vector<float> const &edges, std::vector<float> const &widths,bool mergeOverflow,bool mergeUnderflow)
{
   std::vector<double> binedges=getBinVector(edges, widths);
   return rebinned(h,binedges,mergeOverflow,mergeUnderflow);
}

TH1F hist::rebinned(TH1F const &h, std::vector<double> const &binedges,bool mergeOverflow,bool mergeUnderflow)
{
   TH1F hClone(h);
   std::string name(hClone.GetName());
   name+="_rebinned";
   TH1F *hnew=(TH1F*)hClone.Rebin(binedges.size()-1,name.c_str(),&binedges[0]);
   //if (mergeOverflow) hist::mergeOverflow(*hnew,mergeUnderflow);
   TString yTitle=hClone.GetYaxis()->GetTitle();
   yTitle.ReplaceAll("BIN"," / bin");
   hnew->GetYaxis()->SetTitle(yTitle);
   return *hnew;
}

void hist::divideByBinWidth(TH1& h,bool divideLastBin)
{
   int N=h.GetNbinsX();
   if (h.GetBinContent(N+1) != 0) {
      std::cerr<<"non-emtpy overflow. merge first!";
      throw;
   }
   float w;
   if (divideLastBin) N++;
   //else: not dividing the last bin=merged overflow
   for (int i=0; i<N; i++) {
      w=h.GetBinWidth(i);
      h.SetBinContent(i,h.GetBinContent(i)/w);
      h.SetBinError(i,h.GetBinError(i)/w);
   }
   // labels
   TString yt=h.GetYaxis()->GetTitle();
   yt.ReplaceAll("BIN","BINW");
   yt.ReplaceAll(" / bin","BINW");
   h.GetYaxis()->SetTitle(yt);
}

/*
void hist::mergeOverflow(TH1& h, bool includeUnderflow)
{
   int N=h.GetNbinsX();
   // -- overflow
   float cont=h.GetBinContent(N)+h.GetBinContent(N+1);
   float err2=util::quadSum<double>({h.GetBinError(N),h.GetBinError(N+1)});
   // set content+error of last bin
   h.SetBinContent(N,cont);
   h.SetBinError(N,TMath::Sqrt(err2));
   // clear overflow
   h.SetBinContent(N+1,0);
   h.SetBinError(N+1,0);
   if (!includeUnderflow) return;
   // -- underflow
   cont=h.GetBinContent(0)+h.GetBinContent(1);
   err2=util::quadSum<double>({h.GetBinError(0),h.GetBinError(1)});
   // set content+error of first bin
   h.SetBinContent(1,cont);
   h.SetBinError(1,TMath::Sqrt(err2));
   // clear overflow
   h.SetBinContent(0,0);
   h.SetBinError(0,0);
}

void hist::mergeOverflow(TH2& h, bool includeUnderflow)
{
   int N_X=h.GetNbinsX();
   int N_Y=h.GetNbinsY();
   includeUnderflow = true;
   //loop over Y-bins
   for (int i=0; i<=N_Y+1; i++){
      // -- overflow
      float cont=h.GetBinContent(N_X,i)+h.GetBinContent(N_X+1,i);
      float err2=util::quadSum<double>({h.GetBinError(N_X,i),h.GetBinError(N_X+1,i)});
      // set content+error of last bin
      h.SetBinContent(N_X,i,cont);
      h.SetBinError(N_X,i,TMath::Sqrt(err2));
      // clear overflow
      h.SetBinContent(N_X+1,i,0);
      h.SetBinError(N_X+1,i,0);
      if (!includeUnderflow) return;
      // -- underflow
      cont=h.GetBinContent(0,i)+h.GetBinContent(1,i);
      err2=util::quadSum<double>({h.GetBinError(0,i),h.GetBinError(1,i)});
      // set content+error of first bin
      h.SetBinContent(1,i,cont);
      h.SetBinError(1,i,TMath::Sqrt(err2));
      // clear overflow
      h.SetBinContent(0,i,0);
      h.SetBinError(0,i,0);
   }
   //Loop over X-bins
   for (int i=0; i<=N_X+1; i++){
      // -- overflow
      float cont=h.GetBinContent(i,N_Y)+h.GetBinContent(i,N_Y+1);
      float err2=util::quadSum<double>({h.GetBinError(i,N_Y),h.GetBinError(i,N_Y+1)});
      // set content+error of last bin
      h.SetBinContent(i,N_Y,cont);
      h.SetBinError(i,N_Y,TMath::Sqrt(err2));
      // clear overflow
      h.SetBinContent(i,N_Y+1,0);
      h.SetBinError(i,N_Y+1,0);
      if (!includeUnderflow) return;
      // -- underflow
      cont=h.GetBinContent(i,0)+h.GetBinContent(i,1);
      err2=util::quadSum<double>({h.GetBinError(i,0),h.GetBinError(i,1)});
      // set content+error of first bin
      h.SetBinContent(i,1,cont);
      h.SetBinError(i,1,TMath::Sqrt(err2));
      // clear overflow
      h.SetBinContent(i,0,0);
      h.SetBinError(i,0,0);
   }
}
*/

void hist::setMaximum(TH1& h,std::vector<TH1F> hists,float multiplier)
{
   double max=h.GetMaximum();
   for (TH1 const &hh: hists) max=TMath::Max(max,hh.GetMaximum());
   h.SetMaximum(max*multiplier);
}

void hist::setMinimum(TH1& h,std::vector<TH1F> hists,float multiplier, bool allowNegative)
{
   double min=h.GetMaximum();
   hists.push_back(TH1F(*(TH1F*)&h));
   double v;
   for (TH1 const &hh: hists) {
      for (int i=0; i<=hh.GetNbinsX()+1; i++){
         v=hh.GetBinContent(i);
         if (allowNegative || v>0) min=TMath::Min(min,v);
      }
   }
   h.SetMinimum(min*multiplier);
}

TH1F hist::getRatio(TH1F const &h1,TH1F const &h2,TString title,ErrorType et)
{
   TH1F hRatio(h1);
   hRatio.Divide(&h2);
   float N2,err;
   if (et==ONLY1){
      for (int i=0; i<=h1.GetNbinsX()+1;i++){
         N2 = h2.GetBinContent(i);
         err = (N2==0.0) ? 0.0 : h1.GetBinError(i)/N2;
         hRatio.SetBinError(i,err);
      }
   } else if (et==ONLY2){
      for (int i=0; i<=h1.GetNbinsX()+1;i++){
         N2 = h2.GetBinContent(i);
         err = (N2==0.0) ? 0.0 : h2.GetBinError(i)/N2;
         hRatio.SetBinError(i,err);
      }
   } else assert(et==COMB);

   Double_t min = hRatio.GetBinContent(hRatio.GetMinimumBin());
   Double_t max = hRatio.GetBinContent(hRatio.GetMaximumBin());
   min=TMath::Min(min,0.5);
   max=TMath::Max(max,1.5);
   hRatio.GetYaxis()->SetRangeUser(min,max*1.1);
   hRatio.GetYaxis()->SetTitle(title);

   if (et==ONLY2 || et==COMB){
      hRatio.SetLineColor(kBlack);
      hRatio.SetLineWidth(2);
      hRatio.SetMarkerStyle(1);
      hRatio.SetFillStyle(1001);
      hRatio.SetFillColor(kGray);
   }

   return hRatio;
}

TH1F hist::getResidual(TH1F const &h1,TH1F const &h2,TString title,ErrorType et)
{
   assert(et==COMB); // TODO implement COMB1/2
   TH1F hRes(h1);
   hRes.Add(&h2,-1);
   Double_t min = hRes.GetBinContent(hRes.GetMinimumBin());
   Double_t max = hRes.GetBinContent(hRes.GetMaximumBin());
   hRes.GetYaxis()->SetRangeUser(min,max*1.1);
   hRes.GetYaxis()->SetTitle(title);
   hRes.SetLineColor(kBlack);
   hRes.SetLineWidth(2);
   hRes.SetMarkerStyle(1);
   return hRes;
}

TH1F hist::getPull(TH1F const &h1,TH1F const &h2,TString title,ErrorType et)
{
   TH1F hPull=getResidual(h1,h2,title);
   if (et==COMB) {
     for (int i=0; i<=h1.GetNbinsX()+1;i++){
       float const cont=hPull.GetBinContent(i);
       float const err=hPull.GetBinError(i);
       if (std::abs(cont)<1e-6 && err<1e-6){
          hPull.SetBinContent(i,0);
          hPull.SetBinError(i,0);
       } else {
          if (err<1e-6){
             std::cerr<<"Zero error in non-emtpy bin!";
             throw;
          }
          hPull.SetBinContent(i,cont/err);
          hPull.SetBinError(i,1.);
       }
      }
   } else if (et==ONLY1) {
     for (int i=0; i<=h1.GetNbinsX()+1;i++){
          hPull.SetBinError(i,h1.GetBinError(i));
       }      
   }

   Double_t max = TMath::Max(
      std::abs(hPull.GetBinContent(hPull.GetMinimumBin())),
      std::abs(hPull.GetBinContent(hPull.GetMaximumBin())));
   max=TMath::Min(max,10.);
   max=(int)max+1.4;
   hPull.GetYaxis()->SetRangeUser(-max,max);
   hPull.SetFillStyle(1001);
   hPull.SetLineWidth(1);
   hPull.SetLineColor(kBlack);
   hPull.SetFillColor(kGray);
   return hPull;
}

TH1F hist::getRatio(TH1F const &h1,THStack &h2,TString title,ErrorType et)
{
   return getRatio(h1,*(TH1F*)h2.GetStack()->Last(),title,et);
}

TH1F hist::getPull(TH1F const &h1,THStack &h2,TString title,ErrorType et)
{
   return getPull(h1,*(TH1F*)h2.GetStack()->Last(),title,et);
}

TGraphErrors hist::getPullGraph(TH1F const &h1,TH1F const &h2,TString title,ErrorType et)
{
   TH1F h = *(TH1F*)(h1.Clone());
   for (int i=0; i<=h.GetNbinsX()+1;i++){
     if (h2.GetBinContent(i)) 
       h.SetBinError(i,h1.GetBinContent(i)/h2.GetBinContent(i) );
     h.SetBinContent(i,0);
   }  
   return TGraphErrors(&h);
}

TGraphErrors hist::getRatioGraph(TH1F const &h1,TH1F const &h2,TString title,ErrorType et)
{
   TH1F h = getRatio(h1,h2,title,et);
   return TGraphErrors(&h);
}

TGraphErrors hist::getRatioGraph(TH1F const &h1,THStack &h2,TString title,ErrorType et)
{
   TH1F h = getRatio(h1,h2,title,et);
   return TGraphErrors(&h);
}

THStack hist::stackPrepend(THStack const& stOld, TH1F &h, Option_t *option){
   THStack st;
   st.SetTitle(stOld.GetTitle());
   st.Add(&h,option);
   for (int i=0; i<stOld.GetNhists(); i++){
      st.Add((TH1*)stOld.GetHists()->At(i),option);
   }
   return st;
}
