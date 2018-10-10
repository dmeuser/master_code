#ifndef HIST_HPP__
#define HIST_HPP__

//#include "Config.hpp"
//#include "Dataset.hpp"
#include "gfx.h"

#include <map>

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>

namespace hist
{
   enum ErrorType
   {
      STAT,
      SYST,
      COMB,
      // special for ratios, residuals, ...
      ONLY1,
      ONLY2
   };

   std::vector<double> getBinVector(std::vector<float> edges, std::vector<float> widths);
   TH1F fromWidths(const char *name, const char *title,std::vector<float> edges, std::vector<float> widths);
   TH2F fromWidths_2d(const char *name, const char *title, std::vector<float> edges_x, std::vector<float> widths_x, std::vector<float> edges_y, std::vector<float> widths_y);


   TH1F rebinned(TH1F const &h, std::vector<float> const &edges, std::vector<float> const &widths,bool mergeOverflow=true,bool mergeUnderflow=true);
   TH1F rebinned(TH1F const &h, std::vector<double> const &binedges,bool mergeOverflow=true,bool mergeUnderflow=true);
   void divideByBinWidth(TH1& h,bool divideLastBin=true);
   //void mergeOverflow(TH1& h, bool includeUnderflow=true);
   //void mergeOverflow(TH2& h, bool includeUnderflow=true);

   void setMaximum(TH1& h,std::vector<TH1F> hists,float multiplier=1.1);
   void setMinimum(TH1& h,std::vector<TH1F> hists,float multiplier=0.9,bool allowNegative=true);

   TH1F getRatio   (TH1F const &h1,TH1F const &h2,TString title="ratio",ErrorType et=COMB);
   TH1F getResidual(TH1F const &h1,TH1F const &h2,TString title="residual",ErrorType et=COMB);
   TH1F getPull    (TH1F const &h1,TH1F const &h2,TString title="pull",ErrorType et=COMB);
   TH1F getRatio   (TH1F const &h1,THStack &h2,TString title="ratio",ErrorType et=COMB);
   TH1F getPull    (TH1F const &h1,THStack &h2,TString title="pull",ErrorType et=COMB);

   TGraphErrors getPullGraph(TH1F const &h1,TH1F const &h2,TString title="ratio",ErrorType et=COMB);
   TGraphErrors getRatioGraph(TH1F const &h1,TH1F const &h2,TString title="ratio",ErrorType et=COMB);
   TGraphErrors getRatioGraph(TH1F const &h1,THStack &h2,TString title="ratio",ErrorType et=COMB);

   THStack stackPrepend(THStack const& stOld, TH1F &h, Option_t *option="");
}

#endif /* HIST_HPP__ */
