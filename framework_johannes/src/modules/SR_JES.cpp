#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/io.hpp"
#include "tools/util.hpp"
#include "tools/weighters.hpp"
#include "tools/physics.hpp"

#include <limits>

#include <TFile.h>
#include <TFractionFitter.h>
#include <TGaxis.h>
#include <TMarker.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TStyle.h>

#include <RooAddPdf.h>
#include <RooDataHist.h>
#include <RooEllipse.h>
#include <RooFitResult.h>
#include <RooHistPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>

Config const &cfg=Config::get();
static io::RootFileSaver saver("plots.root","SR_JES");
static io::Logger logFile("SFapplication.log");

int nWeights = 100;
const int nPDF = 91;
float nominalVg1, nominalVg2, nominalVg3, nominalVg4;
float nominalgJ1, nominalgJ2, nominalgJ3, nominalgJ4;
float bin1_Vg[nPDF], bin2_Vg[nPDF], bin3_Vg[nPDF], bin4_Vg[nPDF];
float bin1_gJ[nPDF], bin2_gJ[nPDF], bin3_gJ[nPDF], bin4_gJ[nPDF];
io::Logger SR_Vg_JES_effect(("SR_Vg_JES_variance.tex"));
io::Logger SR_gJ_JES_effect(("SR_gJ_JES_variance.tex"));
io::Logger JES_SR_Vg_graph(("JES_SR_Vg_graph.tex"));
io::Logger JES_SR_gJ_graph(("JES_SR_gJ_graph.tex"));

// Computes the standard deviation
double ComputeStdDev(int n, float *x, float cv)
{
  double sum = 0;
  for (int i = 0; i < n; i++)
    sum += pow(x[i] - cv, 2.0);

  sum /= n-1;
  return sqrt(sum);
}

class Rebinner
{
public:
   // if iRebin==0: use vectors
   Rebinner(int iRebin, bool bDivideByBinWidth=false)
      : iRebin_(iRebin)
      , bDivideByBinWidth_(bDivideByBinWidth)
      {}
   Rebinner(std::vector<float> rebinEdges,std::vector<float> rebinWidths, bool bDivideByBinWidth=false)
      : iRebin_(0)
      , bDivideByBinWidth_(bDivideByBinWidth)
      {
         newBinEdges_=hist::getBinVector(rebinEdges,rebinWidths);
      }
   TH1F operator()(TH1F const &h) const {
      TH1F hNew;
      if (iRebin_!=0){
         hNew=h;
         hNew.Rebin(iRebin_);
      } else {
         hNew=TH1F(h);
         hNew=hist::rebinned(h,newBinEdges_,true,false); // merging overflow, not underflow
      }
      if (bDivideByBinWidth_) {
         hist::divideByBinWidth(hNew);
      }
      return hNew;
   }
   int numberOfBins(){
      if (iRebin_==0) return newBinEdges_.size()-1;
      return 0;
   }
   float binEdge(int i) { // using root bin numbers
      assert(iRebin_==0);
      i--; // using root bin numbers
      if (i==-1) return -std::numeric_limits<float>::infinity();
      if ((unsigned)i==newBinEdges_.size()) return std::numeric_limits<float>::infinity();
      return newBinEdges_.at(i);
   }
private:
   int iRebin_;
   bool bDivideByBinWidth_;
   std::vector<double> newBinEdges_;
};


TH1F combineHists(io::RootFileReader const& histReader, TString sPresel, TString sVar, std::vector<TString> sampleNames)
{
   assert(sampleNames.size()>0);
   TH1F const *hFirst=histReader.read<TH1F>(sPresel+sVar+"/"+sampleNames[0]);
   TH1F hComb(*hFirst);
   hComb.SetLineColor(hFirst->GetLineColor());
   hComb.SetFillColor(hFirst->GetLineColor());
   for (unsigned i=1; i<sampleNames.size(); i++){
      hComb.Add(histReader.read<TH1F>(sPresel+sVar+"/"+sampleNames[i]));
   }
   return hComb;
}

void setupHistForStack(TH1 &h){
   h.SetFillStyle(1001);
   h.SetFillColor(h.GetLineColor());
   h.SetLineColor(kBlack);
   h.SetLineWidth(1);
}

struct SignalBin
{
   struct Component
   {
      Component(){}
      float count=0;
      float estat=0;
      float esyst=0;
      TString formatLine(float nTotalBkg=-1,bool tex=false, bool data=false,float add=-1) const {
         if (nTotalBkg<0) nTotalBkg=count;
         float esy=esyst;
         if (add>0) esy=TMath::Sqrt(esy*esy+add*add*count*count);
         if (!tex) return TString::Format("%6.2f %6.2f %6.2f(%4.1f%%)",
                                          count,estat,esy,esy/nTotalBkg*100);
         TString countFormat=" & %6.2f";
         if (data) countFormat=" & %6.0f";
         if (esy>0 && estat>0) return TString::Format(countFormat+" & %6.2f & %6.2f & %4.1f\\%%\\\\",
                                                        count,estat,esy,esy/nTotalBkg*100);
         if (estat>0) return TString::Format(countFormat+" & %6.2f &  & \\\\",
                                             count,estat);
         return TString::Format(countFormat+" &  &  & \\\\",count);
      }
   };
   float lower,upper; // edges
   std::map<TString,Component> component;

   TString formatBlock(bool tex=false,bool showTotal=false) {
      float nBkgTot=getTotal().count;
      TString block;
      for (TString compName: {"Vg","GJ","TTcomb","efake","diboson","TChiWG","T5Wg"}) {
         TString line;
         if (component.count(compName)) {
            if (compName=="GJ") line=component[compName].formatLine(nBkgTot,tex,false,cfg.sf.uncert_gammaJ);
            if (compName=="Vg") line=component[compName].formatLine(nBkgTot,tex,false,cfg.sf.uncert_Vgamma);           
            else line=component[compName].formatLine(nBkgTot,tex);
         } else {
            continue;
         }
         if (tex) {
            compName.ReplaceAll("5","five");
            compName="\\"+compName;
         }
         block+=TString::Format("   %-10s ",compName.Data());
         block+=line;
         block+="\n";
      }
      if (showTotal) {
         block+=TString::Format("   %-10s ",tex?"\\total":"total");
         block+=getTotal().formatLine(-1,tex);
      }
      if (component.count("data")) {
         block+=TString::Format("\n   %-10s ","data");
         block+=component["data"].formatLine(-1,tex,true);
      }
      return block;
   }
   Component getTotal() const {
      Component tot;
      for (auto const &b: component) {
         if (b.first=="data") continue;
         tot.count+=b.second.count;
         tot.estat+=b.second.estat*b.second.estat;
         tot.esyst+=b.second.esyst*b.second.esyst;
      }
      if (component.count("Vg") && component.count("GJ")) {
         // background: consider correlation
         tot.esyst+=2*cfg.sf.rho
            *component.at("Vg").esyst
            *component.at("GJ").esyst;
         tot.esyst+=(cfg.sf.uncert_gammaJ*component.at("GJ").count*cfg.sf.uncert_gammaJ*component.at("GJ").count);
         tot.esyst+=(cfg.sf.uncert_Vgamma*component.at("Vg").count*cfg.sf.uncert_Vgamma*component.at("Vg").count);         
      }
      tot.estat=TMath::Sqrt(tot.estat);
      tot.esyst=TMath::Sqrt(tot.esyst);
      return tot;
   }
};

std::map<int,float> getVgScale()
{
   std::string path=CMAKE_SOURCE_DIR;
   path += "/output/JES_fitresult.tex";
   unsigned const nCol= 3;
   std::ifstream fstream (path,std::ifstream::in);
   std::string line;
   std::map<int,float> mVg_scale;
   while (fstream.good()) {
      std::getline(fstream, line);
      if (line.find("#") == 0) continue;
      if (line.find("%") == 0) continue;
      line=util::rm_duplicate_spaces(line);
      std::vector<float> values=util::to_vector<float>(line,' ');
      if (values.size() == nCol) {
         mVg_scale[(int)values[0]] = values[1];
      }
   }
   return mVg_scale;
}

std::map<int,float> getgJScale()
{
   std::string path=CMAKE_SOURCE_DIR;
   path += "/output/JES_fitresult.tex";
   unsigned const nCol= 3;
   std::ifstream fstream (path,std::ifstream::in);
   std::string line;
   std::map<int,float> mgJ_scale;
   while (fstream.good()) {
      std::getline(fstream, line);
      if (line.find("#") == 0) continue;
      if (line.find("%") == 0) continue;
      line=util::rm_duplicate_spaces(line);
      std::vector<float> values=util::to_vector<float>(line,' ');
      if (values.size() == nCol) {
         mgJ_scale[(int)values[0]] = values[2];
      }
   }
   return mgJ_scale;
}

enum PlotMode_t
{
   PRE, // preselection
   CR, // control region
   SR_BLIND, // signal region without data
   SR_CHECK, // signal region control distributions with data
   SR, // signal region unblinded
   VR,
};

// General implementation of plot. Call on of the functions below, not this one.
void plot(TString sSelection,TString sVar,int iRebin,
          std::vector<float> rebinEdges,std::vector<float> rebinWidths,
          PlotMode_t plotMode, int weight_index
   )
{
   bool const showData=(plotMode!=SR_BLIND && plotMode!=PRE);
   bool const showSignal=(plotMode!=CR && plotMode!=VR);
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("JES_studies%.1f",cfg.processFraction*100));

   TString saveName=sSelection+sVar;

   bool divideByBinWidth=(plotMode==CR || plotMode==VR || plotMode==PRE);
   Rebinner rebinned(iRebin,divideByBinWidth);
   if (iRebin==0) rebinned=Rebinner(rebinEdges,rebinWidths,divideByBinWidth);

   bool const customBinning = (iRebin==0 && plotMode!=CR && plotMode!=VR);
   if (customBinning) saveName+="_binning";
   if (!showData)     saveName+="_blind";

   std::vector<SignalBin> signalBins;
   std::vector<SignalBin> signalBins_signal;
   SignalBin bkgIntegral; // for the whole signal region
   SignalBin signalIntegral;
   if (plotMode==SR) {
      for (int i=0; i<rebinned.numberOfBins()+2; i++) {
         SignalBin bin;
         bin.lower=rebinned.binEdge(i);
         bin.upper=rebinned.binEdge(i+1);
         signalBins.push_back(bin);
         signalBins_signal.push_back(bin);
      }
      bkgIntegral.lower=bkgIntegral.upper=0;
      signalIntegral.lower=bkgIntegral.upper=0;
   }

   std::vector<TString> vsVgComponents{"WGToLNuG","ZNuNuGJets","ZGTo2LG","ZNuNuJets","WLNuJets"};
   std::vector<TString> vsGJComponents{"GJets_DR"};

   std::map<int,float> mVg_scale=getVgScale();
   std::map<int,float> mgJ_scale=getgJScale();

   TH1F hVg=rebinned(combineHists(histReader,sSelection,sVar,vsVgComponents));
   TH1F hGJ=rebinned(combineHists(histReader,sSelection,sVar,vsGJComponents));

   if (weight_index == 0 ) {
      hVg.Scale(cfg.trigger_eff_Ph*mVg_scale[0]);
      hGJ.Scale(cfg.trigger_eff_Ph*mgJ_scale[0]);
   
      nominalVg1 = hVg.GetBinContent(1);
      nominalVg2 = hVg.GetBinContent(2);
      nominalVg3 = hVg.GetBinContent(3);
      nominalVg4 = hVg.GetBinContent(4);
      nominalgJ1 = hGJ.GetBinContent(1);
      nominalgJ2 = hGJ.GetBinContent(2);
      nominalgJ3 = hGJ.GetBinContent(3);
      nominalgJ4 = hGJ.GetBinContent(4);
   }
   else {
   hVg.Scale(cfg.trigger_eff_Ph*mVg_scale[weight_index]);
   hGJ.Scale(cfg.trigger_eff_Ph*mgJ_scale[weight_index]);
   }

   float VarVg1, VarVg2, VarVg3, VarVg4;
   float VargJ1, VargJ2, VargJ3, VargJ4;
      
   VarVg1 = hVg.GetBinContent(1) - nominalVg1;
   VarVg2 = hVg.GetBinContent(2) - nominalVg2;
   VarVg3 = hVg.GetBinContent(3) - nominalVg3;
   VarVg4 = hVg.GetBinContent(4) - nominalVg4;
   VargJ1 = hGJ.GetBinContent(1) - nominalgJ1;
   VargJ2 = hGJ.GetBinContent(2) - nominalgJ2;
   VargJ3 = hGJ.GetBinContent(3) - nominalgJ3;
   VargJ4 = hGJ.GetBinContent(4) - nominalgJ4;

   TString Vg_variance = std::to_string(weight_index) +
                        " " + std::to_string(VarVg1) + " " + std::to_string(VarVg1/nominalVg1) +
                        " " + std::to_string(VarVg2) + " " + std::to_string(VarVg2/nominalVg2) +
                        " " + std::to_string(VarVg3) + " " + std::to_string(VarVg3/nominalVg3) +
                        " " + std::to_string(VarVg4) + " " + std::to_string(VarVg4/nominalVg4);
   SR_Vg_JES_effect << Vg_variance;

   TString gJ_variance = std::to_string(weight_index) +
                        " " + std::to_string(VargJ1) + " " + std::to_string(VargJ1/nominalgJ1) +
                        " " + std::to_string(VargJ2) + " " + std::to_string(VargJ2/nominalgJ2) +
                        " " + std::to_string(VargJ3) + " " + std::to_string(VargJ3/nominalgJ3) +
                        " " + std::to_string(VargJ4) + " " + std::to_string(VargJ4/nominalgJ4);
   SR_gJ_JES_effect << gJ_variance;

   TString Vg_1 = "700 "+ std::to_string(VarVg1/nominalVg1);
   TString Vg_2 = "900 "+ std::to_string(VarVg2/nominalVg2);
   TString Vg_3 = "1150 "+ std::to_string(VarVg3/nominalVg3);
   TString Vg_4 = "1450 "+ std::to_string(VarVg4/nominalVg4);

   TString gJ_1 = "700 "+ std::to_string(VargJ1/nominalgJ1);
   TString gJ_2 = "900 "+ std::to_string(VargJ2/nominalgJ2);
   TString gJ_3 = "1150 "+ std::to_string(VargJ3/nominalgJ3);
   TString gJ_4 = "1450 "+ std::to_string(VargJ4/nominalgJ4);

   JES_SR_Vg_graph << Vg_1;
   JES_SR_Vg_graph << Vg_2;
   JES_SR_Vg_graph << Vg_3;
   JES_SR_Vg_graph << Vg_4;
   
   JES_SR_gJ_graph << gJ_1;
   JES_SR_gJ_graph << gJ_2;
   JES_SR_gJ_graph << gJ_3;
   JES_SR_gJ_graph << gJ_4;
        
   io::log<<TString::Format("N_Vg=%.2f",hVg.Integral(0,hVg.GetNbinsX()+1,divideByBinWidth ? "width" : ""));
   io::log<<TString::Format("N_GJ=%.2f",hGJ.Integral(0,hGJ.GetNbinsX()+1,divideByBinWidth ? "width" : ""));

   hVg.SetLineColor(kOrange-3);
   hGJ.SetLineColor(kAzure+10);
   setupHistForStack(hVg);
   setupHistForStack(hGJ);

   gfx::LegendEntries le;
   THStack st;
   std::map<TString,TH1F> fixHists;
   for (TString sSample:{"efake","TTJets","TTGJets","diboson","TChiWG","T5Wg"}){
      fixHists[sSample]=rebinned(*(TH1F*)histReader.read<TH1F>(sSelection+sVar+"/"+sSample));
      TH1 &h=fixHists[sSample];
      h.Scale(cfg.trigger_eff_Ph);
      setupHistForStack(h);
      if (sSample == "efake"){
         h.SetFillColor(kMagenta+2);
      }
      else if (sSample == "TTGJets"){
         h.SetFillColor(kBlue+2);
      }
      else if (sSample == "diboson"){
         h.SetFillColor(kMagenta+4);
      }
   }
   fixHists["efake"].Scale(1./cfg.trigger_eff_Ph); // efake is data-driven

   fixHists["TTcomb"]=fixHists["TTGJets"];
   fixHists["TTcomb"].Add(&fixHists["TTJets"]);
   st.Add(&fixHists["diboson"],"hist");
   st.Add(&fixHists["efake"],"hist");
   st.Add(&fixHists["TTcomb"],"hist");
   le.prepend(fixHists["diboson"],cfg.datasets.getLabel("diboson"),"f");
   le.prepend(fixHists["efake"],cfg.datasets.getLabel("efake"),"f");
   le.prepend(fixHists["TTcomb"],"t#bar{t}(+#gamma)","f");
  

   if (plotMode==CR) {
      // CR: more GJet
      st.Add(&hVg,"hist");
      st.Add(&hGJ,"hist");
      le.prepend(hVg,"V(+#gamma)","f");
      le.prepend(hGJ,"#gamma+jets","f");
   } else {
      // SR,VR: more Vg
      st.Add(&hGJ,"hist");
      st.Add(&hVg,"hist");
      le.prepend(hGJ,"#gamma+jets","f");
      le.prepend(hVg,"V(+#gamma)","f");
   }
   TH1F hStackSum(*(TH1F*)st.GetStack()->Last());
   for (TString sSample:{"TChiWG","T5Wg"}) {
      fixHists[sSample+"_stacked"]=fixHists[sSample];
      TH1F &h=fixHists[sSample+"_stacked"];
      h.Add(&hStackSum);
      h.SetLineColor(h.GetFillColor());
      h.SetLineWidth(2);
      h.SetFillStyle(0);
   }

   TH1F hData;
   if (showData){
      hData=rebinned(*histReader.read<TH1F>(sSelection+sVar+"/SinglePhoton"));
      le.prepend(hData,"data","pe");
   }
   TH1F hStatErr(hStackSum);
   hStatErr.SetMarkerSize(0);
   TH1F hSystErr(hStatErr);
   float erri;
   float eVg,eGJ;
   // systematic uncert. from fit
   for (int i=0; i<=hSystErr.GetNbinsX()+1;i++){
      // total uncertainties:
      eVg=hVg.GetBinContent(i)*cfg.sf.e_Vg/cfg.sf.Vg;
      eGJ=hGJ.GetBinContent(i)*cfg.sf.e_GJ/cfg.sf.GJ;
      erri=util::quadSum<double>({eVg,eGJ});
      erri+=2*cfg.sf.rho*eVg*eGJ; // correlation
      hSystErr.SetBinError(i,TMath::Sqrt(erri));
   }
   // syst. uncert. from fixed backgrounds
   for (TString sSample:{"efake","TTcomb","diboson"}){
      TH1F const &h=fixHists[sSample];
      float const uncert=cfg.datasets.getSystUncert(sSample);
      for (int i=0; i<=hSystErr.GetNbinsX()+1;i++){
         erri=h.GetBinContent(i)*uncert;
         erri=util::sqrtQuadSum<double>({erri,hSystErr.GetBinError(i)});
         hSystErr.SetBinError(i,erri);
      }
   }

   // store separate background uncertainties
   if (plotMode==SR) {
      for (int i=0; i<=hSystErr.GetNbinsX()+1;i++){
         std::map<TString,SignalBin::Component> &bini=signalBins[i].component;
         bini["Vg"].count=hVg.GetBinContent(i);
         debug << i << hVg.GetBinContent(i);
         bini["Vg"].estat=hVg.GetBinError(i);
         bini["Vg"].esyst=hVg.GetBinContent(i)*cfg.sf.e_Vg/cfg.sf.Vg;

         bini["GJ"].count=hGJ.GetBinContent(i);
         bini["GJ"].estat=hGJ.GetBinError(i);
         bini["GJ"].esyst=hGJ.GetBinContent(i)*cfg.sf.e_GJ/cfg.sf.GJ;

         bini["data"].count=hData.GetBinContent(i);
         bini["data"].estat=bini["data"].esyst=0;

         for (TString sSample:{"efake","TTcomb","diboson"}){
            bini[sSample].count=fixHists[sSample].GetBinContent(i);
            bini[sSample].estat=fixHists[sSample].GetBinError(i);
            bini[sSample].esyst=fixHists[sSample].GetBinContent(i)*cfg.datasets.getSystUncert(sSample);
         }
         for (TString sSample:{"TChiWG","T5Wg"}) {
            signalBins_signal[i].component[sSample].count=fixHists[sSample].GetBinContent(i);
            signalBins_signal[i].component[sSample].estat=fixHists[sSample].GetBinError(i);
            signalBins_signal[i].component[sSample].esyst=0;
         }
      }
      double integr,integrErr;
      integr=hVg.IntegralAndError(2,hVg.GetNbinsX()+1,integrErr);
      bkgIntegral.component["Vg"].count=integr;
      bkgIntegral.component["Vg"].estat=integrErr;
      bkgIntegral.component["Vg"].esyst=integr*cfg.sf.e_Vg/cfg.sf.Vg;

      integr=hGJ.IntegralAndError(2,hGJ.GetNbinsX()+1,integrErr);
      bkgIntegral.component["GJ"].count=integr;
      bkgIntegral.component["GJ"].estat=integrErr;
      bkgIntegral.component["GJ"].esyst=integr*cfg.sf.e_GJ/cfg.sf.GJ;

      integr=hData.IntegralAndError(2,hData.GetNbinsX()+1,integrErr);
      bkgIntegral.component["data"].count=integr;
      bkgIntegral.component["data"].estat=bkgIntegral.component["data"].esyst=0;

      for (TString sSample:{"efake","TTcomb","diboson"}){
         integr=fixHists[sSample].IntegralAndError(2,fixHists[sSample].GetNbinsX()+1,integrErr);
         bkgIntegral.component[sSample].count=integr;
         bkgIntegral.component[sSample].estat=integrErr;
         bkgIntegral.component[sSample].esyst=integr*cfg.datasets.getSystUncert(sSample);
      }
      for (TString sSample:{"TChiWG","T5Wg"}) {
         integr=fixHists[sSample].IntegralAndError(2,fixHists[sSample].GetNbinsX()+1,integrErr);
         signalIntegral.component[sSample].count=integr;
         signalIntegral.component[sSample].estat=integrErr;
         signalIntegral.component[sSample].esyst=0;
      }
   }

   gfx::SplitCan spcan;
   TCanvas can;
   if (plotMode==PRE) {
      can.cd();
      can.SetLogy();
   } else {
      spcan.cdUp();
      spcan.pU_.SetLogy();
   }
   hStackSum.Draw("axis");

   if (plotMode==SR && sVar=="METS") hStackSum.GetYaxis()->SetRangeUser(0.2,300);
   if (plotMode==SR && sVar=="STg")  hStackSum.GetYaxis()->SetRangeUser(0.4,500);
   if (plotMode==CR && sVar=="METS") hStackSum.GetYaxis()->SetRangeUser(2e-2,200);
   if (plotMode==PRE && sVar=="METS_logx") can.SetLogx();
   st.Draw("same");
   st.GetYaxis()->SetTitle(hStatErr.GetYaxis()->GetTitle());
   st.GetXaxis()->SetTitle(hStatErr.GetXaxis()->GetTitle());
   if (plotMode==CR) logFile<<"-- "+saveName;

   for (TString sSample:{"TChiWG","T5Wg"}) {
      TH1F &h=fixHists[sSample+"_stacked"];
      if (showSignal) {
         if (plotMode==PRE) {
            TH1F &hr=fixHists[sSample];
            hr.SetLineColor(hr.GetFillColor());
            hr.SetLineWidth(2);
            hr.SetFillStyle(0);
            hr.Draw("hist same");
            hr.Draw("hist same");
            le.append(hr,sSample,"l");
         } else {
            h.Draw("hist same");
            le.append(h,sSample,"l");
         }
      }
      if (plotMode==CR) {
         float nS=fixHists[sSample].Integral(0,h.GetNbinsX()+1,"width");
         float nSB=h.Integral(0,h.GetNbinsX()+1,"width");
         logFile<<sSample+TString::Format(": s/s+b = %f/%f = %f",nS,nSB,nS/nSB);
      }
   }
   // hStatErr.SetFillStyle(3354);
   // hStatErr.SetFillColor(kBlack);
   // hStatErr.Draw("same e2");
    hSystErr.SetFillStyle(3354);
    hSystErr.SetFillColor(kGray);
    hSystErr.SetLineColor(kWhite);
    hSystErr.Draw("same e2");
    le.append(hSystErr,"#sigma_{syst}","f");
   if (showData) hData.Draw("same pe1");
   spcan.pU_.RedrawAxis();
   can.RedrawAxis();
   le.buildLegend(.55,.71,-1,-1,2).DrawClone();
   if (plotMode==CR) {
      TString txt="CR";
      if (sSelection.Contains("/0l")) txt+=", 0 leptons";
      if (sSelection.Contains("/1l")) txt+=", 1 lepton";
      if (sSelection.Contains("/0b")) txt+=", 0 b";
      gfx::cornerLabel(txt,1).DrawClone();
   }else if (plotMode==VR) {
      TString txt="VR";
      gfx::cornerLabel(txt,1).DrawClone();
   }
   if (plotMode==PRE) {
      saver.save(can,saveName,!showData);
   } else {
      spcan.cdLow();
      TGraphErrors grRatioStat=hist::getRatioGraph(hStatErr,st,"ratio",hist::ONLY1);
      TH1F hRatioSyst=hist::getRatio(hSystErr,st,"ratio",hist::ONLY1);
      TH1F hRatio;
      if (showData) {
         hRatio=hist::getRatio(hData,st,"ratio",hist::ONLY1);
         hRatio.Draw("axis e0");
         hRatio.SetMaximum(1.9);
         hRatio.SetMinimum(0.1);
      } else {
         hRatioSyst.Draw("axis e0");
         hRatioSyst.SetMaximum(1.9);
         hRatioSyst.SetMinimum(0.1);
      }
      hRatioSyst.SetFillStyle(1001);
      hRatioSyst.SetFillColor(kGray);
      hRatioSyst.SetLineColor(kGray);
      hRatioSyst.Draw("same e2");
      grRatioStat.SetFillStyle(1001);
      grRatioStat.SetFillColor(kGray+1);
      grRatioStat.SetLineColor(kGray+1);
      gfx::scaleXerrors(grRatioStat,.3);
      grRatioStat.Draw("2");
      hRatio.Draw("same e0 axig");
      if (showData) hRatio.Draw("same pe0");
      spcan.pL_.RedrawAxis();
      le.clear();
      le.append(grRatioStat,"#sigma_{stat}","f");
      le.append(hRatioSyst,"#sigma_{syst}","f");
      le.buildLegend(.4,.38,-1,.55,2).DrawClone();
      saver.save(spcan,saveName,!showData);
   }

   if (plotMode==SR && sVar=="STg" && !sSelection.Contains("STg600")) {
      bool tex=false;
      for (unsigned iBin=0; iBin<signalBins.size(); iBin++) {
         logFile*sVar/signalBins[iBin].lower>>signalBins[iBin].upper;
         logFile<<signalBins[iBin].formatBlock(tex,true);
         logFile<<signalBins_signal[iBin].formatBlock(tex);
      }
      logFile*sVar>>"integral";
      logFile<<bkgIntegral.formatBlock(tex,true);
      logFile<<signalIntegral.formatBlock(tex);

      io::Logger texSummary((sVar+"-summary.tex").Data());
      texSummary<<"\\begin{tabular}{LRRRR}";
      texSummary<<"\\ST [\\text{GeV}] &  \\text{prediction} & \\sigma_\\text{stat} &  \\sigma_\\text{syst} & \\text{data}\\\\\\hline";
      tex=true;
      unsigned nBins=signalBins.size();
      for (unsigned iBin=2; iBin<nBins-1; iBin++) {
         float l=signalBins[iBin].lower;
         float u=signalBins[iBin].upper;
         // actually contains overflow:
         if (iBin==nBins-2) u=std::numeric_limits<float>::infinity();

         io::Logger texFile((sVar+TString::Format("-%.0f-%.0f.tex",l,u)).Data());

         texFile<<"\\begin{tabular}{lrrrr}";
         TString binString=sVar+" bin: ";
         binString+=TString::Format("$%.0f-%.0f\\GeV$",l,u);
         binString.ReplaceAll("inf","\\infty");
         binString.ReplaceAll("STg","\\ST{}");
         texFile<<"   \\multicolumn{5}{c}{"+binString+"}\\\\";
         texFile<<"               &  yield & $\\sigma_\\text{stat}$ &  $\\sigma_\\text{syst}$ & $(/$tot. bkg.)\\\\\\hline";
         texFile<<signalBins[iBin].formatBlock(tex,true);
         texFile*signalBins_signal[iBin].formatBlock(tex);
         texFile<<"\\end{tabular}";

         binString=TString::Format("%.0f-%.0f     ",l,u);
         binString.ReplaceAll("inf","\\infty");
         binString.ReplaceAll("STg","\\ST{}");
         texSummary*TString::Format("%-16s &",binString.Data());
         SignalBin::Component tot=signalBins[iBin].getTotal();
    //     float esy=TMath::Sqrt(tot.esyst*tot.esyst+.13*.13*signalBins[iBin].component["GJ"].count*signalBins[iBin].component["GJ"].count);
         float esy=TMath::Sqrt(tot.esyst*tot.esyst);
         texSummary*TString::Format(" %6.1f & \\pm%6.1f & \\pm%6.1f &",tot.count,tot.estat,esy);
         texSummary<<TString::Format("%6.0f \\\\",signalBins[iBin].component["data"].count);
      }
      texSummary<<"\\end{tabular}";
      io::Logger texFile((sVar+"-integral.tex").Data());
      texFile<<"\\begin{tabular}{lrrrr}";
      texFile<<"               &  yield & $\\sigma_\\text{stat}$ &  $\\sigma_\\text{syst}$ & $(/$tot. bkg.)\\\\\\hline";
      texFile<<bkgIntegral.formatBlock(tex,true);
      texFile*signalIntegral.formatBlock(tex);
      texFile<<"\\end{tabular}";

      static io::RootFileSaver saver_yields("yields.root","");
      for (auto comp: signalBins[0].component) {
         TString sBkg=comp.first;
         TH1F h("","",nBins-2,0,1);
         TH1F hSys("","",nBins-2,0,1);
         for (unsigned iBin=1; iBin<nBins-1; iBin++) {
            h.SetBinContent(iBin,signalBins[iBin].component[sBkg].count);
            h.SetBinError(iBin,signalBins[iBin].component[sBkg].estat);
            hSys.SetBinContent(iBin,signalBins[iBin].component[sBkg].esyst);
         }
         saver_yields.save(h,sSelection+sVar+"/"+sBkg);
         saver_yields.save(hSys,sSelection+sVar+"/"+sBkg+"_esyst");
      }
   }
}

// plot functions to call
void plot(TString sSelection,TString sVar,int iRebin,PlotMode_t plotMode, int weight_index)
{
   return plot(sSelection,sVar,iRebin,{},{},plotMode,weight_index);
}
void plot(TString sSelection,TString sVar,
          std::vector<float> rebinEdges,std::vector<float> rebinWidths,
          PlotMode_t plotMode,
          int weight_index
   )
{
   return plot(sSelection,sVar,0,rebinEdges,rebinWidths,plotMode,weight_index);
}

extern "C"
void run()
{
   int check = 0; // 0 is SR STg, 1 is CR STg, 2 is CR photon pT
   
   if (check == 0){
   //Let 0 be nominal, 1 JESu, 2 JESd
      plot("pre_ph165/c_MET300/MT300/STg600/","STg",{600,800,1000,1300,1600},{200,200,300,300},SR,0);
      plot("pre_ph165/c_MET300/MT300/STg600/","STg_JESu",{600,800,1000,1300,1600},{200,200,300,300},SR,1);   
      plot("pre_ph165/c_MET300/MT300/STg600/","STg_JESd",{600,800,1000,1300,1600},{200,200,300,300},SR,2);
   }else if (check == 1){
   //check influence on CR depending on STg
      plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","STg",{600,800,1000,1300,1600},{200,200,300,300},SR,0);
      plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","STg_JESu",{600,800,1000,1300,1600},{200,200,300,300},SR,1);   
      plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","STg_JESd",{600,800,1000,1300,1600},{200,200,300,300},SR,2);
   }else if (check == 2){
   //check influence on CR depending on ph1Pt
      plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","ph1Pt",{0,500,750,1000,1300},{500,250,250,300},SR,0);
      plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","ph1Pt_JESu",{0,500,750,1000,1300},{500,250,250,300},SR,1);   
      plot("pre_ph165/c_MET100/MT100/METl300vMTl300/","ph1Pt_JESd",{0,500,750,1000,1300},{500,250,250,300},SR,2);
   }
}
