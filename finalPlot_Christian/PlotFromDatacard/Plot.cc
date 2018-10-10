#include "Plot.h"
#include "gfx.h"
#include "hist.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <vector>
#include <map>

#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "THStack.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"


enum ErrorType {
  ONLY1, ONLY2, COMB
};

class Background{
 public:
   Background(const std::string& n, double y, double u2, double s2, int c):name(n),yield(y),systematics2(u2),statistics2(s2),color(c){}
   std::string name;
   double yield, systematics2, statistics2;
   int color;
};

class Bin{
 public:
   std::string name, label;
   int data;
   double sig;
   std::vector<Background> bkgd;
   double Systematics() const {double r2=0;for (auto it=bkgd.begin();it!=bkgd.end();++it)r2+=it->systematics2;return sqrt(r2);}
   double Statistics()  const {double r2=0;for (auto it=bkgd.begin();it!=bkgd.end();++it)r2+=it->statistics2;return sqrt(r2);}
   double TotalBkgd(  ) const {double b=0; for (auto it=bkgd.begin();it!=bkgd.end();++it) b+=it->yield;return b;}
};


void SetColorMap(std::map<std::string,int>& c)
{
   Color::reset();
   c["sig"] = 2;
   c["signal"] = 2;
   c["SUSY"] = 2;
   c["t5gg"] = 2;

   c["rare"] = Color::next(); 
   c["efake"] = Color::next();
   c["qcdfakelep"] = Color::next();
   c["Vg"] = Color::next();
   c["gqcd"] =  Color::next();
   c["GJ"] = Color::next();
 
   c["wg"] =  c["Vg"];
   c["zg"] =  c["Vg"];
   c["ele"] = c["efake"];
   c["tg"] = c["rare"]; 
   c["jetfakepho"] = c["GJ"];
   c["elefakepho"] = c["efake"];
   c["VGamma"] = c["Vg"];
   c["ewk"] = c["Vg"];
   c["qcd"] = c["GJ"];
   c["zgg"] =c["rare"];
}

void GetSampleNames(std::map<int,std::string>& c)
{
   std::map<std::string,int> r;
   SetColorMap(r);
   c[1] = "Data";
   c[r["sig"]] = "Signal, GGM M_{1} = M_{2} = 1 TeV";
   c[r["rare"]] = "Rare backgrounds";
   c[r["efake"]] = "e #rightarrow #gamma mis-ID";
   c[r["qcdfakelep"]] = "jet#rightarrow lepton mis-ID";
   c[r["GJ"]] = "#gamma + jet";
   c[r["gqcd"]] = "QCD multi-jet";
   c[r["Vg"]] = "Vector-boson + #gamma";
}


void ParseLine(std::istringstream& iss, std::string& key, int &value)
{
      if (!(iss >> key >> value)) { throw; } // error
}

template<typename T> void ParseLine(std::istringstream& iss, std::string& key, std::vector<T>& value)
{
     if (!(iss >> key)) { throw; } // error
     T temp;
     while (iss >> temp) {value.push_back(temp); }
}

void ParseSystematicsLine(std::istringstream& iss, std::string& key, const std::vector<double>& rate, std::vector<double>& v2, std::vector<double>& s2)
{
   std::vector<double> temp;
   std::string str;
   double syst;
   unsigned index=0;
   if (!(iss >> key >> str)) { throw; } // error
   if (str=="gmN") {
      if (!(iss >> syst)) { throw; } // error
      while (iss >> str) { 
         std::stringstream ss(str);
         if (ss >> syst)
            s2[index]+=(syst-1.0)*(syst-1.0)*rate[index]*rate[index]; 
         //std::cout << key << "  " <<index << "   "<<str<< "   " << rate[index]<<"+="<<sqrt(v2[index])<<std::endl;
         ++index;
      }
   } else {
      while (iss >> str) { 
         std::stringstream ss(str);
         if (ss >> syst)
            v2[index]+=(syst-1.0)*(syst-1.0)*rate[index]*rate[index]; 
         //std::cout << key << "  " <<index << "   "<<str<< "   " << rate[index]<<"+="<<sqrt(v2[index])<<std::endl;
         ++index;      
      }   
   }  
}

template<typename T> void Print(std::ostream& os, std::vector<T>& v)
{
   for (auto it=v.begin(); it!=v.end(); ++it) os << *it << "  "; 
   os << std::endl;
}

void PrintTable(std::ostream& os, const std::vector<Bin>& bins)
{
   os << "";          for (auto it=bins.begin(); it!=bins.end(); ++it) os << "  &   " << it->name; os << std::endl;
   os << "Data  ";    for (auto it=bins.begin(); it!=bins.end(); ++it) os << "  &  $" << std::setw(10) << it->data << "$"; os << std::endl;
   os << "Signal";    for (auto it=bins.begin(); it!=bins.end(); ++it) os << "  &  $" << std::setw(10) << it->sig << "$";  os << std::endl;
   os << "Tot. Bkgd"; for (auto it=bins.begin(); it!=bins.end(); ++it) os << "  &  $" << std::setw(8)  << it->TotalBkgd()<<"\\pm"<<std::setw(8)<<it->Systematics() << "$";  os << std::endl;
}

TH1F getRatio(TH1F const &h1, TH1F const &h2, TString title, ErrorType et)
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


// static
Color& Color::instance_(){
   static Color instance;
   return instance;
}
// static
int Color::next(){
   return instance_().next_();
}
// static
void Color::reset(){
   instance_().pos_=-1;
}
int Color::next_(){
   pos_=(pos_+1)%size_;
   return cols_[pos_];
}



int main(int argc, char* argv[])
{
   
   ///Std Root Magic
   gStyle->SetHistFillColor(0);
   gStyle->SetPalette(1);
   //gStyle->SetFillColor(0);
   gStyle->SetOptStat(kFALSE);
   gStyle->SetCanvasColor(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetPadColor(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);

   gStyle->SetTitleFillColor(0);
   gStyle->SetTitleBorderSize(0);
   gStyle->SetTitleX(0.10);
   gStyle->SetTitleY(0.98);
   gStyle->SetTitleW(0.8);
   gStyle->SetTitleH(0.06);

   gStyle->SetErrorX(0.5);
   gStyle->SetStatColor(0);
   gStyle->SetStatBorderSize(0);
   gStyle->SetStatX(0);
   gStyle->SetStatY(0);
   gStyle->SetStatW(0);
   gStyle->SetStatH(0);

   gStyle->SetTitleFont(42);
   gStyle->SetTitleSize(0.048, "XYZ");
   gStyle->SetLabelFont(42,"X");
   gStyle->SetLabelFont(42,"Y");
   gStyle->SetLabelFont(42,"Z");
   gStyle->SetLabelSize(0.03,"X");
   gStyle->SetLabelSize(0.03,"Y");
   gStyle->SetLabelSize(0.03,"Z");

   ///Read input datacard file
   if (argc!=2) {
      std::cerr << "Error: Expected input datacard file as argument! Exiting..." <<std::endl;
      exit(1);
   }
   std::ifstream infile( argv[1] );
   std::string line;
   int linenumber=0;

   ///Variables to store input datacard file values
   int imax, jmax, kmax = 0;
   std::vector<std::string> bin_labels, systbin_labels, process_labels;
   std::vector<int> data, bin_IDs, process_IDs;
   std::vector<double> rates;

   ///Read input datacard line by line
   while (std::getline(infile, line))
   {
      ++linenumber;
      std::istringstream iss(line);
      if      (linenumber== 1) ParseLine(iss, line, imax);
      else if (linenumber== 2) ParseLine(iss, line, jmax);
      else if (linenumber== 3) ParseLine(iss, line, kmax);
      else if (linenumber== 5) ParseLine(iss, line, bin_labels);
      else if (linenumber== 6) ParseLine(iss, line, data);
      else if (linenumber== 8) ParseLine(iss, line, systbin_labels);
      else if (linenumber== 9) ParseLine(iss, line, process_labels);
      else if (linenumber==10) ParseLine(iss, line, process_IDs);
      else if (linenumber==11) ParseLine(iss, line, rates);
      else if (linenumber==12) break;
   }

   ///Print some of the read stuff
   std::cout<<imax<<"  "<<jmax<<"  "<<kmax<<std::endl;
   //Print(std::cout, bin_labels);
   //Print(std::cout, data);
   //Print(std::cout, process_labels);
   //Print(std::cout, rates);


   ///Read systematics from datacard
   std::vector<double> systematics2(rates.size(), 0.0), statistics2(rates.size(), 0.0);
   while (std::getline(infile, line))
   {
      std::istringstream iss(line);
      ParseSystematicsLine(iss, line, rates, systematics2, statistics2);
   }

   std::map<std::string,int> ColorMap;
   SetColorMap(ColorMap);

   ///Fill bin-wise vector
   std::vector<Bin> bins(imax);
   unsigned sample = 0;
   auto sID=process_IDs.begin(); 
   for (unsigned b=0; b<imax; ++b){
      bins[b].name = bin_labels[b];
      bins[b].label = std::to_string(b);
      bins[b].data = data[b];
      bins[b].sig = rates[sample];
      ++sample;
      ++sID;
      //bins[b].bkgd.clear();
      for (; *sID>0; ++sID,++sample) {
         bins[b].bkgd.push_back( Background( process_labels[sample],
                                             rates[sample],
                                             systematics2[sample],
                                             statistics2[sample],
                                             ColorMap[ process_labels[sample] ]
                                           ) ); 
      }
   }

   PrintTable(std::cout, bins);

   //TCanvas * c1 = new TCanvas("c1","c1", 1200, 600);
   //c1->SetLogy(1);
   gfx::SplitCan can;
   can.pU_.SetLogy();


   TH1F * a = new TH1F("template",";Signal bin;Events/bin",imax,+0.5,imax+0.5);
   TH1F * d = new TH1F("data",    "",imax,0.5,imax+0.5);
   TH1F * s = new TH1F("signal",  "",imax,0.5,imax+0.5);
   TH1F * tot = new TH1F("totalb",  "",imax,0.5,imax+0.5);
   TH1F * syst = new TH1F("syst",  "",imax,0.5,imax+0.5);
   TH1F * onlysyst = new TH1F("onlysyst",  "",imax,0.5,imax+0.5);
   TH1F * onlytotal = new TH1F("onlytotal",  "",imax,0.5,imax+0.5);
   TH1F * stat = new TH1F("stat",  "",imax,0.5,imax+0.5);
   TH1F * onlystat = new TH1F("onlystat",  "",imax,0.5,imax+0.5);
   THStack * st = new THStack("thst","");
   std::map<int,TH1F*> b;
   gfx::LegendEntries * leg = new gfx::LegendEntries();

   for (unsigned i=0;i<imax;++i) {
      d->SetBinContent(        i+1, bins[i].data);
      s->SetBinContent(        i+1, bins[i].sig);
      tot->SetBinContent(      i+1, bins[i].TotalBkgd());
      tot->SetBinError(        i+1, sqrt(pow(bins[i].Systematics(),2)+pow(bins[i].Statistics(),2)));
      syst->SetBinContent(     i+1, bins[i].TotalBkgd());
      syst->SetBinError(       i+1, sqrt(pow(bins[i].Systematics(),2)+pow(bins[i].Statistics(),2)));
      onlytotal->SetBinContent(i+1, sqrt(pow(bins[i].Systematics(),2)+pow(bins[i].Statistics(),2)));
      onlytotal->SetBinError(  i+1, sqrt(pow(bins[i].Systematics(),2)+pow(bins[i].Statistics(),2)));
      //~ onlysyst->SetBinContent( i+1, bins[i].Systematics());
      onlysyst->SetBinContent( i+1, bins[i].TotalBkgd());
      onlysyst->SetBinError(   i+1, bins[i].Systematics());
      stat->SetBinContent(     i+1, bins[i].TotalBkgd());
      stat->SetBinError(       i+1, bins[i].Statistics());
      onlystat->SetBinContent( i+1, bins[i].Statistics());
      onlystat->SetBinError(   i+1, bins[i].Statistics());
      for (auto bit = bins[i].bkgd.begin(); bit!=bins[i].bkgd.end(); ++bit){
         if (b.find(bit->color) == b.end()) { 
            b[bit->color] = new TH1F(("b"+bit->name).c_str(),  "",imax,0.5,imax+0.5);
            b[bit->color]->SetFillColor(bit->color);
            b[bit->color]->SetLineColor(bit->color+1);
         }   
         b[bit->color]->SetBinContent(i+1, b[bit->color]->GetBinContent(i+1) + bit->yield );
         //b[bit->color]->SetBinError(i+1, sqrt(b[bit->color]->GetBinError(i+1)*b[bit->color]->GetBinError(i+1) + bit->systematics2 + bit->statistics2) );
         //if (i==0) std::cout << i << "  "<< bit->name<< "  c:"<<bit->color <<",  yield:"<<b[bit->color]->GetBinContent(i+1) <<", tot: "<< bins[i].TotalBkgd()<<std::endl;
      }   
   }  
   TGraphErrors h_syst(syst); 
   h_syst.SetFillStyle(3354);
   h_syst.SetFillColor(kRed+1);
   h_syst.SetLineColor(kWhite);
 
   std::map<int,std::string> BNames;
   GetSampleNames(BNames);
   for (auto it=b.rbegin(); it!=b.rend(); ++it) {
         st->Add( it->second );
         for (int ibs=0; ibs<bins.size(); ++ibs){
            std::vector<Background>::iterator this_b = bins[ibs].bkgd.begin();
            for (;this_b->color!=it->first && this_b!=bins[ibs].bkgd.end(); ++this_b){}
            if (this_b==bins[ibs].bkgd.end()) continue;
            leg->prepend( *(it->second), BNames[it->first].c_str(),"f");
            break;
         }
   }      
   leg->append(*s,"GGM M_{1} = 1000 GeV M_{2} = 750 GeV","l");
   leg->prepend(h_syst,"Total bkg. unc.","fe");
   leg->prepend(*d,"Data","pe");

   a->SetMinimum(0.5);
   a->SetMaximum( 60*d->GetMaximum() + sqrt(d->GetMaximum()));
   d->SetLineColor( kBlack ); 
   d->SetMarkerStyle( 8 ); 
   d->SetMarkerSize( 0.5 ); 
   s->SetLineColor( kRed );
   s->SetLineWidth( 2 );
   
   a->Draw();
   st->Draw("hf,same");
   h_syst.Draw("2");
   s->Draw("h,same");
   d->Draw("pe0,same");
   tot->Draw("h,same");
   TLegend * dleg = new TLegend(leg->buildLegend(0.4,0.7,0.9,0.9));
   dleg->SetNColumns(3);
   dleg->Draw("same");

   TLine * aline = new TLine();
   TLatex * atext = new TLatex();
   atext->SetTextSize(0.025);
   atext->SetTextFont();
   aline->SetLineWidth(2);
   aline->DrawLine(4.5,a->GetMinimum(),4.5,1.5*syst->GetMaximum());
   aline->DrawLine(7.5,a->GetMinimum(),7.5,1.5*syst->GetMaximum());
   aline->DrawLine(43.5,a->GetMinimum(),43.5,1.5*syst->GetMaximum());
   atext->DrawLatex(2.0,1.1*syst->GetMaximum(),"S_{#scale[.8]{T}}^{#scale[.8]{#gamma}}");
   atext->DrawLatex(5.5,1.1*syst->GetMaximum(),"H_{#scale[.8]{T}}^{#scale[.8]{#gamma}}");
   atext->DrawLatex(13.5,1.1*syst->GetMaximum(),"Lepton");
   atext->DrawLatex(44.5,1.1*syst->GetMaximum(),"Diphoton");

   gPad->RedrawAxis();

   can.cdLow();

   TLegend * sleg = new TLegend(0.65,0.89,0.85,0.99);


   TH1F hRatio=hist::getRatio(*d,*tot,"Data/Pred.",hist::ONLY1);
   TGraphErrors grRatioStat=hist::getRatioGraph(*stat,*tot,"Ratio",hist::ONLY1);
   TGraphErrors grRatioSyst=hist::getRatioGraph(*onlysyst,*tot,"Ratio",hist::ONLY1);
   TGraphErrors grRatioTot=hist::getRatioGraph(*tot,*tot,"Ratio",hist::ONLY1);
   hRatio.GetXaxis()->SetTitle("Search bin");
   hRatio.GetXaxis()->SetTitleFont(42);   
   //~ hRatio.SetMaximum(2.7);
   hRatio.SetMaximum(4.2);
   hRatio.SetMinimum(0.);
   hRatio.Draw("axis");
   //~ grRatioStat.SetFillStyle(1001);
   //~ grRatioStat.SetFillColor(kGray+1);
   //~ grRatioStat.SetLineColor(kGray+1);
   grRatioTot.SetFillStyle(1001);
   grRatioTot.SetFillColor(kGray+1);
   grRatioTot.SetLineColor(kGray+1);
   grRatioSyst.SetFillStyle(1001);
   grRatioSyst.SetFillColor(kGray);
   grRatioSyst.SetLineColor(kGray);
   //~ gfx::scaleXerrors(grRatioTot,.3);
           
   
   //~ grRatioSyst.Draw("2");
   grRatioTot.Draw("2");     
   hRatio.Draw("pe,same");
   aline->SetLineWidth(1);   
   aline->DrawLine(0,1,49.5,1);  
   sleg->AddEntry(&grRatioTot,"Total bkg. unc.","f");
   //~ sleg->AddEntry(&grRatioSyst,"Syst.","f");

   
/*

   TH1F hPull=hist::getRatio(hist::getPull(*d,*tot,"Pull (d-b)/#sigma_{tot}",hist::ONLY1), *onlytotal, "Pull (d-b)/#sigma_{tot}",hist::ONLY1);
   TGraphErrors grPullStat=hist::getPullGraph(*onlystat,*onlytotal,"Ratio",hist::ONLY1);
   TGraphErrors grPullSyst=hist::getPullGraph(*onlysyst,*onlytotal,"Ratio",hist::ONLY1);
   grPullStat.SetFillStyle(3345);
   grPullStat.SetFillColor(kMagenta+4);
   grPullSyst.SetFillStyle(3354);
   grPullSyst.SetFillColor(kRed+1);
   gfx::scaleXerrors(grPullStat,.99);
   gfx::scaleXerrors(grPullSyst,.85);

   grPullStat.Draw("2");
   grPullSyst.Draw("2");

   hPull.GetXaxis()->SetTitle("Signal region");
   hPull.GetXaxis()->SetTitleFont(42);  
   hPull.SetMarkerStyle(8); 
   hPull.SetMaximum(2.95);
   hPull.SetMinimum(-2.95);
   hPull.Draw("axis");
   hPull.Draw("pe,same");
   grPullStat.Draw("2");
   aline->SetLineWidth(1);   
   aline->DrawLine(0,0,49,0);
   sleg->AddEntry(&grPullStat,"Stat.","f");
   sleg->AddEntry(&grPullSyst,"Syst.","f");
*/

   sleg->SetNColumns(2); 
   sleg->SetBorderSize(0);
   sleg->SetFillStyle(0);
   sleg->Draw("same");

   gPad->RedrawAxis();
   setupDrawnAxes(can);


   can.draw_lumi();
   can.SaveAs("plot.pdf");
   can.SaveAs("plot.root");




   return 0;
}


