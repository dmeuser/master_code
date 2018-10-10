#ifndef PLOT_H
#define PLOT_H

/* Singleton-like color class for automatic colors*/

#include <vector>

#include "TColor.h"

class Color
{
public:
   static int next();
   static void reset();
private:
   std::vector<int> cols_={kGray,kGray+2,kAzure-8,kOrange-3,kGreen-6,kRed-6,kCyan-3,kMagenta-6,kGreen-5,kBlack};
   //std::vector<int> cols_={1,2,3,4,5,6,7,8,11,12};
   int size_=cols_.size();
   int pos_=-1;
   Color(){}
   static Color& instance_();
   int next_();
};

#endif
