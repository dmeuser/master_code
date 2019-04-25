#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

from finalPredictions import metHist

def getAcceptance(scanName, point, fullRange):
	
    binning = [350, 450, 600]
    
    dataset1 = Dataset(scanName, 0)
    hL = metHist(dataset1, point+"/original/nominal", binning+[800], True)
    hH = metHist(dataset1, point+"/original/nominal", binning+[800], False)
    if fullRange:
		acceptances = [hL.GetBinContent(i+1) for i in range(3)] + [hH.GetBinContent(i+1) for i in range(3)]
    else:
		acceptances = [hH.GetBinContent(i+1) for i in range(3)]
    return sum(acceptances)

def diffPrefire(scanName,dset,fullRange):
	dirs = aux.getDirNames(dset.files[0])
	
	graphAcc = ROOT.TGraph2D();
	
	for d in dirs:
		mass = d.split("_")
		if (getAcceptance(scanName,d,fullRange)!=0):
			graphAcc.SetPoint(graphAcc.GetN(), int(mass[0]), int(mass[1]), (getAcceptance(scanName,d,fullRange)-getAcceptance(scanName+"_prefire2016",d,fullRange))*100/(getAcceptance(scanName,d,fullRange)))
			if ((getAcceptance(scanName,d,fullRange)-getAcceptance(scanName+"_prefire2016",d,fullRange))*100/(getAcceptance(scanName,d,fullRange)))>1.7 :
				print int(mass[0]), int(mass[1])
		else:
			graphAcc.SetPoint(graphAcc.GetN(), int(mass[0]), int(mass[1]), 0)
	
	
	graphAcc.SetTitle("")
	style.style2d()
	histAcc = graphAcc.GetHistogram("old");
	#~ histAcc = graphAcc.GetHistogram();
	c = ROOT.TCanvas("","",600,500)
	
	if dset==ggm1:
		histAcc.GetXaxis().SetTitle("M_{1} (GeV)");
		histAcc.GetYaxis().SetTitle("M_{2} (GeV)");
	else:
		histAcc.GetXaxis().SetTitle("m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV)");
		histAcc.GetYaxis().SetTitle("m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
		histAcc.GetXaxis().SetRangeUser(1400,2100)
	histAcc.GetZaxis().SetTitle("relative efficience difference (%)");
	histAcc.SetStats(False);
	histAcc.Draw("colz")
      
	aux.Label2D(sim=True, status="Private work")
	
	if fullRange:
		c.SaveAs("plots/prefiring/"+scanName+".pdf")
		c.SetLogz();
		c.SaveAs("plots/prefiring/"+scanName+"_log.pdf")

	else:
		c.SaveAs("plots/prefiring/"+scanName+"_highHTG.pdf")
		c.SetLogz();
		c.SaveAs("plots/prefiring/"+scanName+"_highHTG_log.pdf")
	return 0
	
def nominalAcceptance(scanName,dset,fullRange):
	dirs = aux.getDirNames(dset.files[0])
	
	graphAcc = ROOT.TGraph2D();
	
	for d in dirs:
		mass = d.split("_")
		if (getAcceptance(scanName,d,fullRange)!=0):
			graphAcc.SetPoint(graphAcc.GetN(), int(mass[0]), int(mass[1]), getAcceptance(scanName,d,fullRange))
		else:
			graphAcc.SetPoint(graphAcc.GetN(), int(mass[0]), int(mass[1]), 0)
	
	
	graphAcc.SetTitle("")
	style.style2d()
	histAcc = graphAcc.GetHistogram("old");
	c = ROOT.TCanvas("","",600,500)
	
	if dset==ggm1:
		histAcc.GetXaxis().SetTitle("M_{1} (GeV)");
		histAcc.GetYaxis().SetTitle("M_{2} (GeV)");
		histAcc.SetMaximum(0.2)
		histAcc.SetMinimum(0.00001)
	else:
		histAcc.GetXaxis().SetTitle("m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV)");
		histAcc.GetYaxis().SetTitle("m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)");
		histAcc.SetMaximum(0.2)
	histAcc.GetZaxis().SetTitle("Acceptance");
	histAcc.SetStats(False);
	histAcc.Draw("colz")
      
	aux.Label2D(sim=True, status="Private work")
	
	if fullRange:
		c.SaveAs("plots/acceptance/"+scanName+".pdf")
		c.SetLogz();
		c.SaveAs("plots/acceptance/"+scanName+"_log.pdf")

	else:
		c.SaveAs("plots/acceptance/"+scanName+"_highHTG.pdf")
		c.SetLogz();
		c.SaveAs("plots/acceptance/"+scanName+"_highHTG_log.pdf")
	return 0


if __name__ == "__main__":
	diffPrefire("SMS-T5Wg",t5wg,True)
	#~ diffPrefire("SMS-T5Wg",t5wg,False)
	#~ diffPrefire("GGM_GravitinoLSP_M1-200to1500_M2-200to1500",ggm1,True)
	#~ diffPrefire("GGM_GravitinoLSP_M1-200to1500_M2-200to1500",ggm1,False)
	
	#~ nominalAcceptance("SMS-T5Wg",t5wg,True)
	#~ nominalAcceptance("SMS-T5Wg_prefire2016",t5wg,True)
	#~ nominalAcceptance("GGM_GravitinoLSP_M1-200to1500_M2-200to1500",ggm1,True)
	#~ nominalAcceptance("GGM_GravitinoLSP_M1-200to1500_M2-200to1500_prefire2016",ggm1,True)

