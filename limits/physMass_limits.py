#Script to create limit for the physical masses in the GGM scan

import os, re, subprocess as sp, math
from distutils import spawn
import ROOT as rt
import pyslha


def toPhysMass(hs,cont):
	graph=rt.TGraph2D()
	quant=rt.TH2F("",";M_Neutralino (GeV);M_Chargino (GeV);Number of signal points",100,100,725,100,135,1240)
	double=rt.TGraph()
	double.SetTitle(";M_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV);signal strength")
	k=0
	l=0
	for i in range(1,hs.GetNbinsX()+1):
		for j in range(1,hs.GetNbinsY()+1):
			
			m1=hs.GetXaxis().GetBinCenter(i)
			m2=hs.GetXaxis().GetBinCenter(j)
			
			data = pyslha.readSLHAFile("input/SLHA_M1_M2/out_M1_%iM2_%i.slha"%(m1,m2))
			massblock = data.blocks["MASS"]
			if hs.GetBinContent(i,j)!=0 and m1>200 and m2>200:
				#~ print m1,m2
				graph.SetPoint(k,massblock[1000022],massblock[1000024],hs.GetBinContent(i,j))
				quant.Fill(massblock[1000022],massblock[1000024])
				k+=1
				
				if (massblock[1000024]-massblock[1000022])<1:
					print massblock[1000022],m1,m2
					double.SetPoint(l,massblock[1000022],hs.GetBinContent(i,j))
					l+=1
	f=rt.TFile("test/physmass.root","update")
	if (f.cd(cont)!=1):
		f.mkdir(cont)
	f.cd(cont)
	c=rt.TCanvas()
	c.cd()
	double.SetMarkerStyle(20)
	double.SetMarkerSize(0.5)
	double.Draw("AP1")
	c.Write("massDouble",rt.TObject.kOverwrite)
	c3=rt.TCanvas()
	c3.cd()
	quant.Draw("colz")
	c3.Write("quant",rt.TObject.kOverwrite)
	f.Close()
	return graph


def getContour(gr2):
	""" inspired by https://hypernews.cern.ch/HyperNews/CMS/get/susy/2122/2.html
	"""
	c=rt.TCanvas()
	gr2.Draw("tri1")
	c.Update()
	cont=gr2.GetContourList(1.)
	print cont
	if cont and len(cont)>0:
		N=0
		c=None
		for ci in cont:
			if ci.GetN() > N:
				N=ci.GetN()
				c=ci
		return c
	else:
		print "GetContourList failed"
		return None

def getContourHS(hist,flag=False):
	
	print hist
	if flag:
		for i in xrange(hist.GetNbinsX()+1):
			for j in xrange(hist.GetNbinsY()+1):
				if hist.GetBinContent(i,j)<1 :
					hist.SetBinContent(i,j,0)
				else :
					hist.SetBinContent(i,j,1.5)
	
	c=rt.TCanvas()
	c.cd()
	hist.SetContour(1)
	hist.SetContourLevel(0,1.)
	hist.Draw("CONT Z LIST")
	c.Update()
	conts=rt.gROOT.GetListOfSpecials().FindObject("contours")
	cont=conts.At(0)
	
	i=0
	if cont and len(cont)>0:
		N=0
		c=None
		for ci in cont:
			ci.Write("cont"+str(i),rt.TObject.kOverwrite)
			i+=1
			if ci.GetN() > N:
				N=ci.GetN()
				c=ci
		return c
	else:
		print "GetContourList failed"
		return None
		
		
#~ def getContours():
	#~ f=rt.TFile(outdir+"limits_%s_"%sScan+selection+".root","update")
	#~ for lvl in ["obs","obs+1","obs-1","exp","exp+1","exp-1","exp+2","exp-2"]:
		#~ 
		#~ hs=f.Get("h_"+lvl)
		#~ grPhys=toPhysMass(hs)
		#~ grC=getContour(grPhys)
		#~ 
		#~ if grC:
			#~ if(f.cd("physmass")!=1):
				#~ f.mkdir("physmass")
			#~ f.cd("physmass")
			#~ grC.Write("gr_"+lvl+"C",rt.TObject.kOverwrite)
		#~ else:
			#~ print "could not get contour for",lvl
	#~ f.Close()
#~ 
#~ 
#~ def smoothContours():
	#~ f=rt.TFile(outdir+"limits_%s_"%sScan+selection+".root","update")
	#~ gs=rt.TGraphSmooth()
	#~ f.cd("physmass")
	#~ for lvl in ["obs","obs+1","obs-1","exp","exp+1","exp-1","exp+2","exp-2"]:
		#~ gr=f.Get("physmass/gr_"+lvl+"C")
		#~ print lvl
        #~ gr_sm=gs.Approx(gr)
        #~ if gr_sm:
			#~ print "okay"
        #~ # gr_sm=gs.Approx(gr)
        #~ gr_sm=gs.SmoothSuper(gr)
        #~ gr_sm.Write("gr_"+lvl+"C_sm",rt.TObject.kOverwrite)
	#~ f.Close()
#~ 
#~ 
#~ outdir="output/"
#~ sScan="GGM_M1_M2"
#~ selection="inclusiv"
#~ getContours()
#~ smoothContours()


###########################
###########Testing#########
###########################

#~ for cont in ["h_obs","h_exp","h_exp+1","h_exp-1","h_exp+2","h_exp-2"]:
for cont in ["h_exp"]:
	f=rt.TFile("output/limits_GGM_M1_M2_inclusiv.root","read")
	hist=f.Get(cont)
	gr=toPhysMass(hist,cont)
	f.Close()

	f=rt.TFile("test/physmass.root","update")
	f.cd(cont)
	
	x="M_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV)"
	y="M_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV)"
	hsInter=rt.TH2F("",";"+x+";"+y+";signal strength",100,100,725,100,135,1240)
	hsInter.GetYaxis().SetTitleOffset(1.3)
	for i in range(1,hsInter.GetXaxis().GetNbins()+1):
		for j in range(1,hsInter.GetYaxis().GetNbins()+1):
			hsInter.SetBinContent(i,j,gr.Interpolate(hsInter.GetXaxis().GetBinCenter(i),hsInter.GetYaxis().GetBinCenter(j)))
			if (hsInter.GetXaxis().GetBinCenter(i)>hsInter.GetYaxis().GetBinCenter(j)):
				hsInter.SetBinContent(i,j,0)
	hsInter.Write("hist_Inter",rt.TObject.kOverwrite)
	
	#~ if cont=="h_obs":
		#~ for i in range(100,725):
			#~ hsInter.SetBinContent(hsInter.GetXaxis().FindBin(i),hsInter.GetYaxis().FindBin(i),2)
		#~ hsInter.Write("hist_Inter_SetDia",rt.TObject.kOverwrite)
	
	contHSinter=getContourHS(hsInter)
	contHSinter.Write("contHSinter",rt.TObject.kOverwrite)

	f.Close()

