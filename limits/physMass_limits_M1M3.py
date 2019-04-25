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
			m3=hs.GetYaxis().GetBinCenter(j)
			
			data = pyslha.readSLHAFile("input/SLHA_M1_M3/out_M1_%i_M3_%i.slha"%(m1,m3))
			massblock = data.blocks["MASS"]
			if hs.GetBinContent(i,j)!=0:
				graph.SetPoint(k,massblock[1000022],massblock[1000021],hs.GetBinContent(i,j))
				quant.Fill(massblock[1000022],massblock[1000021])
				k+=1
				
				if (massblock[1000021]-massblock[1000022])<1:
					print massblock[1000022],m1,m3
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
	#~ c.Write("massDouble",rt.TObject.kOverwrite)
	c3=rt.TCanvas()
	c3.cd()
	quant.Draw("colz")
	#~ c3.Write("quant",rt.TObject.kOverwrite)
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

def smoothContour_knut(gr, neighbors=5, sigma=0.5):
    fgaus = rt.TF1("fgaus", "gaus", -10, 10)
    fgaus.SetParameters(1,0,1)
    weights = [fgaus.Eval(i*sigma) for i in range(neighbors)]
    #~ out = gr.Clone(aux.randomName())
    if gr: out = gr.Clone()
    else: return rt.TGraph()
    out.Set(0)
    n = gr.GetN()
    Xs = [gr.GetX()[i] for i in range(n)]
    Ys = [gr.GetY()[i] for i in range(n)]
    #~ n1 = Ys.index(max(Ys))+1
    #~ n2 = Ys.index(min(Ys))+1
    #~ nY=max(n1,n2)
    #~ n1 = Xs.index(max(Xs))+1
    #~ n2 = Xs.index(min(Xs))+1
    #~ nX=max(n1,n2)
    #~ n=max(nX,nY)
    Xs = Xs[0:n]
    Ys = Ys[0:n]
    for i, (x, y) in enumerate(zip(Xs,Ys)):
        pNeigh = min(neighbors, i+1, n-i)
        newX, ws, newY = 0, 0, 0
        for j in range(pNeigh):
            if j:
                newX += (Xs[i-j]+Xs[i+j])*weights[j]
                newY += (Ys[i-j]+Ys[i+j])*weights[j]
                ws += 2*weights[j]
            else:
                newX += x*weights[0]
                newY += y*weights[0]
                ws += weights[0]
        out.SetPoint(i, newX/ws, newY/ws)
    return out


###########################
###########Testing#########
###########################


#~ for cont in ["h_obs","h_exp","h_exp+1","h_exp-1","h_exp+2","h_exp-2"]:
#~ for analysis in ["inclusiv","htg","lepton","diphoton","allCombined_FullST"]:
#~ for analysis in ["inclusivNN","htgNN","lepton","diphoton","allCombined_highHtgNN"]:
#~ for analysis in ["inclusivFinal","htgFinal","lepton_final","diphoton_final","allCombined_final"]:
#~ for analysis in ["inclusivFinal","lepton_final","diphoton_final","allCombined_finalPre"]:
for analysis in ["allCombined_finalPre"]:
	print analysis
	#~ cont="h_exp"
	#~ cont="h_obs"
	cont="h_exp_xs"
	#~ cont="h_obs_xs"
	f=rt.TFile("output/limits_GGM_M1_M3_"+analysis+".root","read")
	hist=f.Get(cont)
	gr=toPhysMass(hist,cont)
	f.Close()

	#~ f=rt.TFile("output/physmass_GGM_M1_M3.root","update")
	#~ f=rt.TFile("output/physmass_GGM_M1_M3_NN.root","update")
	#~ f=rt.TFile("output/physmass_GGM_M1_M3_final.root","update")
	f=rt.TFile("output/physmass_GGM_M1_M3_finalPre.root","update")
	if (f.cd(analysis)!=1):
		f.mkdir(analysis)
	f.cd(analysis)
	
	x="M_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV)"
	y="m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV)"
	#~ hsInter=rt.TH2F("",";"+x+";"+y+";signal strength",100,75,750,100,2350,5300)
	#~ hsInter=rt.TH2F("",";"+x+";"+y+";signal strength",30,30,750,100,2350,5300)
	hsInter=rt.TH2F("",";"+x+";"+y+";",30,30,750,100,2350,5300)
	hsInter.GetYaxis().SetTitleOffset(1.3)
	for i in range(1,hsInter.GetXaxis().GetNbins()+1):
		for j in range(1,hsInter.GetYaxis().GetNbins()+1):
			hsInter.SetBinContent(i,j,gr.Interpolate(hsInter.GetXaxis().GetBinCenter(i),hsInter.GetYaxis().GetBinCenter(j)))
			if (hsInter.GetXaxis().GetBinCenter(i)>hsInter.GetYaxis().GetBinCenter(j)):
				hsInter.SetBinContent(i,j,0)
	
	if cont=="h_obs":
		hsInter.Write("hist_Inter_obs",rt.TObject.kOverwrite)
	elif cont=="h_obs_xs":
		hsInter.Write("hist_Inter_obs_xs",rt.TObject.kOverwrite)
		continue
	elif cont=="h_exp_xs":
		hsInter.Write("hist_Inter_exp_xs",rt.TObject.kOverwrite)
		continue
	else:
		hsInter.Write("hist_Inter",rt.TObject.kOverwrite)
    
	contour=getContourHS(hsInter)
	cont_sm=smoothContour_knut(contour)
	
	if cont=="h_obs":
		cont_sm.Write("cont_obs_sm",rt.TObject.kOverwrite)
	else:
		cont_sm.Write("cont_sm",rt.TObject.kOverwrite)

	f.Close()

"""
# Uncertainty bands for combined
#~ for analysis in ["allCombined_final"]:
for analysis in ["allCombined_finalPre"]:
	print analysis
	for cont in ["h_exp+1","h_exp-1","h_obs+1","h_obs-1"]:
		print cont
		f=rt.TFile("output/limits_GGM_M1_M3_"+analysis+".root","read")
		hist=f.Get(cont)
		gr=toPhysMass(hist,cont)
		f.Close()

		#~ f=rt.TFile("output/physmass_GGM_M1_M3.root","update")
		#~ f=rt.TFile("output/physmass_GGM_M1_M3_NN.root","update")
		#~ f=rt.TFile("output/physmass_GGM_M1_M3_final.root","update")
		f=rt.TFile("output/physmass_GGM_M1_M3_finalPre.root","update")
		f.mkdir(analysis)
		f.cd(analysis)
		
		x="M_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}} (GeV)"
		y="m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV)"
		hsInter=rt.TH2F("",";"+x+";"+y+";signal strength",30,30,750,100,2350,5300)
		hsInter.GetYaxis().SetTitleOffset(1.3)
		for i in range(1,hsInter.GetXaxis().GetNbins()+1):
			for j in range(1,hsInter.GetYaxis().GetNbins()+1):
				hsInter.SetBinContent(i,j,gr.Interpolate(hsInter.GetXaxis().GetBinCenter(i),hsInter.GetYaxis().GetBinCenter(j)))
				if (hsInter.GetXaxis().GetBinCenter(i)>hsInter.GetYaxis().GetBinCenter(j)):
					hsInter.SetBinContent(i,j,0)

		hsInter.Write("hist_Inter_"+cont,rt.TObject.kOverwrite)
	    
		contour=getContourHS(hsInter)
		cont_sm=smoothContour_knut(contour)
		
		cont_sm.Write("cont_"+cont+"_sm",rt.TObject.kOverwrite)
	    
		f.Close()
"""
