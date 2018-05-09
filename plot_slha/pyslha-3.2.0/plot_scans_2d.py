#script to create 2d plots of scans given by slha files
#Use CMSSW810 to be able to import matplotlib

import pyslha
import numpy as np
import os,re
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
import scipy.interpolate as inter
import ROOT as rt
import scipy.constants as const

def plot_2d_M1M2(grid,label,title,output):
	M2,M1 = np.mgrid[25:1525:31j,25:1525:31j]
	plt.figure()
	if output == "Width NLSP" or output == "Length NLSP" or output == "NLO_ms_prospino" or output=="MassDiff":
		plt.pcolor(M1, M2, grid, cmap=plt.get_cmap('jet'), norm=LogNorm(vmin=grid.min(), vmax=grid.max()),linewidth=0, rasterized=True)
	else:
		plt.pcolor(M1, M2, grid, cmap=plt.get_cmap('jet'), vmin=grid.min(), vmax=grid.max(),linewidth=0, rasterized=True)
	if label.find("TeV") == -1:
		cbar = plt.colorbar()
	else : cbar = plt.colorbar(format="%.2f")
	cbar.set_label(label,fontsize = 20)
	plt.title(title)
	plt.xlim(25,1525)
	plt.ylim(25,1525)
	plt.xlabel("$M_1$ (GeV)",fontsize = 20)
	plt.ylabel("$M_2$ (GeV)",fontsize = 20)
	plt.savefig("output/M1_M2/"+output+".pdf",format="pdf")

def plot_2d_M1M3(grid,label,title,output,mctest):
	M3,M1 = np.mgrid[975:2525:32j,25:1525:31j]
	if mctest == True : 
		M3,M1 = np.mgrid[975:2475:31j,25:1525:31j]
	plt.figure()
	if output == "NLO_ms_prospino":
		plt.pcolor(M1, M3, grid, cmap=plt.get_cmap('jet'), norm=LogNorm(vmin=1E-5, vmax=grid.max()),linewidth=0, rasterized=True)
	else:
		plt.pcolor(M1, M3, grid, cmap=plt.get_cmap('jet'), vmin=grid.min(), vmax=grid.max(),linewidth=0, rasterized=True)
	if label.find("TeV") == -1:
		cbar = plt.colorbar()
	else : cbar = plt.colorbar(format="%.2f")
	cbar.set_label(label,fontsize = 20)
	plt.title(title)
	plt.xlim(25,1525)
	plt.ylim(975,2525)
	plt.xlabel("$M_1$ (GeV)",fontsize = 20)
	plt.ylabel("$M_3$ (GeV)",fontsize = 20)
	plt.savefig("output/M1_M3/"+output+".pdf",format="pdf")

def run_M1M2():
	readout = []

	for f in os.listdir("../../susyhit/output/Bench_M1_M2/merged/"):
		if f[-5:] == ".slha":
			data = pyslha.readSLHAFile("../../susyhit/output/Bench_M1_M2/merged/"+f)
			massblock = data.blocks["MASS"]
			extparblock = data.blocks["EXTPAR"]
			uchmixblock = data.blocks["UMIX"]
			vchmixblock = data.blocks["VMIX"]
			neutrmixblock = data.blocks["NMIX"]
			decays_nlsp = data.decays[1000022]
			width = decays_nlsp.totalwidth
			print f
			#if width>1e-20: width = 1
			if len(decays_nlsp.decays)==1:
				readout.append([extparblock[1],extparblock[2],massblock[1000022],massblock[1000024],massblock[1000021],massblock[1000006],massblock[35],decays_nlsp.decays[0].br,0.,width])
			else:
				if decays_nlsp.decays[0].ids[1]==22:
					brg = decays_nlsp.decays[0].br
					brz = decays_nlsp.decays[1].br
				else :
					brz = decays_nlsp.decays[0].br
					brg = decays_nlsp.decays[1].br
				readout.append([extparblock[1],extparblock[2],massblock[1000022],massblock[1000024],massblock[1000021],massblock[1000006],massblock[35],brg,brz,width])
	def getKey(item):
		return (item[0],item[1])
		
	readout = sorted(readout,key=getKey)
	
	M1 = np.array([item[0] for item in readout])
	M2 = np.array([item[1] for item in readout])
	M_chi_10 = np.array([item[2] for item in readout])
	M_chi_1p = np.array([item[3] for item in readout])
	M_gl = np.array([item[4] for item in readout])
	M_st = np.array([item[5] for item in readout])
	M_H = np.array([item[6] for item in readout])
	BR_g = np.array([item[7] for item in readout])
	BR_Z = np.array([item[8] for item in readout])
	width_nlsp = np.array([item[9] for item in readout])
	diffMass = M_chi_1p-M_chi_10
	decay_length = const.hbar/(const.e*width_nlsp)*const.c*1E-7
	
	values = [M_chi_10,M_chi_1p,M_gl,M_st,M_H,BR_g,BR_Z,width_nlsp,decay_length,diffMass]
	names = ["massNLSP","massChargino","massGluino","massStop","massHiggsH","BRtoPhoton","BRtoZ","widthNLSP","length_nlsp","massDiff"]
	j = 0
	rfile = rt.TFile("output/M1_M2/hist2d.root","UPDATE")
	
	
	for temp in values:
		name = names[j]
		histo = rt.TH2F(name,";M1 (GeV);M2 (GeV);"+name,30,25,1525,30,25,1525)
		for i in xrange(len(temp)):
			histo.Fill(M1[i],M2[i],temp[i])
		
		histo.GetYaxis().SetTitleOffset(1.3);
		histo.Write(name,rt.TObject.kOverwrite)
		j +=1
		
	rfile.Close()

	M_chi_10 = (M_chi_10.reshape((30,30))).T
	M_chi_1p = (M_chi_1p.reshape((30,30))).T
	M_gl = (M_gl.reshape((30,30))).T
	M_st = (M_st.reshape((30,30))).T
	M_H = (M_H.reshape((30,30))).T
	BR_g = (BR_g.reshape((30,30))).T
	BR_Z = (BR_Z.reshape((30,30))).T
	width_nlsp = (width_nlsp.reshape((30,30))).T
	diffMass = (diffMass.reshape((30,30))).T
	decay_length = (decay_length.reshape((30,30))).T
	
	
	decay_length_thresh = const.hbar/(const.e*width_nlsp)*const.c*1E-7
	for i in xrange(30):
		for j in xrange(30):
			if decay_length_thresh[i,j] > 2:
				decay_length_thresh[i,j] = 10

	plot_2d_M1M2(M_chi_10,"$M_{\\tilde{\\chi}^0_1}$ (GeV)","Mass $\\tilde{\\chi}^0_1$","Mass_neutralino")
	plot_2d_M1M2(M_chi_1p,"$M_{\\tilde{\\chi}^{\\pm}_1}$ (GeV)","Mass $\\tilde{\\chi}^{\\pm}_1$","Mass_chargino")
	plot_2d_M1M2(M_gl,"$M_{\\tilde{g}}$ (GeV)","Mass $\\tilde{g}$","Mass_gluino")
	plot_2d_M1M2(M_st,"$M_{\\tilde{t}_1}$ (GeV)","Mass $\\tilde{t}_1$","Mass_stop")
	plot_2d_M1M2(M_H,"$M_{H}$ (GeV)","Mass $H$","Mass_Higgs")
	plot_2d_M1M2(BR_g,"BR($\\tilde{\\chi}^0_1\\rightarrow\\gamma$)","Branching Ratio $\\tilde{\\chi}^0_1\\rightarrow\\gamma$","BR_gamma")
	plot_2d_M1M2(BR_Z,"BR($\\tilde{\\chi}^0_1\\rightarrow Z$)","Branching Ratio $\\tilde{\\chi}^0_1\\rightarrow Z$","BR_z")
	plot_2d_M1M2(width_nlsp,"Decay width $\\tilde{\\chi}^0_1$ (GeV)","Decay width $\\tilde{\\chi}^0_1$","Width NLSP")
	plot_2d_M1M2(decay_length,"Decay length $\\tilde{\\chi}^0_1$ (cm)","Decay length $\\tilde{\\chi}^0_1$","Length NLSP")
	plot_2d_M1M2(decay_length_thresh,"Decay length $\\tilde{\\chi}^0_1$ (cm)","Decay length tresh $\\tilde{\\chi}^0_1$","Length tresh NLSP")
	plot_2d_M1M2(diffMass,"$\\Delta M(\\tilde{\\chi}^{\\pm}_1,\\tilde{\\chi}^0_1$) (GeV)","$\\Delta M(\\tilde{\\chi}^{\\pm}_1,\\tilde{\\chi}^0_1$)","MassDiff")
	
def run_M1M3():
	readout = []

	for f in os.listdir("../../susyhit/output/Bench_M1_M3/merged/"):
		if f[-5:] == ".slha":
			data = pyslha.readSLHAFile("../../susyhit/output/Bench_M1_M3/merged/"+f)
			massblock = data.blocks["MASS"]
			extparblock = data.blocks["EXTPAR"]
			uchmixblock = data.blocks["UMIX"]
			vchmixblock = data.blocks["VMIX"]
			neutrmixblock = data.blocks["NMIX"]
			decays_nlsp = data.decays[1000022]
			print f
			if len(decays_nlsp.decays)==1:
				readout.append([extparblock[1],extparblock[3],massblock[1000022],massblock[1000024],massblock[1000021],massblock[1000006],massblock[35],decays_nlsp.decays[0].br,0.])
			else:
				if decays_nlsp.decays[0].ids[1]==22:
					brg = decays_nlsp.decays[0].br
					brz = decays_nlsp.decays[1].br
				elif decays_nlsp.decays[0].ids[1]==23:
					brz = decays_nlsp.decays[0].br
					brg = decays_nlsp.decays[1].br
				readout.append([extparblock[1],extparblock[3],massblock[1000022],massblock[1000024],massblock[1000021],massblock[1000006],massblock[35],brg,brz])
	def getKey(item):
		return (item[0],item[1])
		
	readout = sorted(readout,key=getKey)
	
	M1 = np.array([item[0] for item in readout])
	M3 = np.array([item[1] for item in readout])
	M_chi_10 = np.array([item[2] for item in readout])
	M_chi_1p = np.array([item[3] for item in readout])
	M_gl = np.array([item[4] for item in readout])
	M_st = np.array([item[5] for item in readout])
	M_H = np.array([item[6] for item in readout])
	BR_g = np.array([item[7] for item in readout])
	BR_Z = np.array([item[8] for item in readout])
	
	values = [M_chi_10,M_chi_1p,M_gl,M_st,M_H,BR_g,BR_Z]
	names = ["massNLSP","massChargino","massGluino","massStop","massHiggsH","BRtoPhoton","BRtoZ"]
	j = 0
	rfile = rt.TFile("output/M1_M3/hist2d_scan2.root","UPDATE")
	
	
	for temp in values:
		name = names[j]
		histo = rt.TH2F(name,";M1 (GeV);M3 (GeV);"+name,30,25,1525,31,975,2525)
		for i in xrange(len(temp)):
			histo.Fill(M1[i],M3[i],temp[i])
		
		histo.GetYaxis().SetTitleOffset(1.3);
		histo.Write(name,rt.TObject.kOverwrite)
		j +=1
		
	rfile.Close()
	
	M_chi_10 = (M_chi_10.reshape((30,31))).T
	M_chi_1p = (M_chi_1p.reshape((30,31))).T
	M_gl = (M_gl.reshape((30,31))).T
	M_st = (M_st.reshape((30,31))).T
	M_H = (M_H.reshape((30,31))).T
	BR_g = (BR_g.reshape((30,31))).T
	BR_Z = (BR_Z.reshape((30,31))).T
	
	plot_2d_M1M3(M_chi_10,"$M_{\\tilde{\\chi}^0_1}$ (GeV)","Mass $\\tilde{\\chi}^0_1$","Mass_neutralino",False)
	plot_2d_M1M3(M_chi_1p,"$M_{\\tilde{\\chi}^{\\pm}_1}$ (GeV)","Mass $\\tilde{\\chi}^{\\pm}_1$","Mass_chargino",False)
	plot_2d_M1M3(M_gl,"$M_{\\tilde{g}}$ (GeV)","Mass $\\tilde{g}$","Mass_gluino",False)
	plot_2d_M1M3(M_st,"$M_{\\tilde{t}_1}$ (GeV)","Mass $\\tilde{t}_1$","Mass_stop",False)
	plot_2d_M1M3(M_H,"$M_{H}$ (GeV)","Mass $H$","Mass_Higgs",False)
	plot_2d_M1M3(BR_g,"BR($\\tilde{\\chi}^0_1\\rightarrow\\gamma$)","Branching Ratio $\\tilde{\\chi}^0_1\\rightarrow\\gamma$","BR_gamma",False)
	plot_2d_M1M3(BR_Z,"BR($\\tilde{\\chi}^0_1\\rightarrow Z$)","Branching Ratio $\\tilde{\\chi}^0_1\\rightarrow Z$","BR_z",False)

def M1M2_BRphys():
	readout = []

	for f in os.listdir("../../susyhit/output/Bench_M1_M2/merged/"):
		if f[-5:] == ".slha":
			print f
			data = pyslha.readSLHAFile("../../susyhit/output/Bench_M1_M2/merged/"+f)
			massblock = data.blocks["MASS"]
			extparblock = data.blocks["EXTPAR"]
			decays_nlsp = data.decays[1000022]
			if len(decays_nlsp.decays)==1:
				readout.append([massblock[1000022],massblock[1000024],decays_nlsp.decays[0].br,0.,extparblock[1],extparblock[2]])
			else:
				if decays_nlsp.decays[0].ids[1]==22:
					brg = decays_nlsp.decays[0].br
					brz = decays_nlsp.decays[1].br
				else :
					brz = decays_nlsp.decays[0].br
					brg = decays_nlsp.decays[1].br
				readout.append([massblock[1000022],massblock[1000024],brg,brz,extparblock[1],extparblock[2]])
		
		#~ if (len(readout)>40): break

	def getKey(item):
		return (item[0],item[1])
		
	readout = sorted(readout,key=getKey)
	
	M_chi1 = np.array([item[0] for item in readout])
	M_char1 = np.array([item[1] for item in readout])
	BRg = np.array([item[2] for item in readout])
	BRz = np.array([item[3] for item in readout])
	M1 = np.array([item[4] for item in readout])
	M2 = np.array([item[5] for item in readout])
	
	Graph = rt.TGraph2D()
	for i in xrange(len(BRg)):
		Graph.SetPoint(i,M_chi1[i],M_char1[i],BRg[i])
	
	
	cd = rt.TCanvas()
	chi1 = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
	char1 = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
	Graph.GetHistogram().SetTitle(";m#kern[0.1]{_{"+chi1+"}} (GeV);m#kern[0.1]{_{"+char1+"}} (GeV);BR("+chi1+"#rightarrow#gamma)")
	Graph.GetHistogram().GetXaxis().SetRangeUser(100,800)
	Graph.GetHistogram().GetYaxis().SetRangeUser(200,1300)
	rt.gPad.SetRightMargin(0.15)
	Graph.Draw("hist colz")
	cd.SaveAs("output/M1_M2/BRg_physmass.pdf")
	
	cd2 = rt.TCanvas()
	diff_data = abs(M_char1)-abs(M_chi1)
	diff_wino = rt.TGraph()
	diff_bino = rt.TGraph()
	diff_same = rt.TGraph()
	j = 0
	k = 0
	l = 0
	for i in xrange(len(diff_data)):
		if(M1[i]>M2[i]):
			diff_wino.SetPoint(j,diff_data[i],BRg[i])
			j += 1
		elif(M1[i]<M2[i]):
			diff_bino.SetPoint(k,diff_data[i],BRg[i])
			k += 1
		elif(M1[i]==M2[i]):
			diff_same.SetPoint(l,diff_data[i],BRg[i])
			l += 1
	
	diff_bino.SetTitle(";m#kern[0.1]{_{"+char1+"}}-m#kern[0.1]{_{"+chi1+"}} (GeV);BR("+chi1+"#rightarrow#gamma)")
	diff_bino.SetMarkerColor(rt.kRed)
	diff_bino.GetYaxis().SetRangeUser(0,1.05)
	diff_bino.Draw("A*")
	diff_wino.Draw("* same")
	diff_same.SetMarkerColor(rt.kGreen)
	diff_same.Draw("* same")
	leg = rt.TLegend(0.65,0.15,0.85,0.4)
	leg.AddEntry(diff_bino,"Bino like (M1<M2)","P")
	leg.AddEntry(diff_wino,"Wino like (M1>M2)","P")
	leg.AddEntry(diff_same,"M1=M2","P")
	leg.SetBorderSize(0)
	leg.Draw("same")
	cd2.SaveAs("output/M1_M2/BRg_diff.pdf")
	cd2.SaveAs("output/M1_M2/BRg_diff.root")

def M1M3_BRphys():
	readout = []

	for f in os.listdir("../../susyhit/output/Bench_M1_M3/merged/"):
		if f[-5:] == ".slha":
			print f
			data = pyslha.readSLHAFile("../../susyhit/output/Bench_M1_M3/merged/"+f)
			massblock = data.blocks["MASS"]
			extparblock = data.blocks["EXTPAR"]
			decays_nlsp = data.decays[1000022]
			if len(decays_nlsp.decays)==1:
				readout.append([massblock[1000022],massblock[1000024],decays_nlsp.decays[0].br,0.,extparblock[1],extparblock[2]])
			else:
				if decays_nlsp.decays[0].ids[1]==22:
					brg = decays_nlsp.decays[0].br
					brz = decays_nlsp.decays[1].br
				else :
					brz = decays_nlsp.decays[0].br
					brg = decays_nlsp.decays[1].br
				readout.append([massblock[1000022],massblock[1000024],brg,brz,extparblock[1],extparblock[3]])
		
		#~ if (len(readout)>40): break

	def getKey(item):
		return (item[0],item[1])
		
	readout = sorted(readout,key=getKey)
	
	M_chi1 = np.array([item[0] for item in readout])
	M_char1 = np.array([item[1] for item in readout])
	BRg = np.array([item[2] for item in readout])
	BRz = np.array([item[3] for item in readout])
	M1 = np.array([item[4] for item in readout])
	M3 = np.array([item[5] for item in readout])
	
	Graph = rt.TGraph2D()
	for i in xrange(len(BRg)):
		Graph.SetPoint(i,M_chi1[i],M_char1[i],BRg[i])
	

	cd = rt.TCanvas()
	chi1 = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
	char1 = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
	Graph.GetHistogram().SetTitle(";m#kern[0.1]{_{"+chi1+"}} (GeV);m#kern[0.1]{_{"+char1+"}} (GeV);BR("+chi1+"#rightarrow#gamma)")
	#~ Graph.GetHistogram().GetXaxis().SetRangeUser(100,800)
	#~ Graph.GetHistogram().GetYaxis().SetRangeUser(200,1300)
	rt.gPad.SetRightMargin(0.15)
	Graph.Draw("hist colz")
	cd.SaveAs("output/M1_M3/BRg_physmass.pdf")	

def M1M2_BR_MCtest():
	data = np.loadtxt("input/BranchingRatios.txt",dtype={"names":("point","z","gamma","tot","mass"),"formats":("S50","i","i","i","f")},delimiter=",")
	pattern = "SUS-RunIISummer16FSPremix-M1(.*)_M2(.*).root"
	
	BR = []
	
	for i in xrange(len(data)):
		m = re.search(pattern,data[i][0])
		if (data[i][3] != 0):
			BR.append([int(m.groups()[0]),int(m.groups()[1]),data[i][1],data[i][2],data[i][3],data[i][4]])
		else:
			BR.append([int(m.groups()[0]),int(m.groups()[1]),0.8,0.2,1,data[i][4]])
	
	def getKey(item):
		return (item[0],item[1])
		
	BR = sorted(BR,key=getKey)
	
	BR_g = np.array([item[3] for item in BR])/(1.0*np.array([item[4] for item in BR]))
	BR_Z = np.array([item[2] for item in BR])/(1.0*np.array([item[4] for item in BR]))
	nNlsp_decay = np.array([item[4] for item in BR])
	mass_nlsp = np.array([item[5] for item in BR])
	
	for i in xrange(len(nNlsp_decay)):
		if (nNlsp_decay[i] == 1):
			nNlsp_decay[i] = 0
	
	
	BR_g = (BR_g.reshape((30,30))).T
	BR_Z = (BR_Z.reshape((30,30))).T
	nNlsp_decay = (nNlsp_decay.reshape((30,30))).T
	mass_nlsp = (mass_nlsp.reshape((30,30))).T
	
	print BR_g
	print BR_Z
	print nNlsp_decay
	print mass_nlsp
	
	plot_2d_M1M2(BR_g,"BR($\\tilde{\\chi}^0_1\\rightarrow\\gamma$)","Branching Ratio $\\tilde{\\chi}^0_1\\rightarrow\\gamma$","BR_gamma_MCtest")
	plot_2d_M1M2(BR_Z,"BR($\\tilde{\\chi}^0_1\\rightarrow Z$)","Branching Ratio $\\tilde{\\chi}^0_1\\rightarrow Z$","BR_z-MCtest")
	plot_2d_M1M2(nNlsp_decay/1000.,"N$_{NLSP-decays}$ per Event","N$_{NLSP-decays}$ per Event","nNLSP-MCtest")
	plot_2d_M1M2(mass_nlsp,"$M_{\\tilde{\\chi}^0_1}$ (GeV)","Mass $\\tilde{\\chi}^0_1$","Mass_neutralino_MCtest")
	
def M1M3_BR_MCtest():
	data = np.loadtxt("input/BranchingRatios_bench2.txt",dtype={"names":("point","z","gamma","tot","mass"),"formats":("S50","i","i","i","f")},delimiter=",")
	pattern = "SUS-RunIISummer16FSPremix-M1(.*)_M3(.*).root"
	
	BR = []
	
	for i in xrange(len(data)):
		m = re.search(pattern,data[i][0])
		if (data[i][3] != 0):
			BR.append([int(m.groups()[0]),int(m.groups()[1]),data[i][1],data[i][2],data[i][3],data[i][4]])
		else:
			BR.append([int(m.groups()[0]),int(m.groups()[1]),0.8,0.2,1,data[i][4]])
	
	def getKey(item):
		return (item[0],item[1])
		
	BR = sorted(BR,key=getKey)
	print BR
	
	BR_g = np.array([item[3] for item in BR])/(1.0*np.array([item[4] for item in BR]))
	BR_Z = np.array([item[2] for item in BR])/(1.0*np.array([item[4] for item in BR]))
	nNlsp_decay = np.array([item[4] for item in BR])
	mass_nlsp = np.array([item[5] for item in BR])
	
	for i in xrange(len(nNlsp_decay)):
		if (nNlsp_decay[i] == 1):
			nNlsp_decay[i] = 0
	
	
	BR_g = (BR_g.reshape((30,30))).T
	BR_Z = (BR_Z.reshape((30,30))).T
	nNlsp_decay = (nNlsp_decay.reshape((30,30))).T
	mass_nlsp = (mass_nlsp.reshape((30,30))).T
	
	print BR_g
	print BR_Z
	print nNlsp_decay
	print mass_nlsp
	
	plot_2d_M1M3(BR_g,"BR($\\tilde{\\chi}^0_1\\rightarrow\\gamma$)","Branching Ratio $\\tilde{\\chi}^0_1\\rightarrow\\gamma$","BR_gamma_MCtest",True)
	plot_2d_M1M3(BR_Z,"BR($\\tilde{\\chi}^0_1\\rightarrow Z$)","Branching Ratio $\\tilde{\\chi}^0_1\\rightarrow Z$","BR_z-MCtest",True)
	plot_2d_M1M3(nNlsp_decay/1000.,"N$_{NLSP-decays}$ per Event","N$_{NLSP-decays}$ per Event","nNLSP-MCtest",True)
	plot_2d_M1M3(mass_nlsp,"$M_{\\tilde{\\chi}^0_1}$ (GeV)","Mass $\\tilde{\\chi}^0_1$","Mass_neutralino_MCtest",True)

def M1M2_cross_prosp():
	data = np.loadtxt("input/Xsec_Scan1_M1_M2_merged.txt",dtype={"names":("M1","M2","LO","RelErr_LO","NLO","RelErr_NLO","LO_ms","Err_LO_ms","NLO_ms","Err_NLO_ms"),"formats":("i","i","f","f","f","f","f","S50","f","S50")})
	
	f = open("output/failed_M1_M2.txt","w")
	
	for temp in data:
		if temp[7] == "-nan":
			if temp[0] > 150 and temp[1] > 200:
				print temp[0],temp[1]
				f.write(str(temp[0])+" "+str(temp[1])+"\n")
			temp[7] = 0
		if temp[9] == "-nan": temp[9] = 0
		
	cross = []
	
	for i in xrange(len(data)):
		cross.append([data[i][0],data[i][1],np.maximum(1E-10,data[i][8]),data[i][9]])
	
	def getKey(item):
		return (item[0],item[1])
		
	cross = sorted(cross,key=getKey)
	
	xsec = open("output/M1_M2/xsec_GGM_M1_M2.txt","w")
	xsec.write("#M1 M2 xsec[pb] uncertainty[pb]\n")
	for i in xrange(len(cross)):
		if cross[i][0]>=200 and cross[i][1]>=200:
			if cross[i][2]!=1E-10:
				xsec.write(str(cross[i][0])+" "+str(cross[i][1])+" "+str(cross[i][2])+" "+str(cross[i][3])+"\n")
			else:
				xsec.write(str(cross[i][0])+" "+str(cross[i][1])+" "+str(0)+" "+str(cross[i][3])+"\n")
	
	xsec.close()
	
	M1 = np.array([item[0] for item in cross])
	M2 = np.array([item[1] for item in cross])
	NLO_ms = np.array([item[2] for item in cross])
	
	rfile = rt.TFile("output/M1_M2/hist2d.root","UPDATE")
	histo = rt.TH2F("NLO_ms_prospino",";M1 (GeV);M2 (GeV);xsec(pb)",30,25,1525,30,25,1525)
	for i in xrange(len(NLO_ms)):
		histo.Fill(M1[i],M2[i],NLO_ms[i])

	histo.Write("NLO_ms_prospino",rt.TObject.kOverwrite)	
	rfile.Close()
	
	NLO_ms = (NLO_ms.reshape((30,30))).T
	
	plot_2d_M1M2(NLO_ms,"NLO_ms (pb)","NLO_ms","NLO_ms_prospino")

def M1M3_cross_prosp():
	data = np.loadtxt("input/Xsec_Scan2_M1_M3_merged.txt",dtype={"names":("M1","M3","LO","RelErr_LO","NLO","RelErr_NLO","LO_ms","Err_LO_ms","NLO_ms","Err_NLO_ms"),"formats":("i","i","f","f","f","f","f","S50","f","S50")})
	
	f = open("output/failed_M1_M3.txt","w")
	
	for temp in data:
		if temp[7] == "-nan":
			print temp[0],temp[1]
			f.write(str(temp[0])+" "+str(temp[1])+"\n")
			temp[7] = 0
		if temp[9] == "-nan": temp[9] = 0
		
	f.close()
	cross = []
	
	for i in xrange(len(data)):
		cross.append([data[i][0],data[i][1],np.maximum(1E-10,data[i][8]),data[i][9]])
	
	def getKey(item):
		return (item[0],item[1])
		
	cross = sorted(cross,key=getKey)
	
	xsec = open("output/M1_M3/xsec_GGM_M1_M3.txt","w")
	xsec.write("#M1 M3 xsec[pb] uncertainty[pb]\n")
	for i in xrange(len(cross)):
		if cross[i][2]!=1E-10:
			xsec.write(str(cross[i][0])+" "+str(cross[i][1])+" "+str(cross[i][2])+" "+str(cross[i][3])+"\n")
		else:
			xsec.write(str(cross[i][0])+" "+str(cross[i][1])+" "+str(0)+" "+str(cross[i][3])+"\n")
	
	xsec.close()
	
	M1 = np.array([item[0] for item in cross])
	M3 = np.array([item[1] for item in cross])
	NLO_ms = np.array([item[2] for item in cross])
	
	rfile = rt.TFile("output/M1_M3/hist2d_scan2.root","UPDATE")
	histo = rt.TH2F("NLO_ms_prospino",";M1 (GeV);M3 (GeV);xsec(pb)",30,25,1525,31,975,2525)
	for i in xrange(len(NLO_ms)):
		histo.Fill(M1[i],M3[i],NLO_ms[i])

	histo.Write("NLO_ms_prospino",rt.TObject.kOverwrite)	
	rfile.Close()
	
	NLO_ms = (NLO_ms.reshape((30,31))).T
	
	plot_2d_M1M3(NLO_ms,"NLO_ms (pb)","NLO_ms","NLO_ms_prospino",False)

def M1M2_cross_physMass():
	readout = []

	for f in os.listdir("../../susyhit/output/Bench_M1_M2/merged/"):
		if f[-5:] == ".slha":
			print f
			data = pyslha.readSLHAFile("../../susyhit/output/Bench_M1_M2/merged/"+f)
			massblock = data.blocks["MASS"]
			extparblock = data.blocks["EXTPAR"]
			if extparblock[1]>=200 and extparblock[2]>=200:
				readout.append([massblock[1000022],massblock[1000024],extparblock[1],extparblock[2],0.])		
		#~ if (len(readout)>40): break

	def getKey(item):
		return (item[2],item[3])
		
	readout = sorted(readout,key=getKey)
	
	data = np.loadtxt("output/M1_M2/xsec_GGM_M1_M2.txt",dtype={"names":("M1","M2","Xsec"),"formats":("i","i","f")})
	
	i = 0
	for entry in data:
		readout[i][4] = entry[2]
		i += 1
	
	M_chi1 = np.array([item[0] for item in readout])
	M_char1 = np.array([item[1] for item in readout])
	M1 = np.array([item[2] for item in readout])
	M2 = np.array([item[3] for item in readout])
	xsec = np.array([item[4] for item in readout])
	
	Graph = rt.TGraph2D()
	for i in xrange(len(xsec)):
		Graph.SetPoint(i,M_chi1[i],M_char1[i],np.abs(xsec[i]))
	
	cd = rt.TCanvas()
	hist = Graph.GetHistogram()

	cd.SetLogz()
	chi1 = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
	char1 = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
	hist.SetTitle("Xsec GGM_M1_M2;m#kern[0.1]{_{"+chi1+"}} (GeV);m#kern[0.1]{_{"+char1+"}} (GeV);xsec[pb]")
	hist.GetZaxis().SetRangeUser(0.00001,10)
	rt.gPad.SetRightMargin(0.15)
	hist.Draw("colz")
	cd.SaveAs("output/M1_M2/Xsec_physmass.pdf")
	cd.SaveAs("output/M1_M2/Xsec_physmass.root")
	
def M1M3_cross_physMass():
	readout = []

	for f in os.listdir("../../susyhit/output/Bench_M1_M3/merged/"):
		if f[-5:] == ".slha":
			print f
			data = pyslha.readSLHAFile("../../susyhit/output/Bench_M1_M3/merged/"+f)
			massblock = data.blocks["MASS"]
			extparblock = data.blocks["EXTPAR"]
			readout.append([massblock[1000022],massblock[1000024],extparblock[1],extparblock[3],0.])		
		#~ if (len(readout)>40): break

	def getKey(item):
		return (item[2],item[3])
		
	readout = sorted(readout,key=getKey)
	
	data = np.loadtxt("output/M1_M3/xsec_GGM_M1_M3.txt",dtype={"names":("M1","M3","Xsec"),"formats":("i","i","f")})
	
	i = 0
	for entry in data:
		readout[i][4] = entry[2]
		i += 1
	
	M_chi1 = np.array([item[0] for item in readout])
	M_char1 = np.array([item[1] for item in readout])
	M1 = np.array([item[2] for item in readout])
	M2 = np.array([item[3] for item in readout])
	xsec = np.array([item[4] for item in readout])
	
	Graph = rt.TGraph2D()
	for i in xrange(len(xsec)):
		Graph.SetPoint(i,M_chi1[i],M_char1[i],np.abs(xsec[i]))
	
	hist = Graph.GetHistogram()
	cd = rt.TCanvas()
	
	cd.SetLogz()
	chi1 = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
	char1 = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
	hist.SetTitle("Xsec GGM_M1_M3;m#kern[0.1]{_{"+chi1+"}} (GeV);m#kern[0.1]{_{"+char1+"}} (GeV);xsec[pb]")
	hist.GetZaxis().SetRangeUser(0.00001,10)
	rt.gPad.SetRightMargin(0.15)
	Graph.Draw("hist colz")
	cd.SaveAs("output/M1_M3/Xsec_physmass.pdf")
	hist.SaveAs("output/M1_M3/Xsec_physmass.root")


#~ run_M1M3()
run_M1M2()
#~ M1M2_BRphys()
#~ M1M3_BRphys()
#~ M1M2_BR_MCtest()
#~ M1M3_BR_MCtest()
#~ M1M2_cross_prosp()
#~ M1M3_cross_prosp()
#~ M1M2_cross_physMass()
#~ M1M3_cross_physMass()






