#script to read in the prospino output and derive mean xsec and uncertainties
# Use CMSSW810 to be able to import matplotlib

import numpy as np
import ROOT as rt
import os,re
import math

def readOutput(M1,M3,rep,add=""):
	fileName="input/M1_M3/out_M1_"+str(M1)+"_M3_"+str(M3)+".slha/prospino_out_M1_"+str(M1)+"_M3_"+str(M3)+".slha_"+str(rep)+add+".dat"
	data=np.loadtxt(fileName,usecols=(1,2,5,pos,7,6,13),dtype={"names":("n1","n2","mu","xsec","chMass","nMass","k"),"formats":("i","i","d","f","d","d","S6")},comments="i1")
	return data

def xsecUnc_alphaS(M1,M3):
	xsecs=[]
	xsec_up=0
	xsec_down=0
	k = [[0. for i in range(101)] for j in range(13)]
	for rep in xrange(103):
		#~ if os.path.isfile("input/M1_M3/out_M1_"+str(M1)+"_M3_"+str(M3)+".slha/prospino_out_M1_"+str(M1)+"_M3_"+str(M3)+".slha_"+str(rep)+""+".dat")==False: continue
		data=readOutput(M1,M3,rep)
		if len(data)==0: continue
		tempxsec=0
		if rep<101:
			for process in data:
				if math.isnan(process[3]): continue
				if process[2]==1.0 : tempxsec+=process[3]
			xsecs.append(tempxsec)
			
			for n1 in [1,2,5]:
				for n2 in [5,7]:
					if n1==n2: continue
					for process in data:
						if process[0]==n1 and process[1]==n2:
							if process[2]==1.0 : 
								try:
									k[n1+n2][rep]=float(process[6])
								except ValueError:
									print "Not a float"
		
		elif rep==101:
			for process in data:
				if math.isnan(process[3]): continue
				if process[2]==1.0 : xsec_up+=process[3]
		elif rep==102:
			for process in data:
				if math.isnan(process[3]): continue
				if process[2]==1.0 : xsec_down+=process[3]	
		
		#~ print rep, tempxsec
	xsecs=np.array(xsecs)
	
	k_mean=np.zeros(13)
	i=0
	for process in k:
		process=np.array(process)
		k_mean[i]=process.mean()
		i+=1
	
	rootfile = rt.TFile("output/xsecs_distributions_M1_M3.root","UPDATE")
	histo = rt.TH1F(str(M1)+"_"+str(M3),";cross section (fb);nReplicas",10,xsecs.min(),xsecs.max())
	for value in xsecs:
		histo.Fill(value)
	histo.Write(str(M1)+"_"+str(M3),rt.TObject.kOverwrite)
	rootfile.Close()
	
	data=readOutput(M1,M3,0,"_scale")
	scale_up=0
	scale_down=0
	temp_down=0
	temp_up=0
	temp_mean=0
	n1_mass=0
	n2_mass=0
	#~ k=np.zeros(13)
	for n1 in [1,2,5]:
		for n2 in [5,7]:
			if n1==n2: continue
			for process in data:
				if process[0]==n1 and process[1]==n2:
					if process[2]==0.5 : temp_down=process[3]
					if process[2]==2.0 : temp_up=process[3]
					if process[2]==1.0 : 
						temp_mean=process[3]
						#~ try:
							#~ k[n1+n2]=float(process[6])
						#~ except ValueError:
							#~ print "Not a float"
					if n1==1 and n1_mass==0: n1_mass=process[5]
					elif n1==2 and n2_mass==0: n2_mass=process[5]
			scale_up=scale_up+(temp_up-temp_mean)**2
			scale_down=scale_down+(temp_down-temp_mean)**2
	
	return (xsecs.mean(), xsecs.std(ddof=1), np.abs(xsec_up-xsec_down)/2.0, np.median(xsecs), xsec_up, xsec_down, np.sqrt(scale_up), np.sqrt(scale_down),process[4],n1_mass,n2_mass,k_mean)

def getMasses(folder):
	pattern="out_M1_(.*)_M3_(.*).slha"
	m=re.search(pattern,folder)
	masses=[]
	if m and (len(m.groups())==2):
		for s in m.groups():
			masses.append(int(s))
	else:
		print "don't know what this is:",point
		exit(-1)
	return masses

if __name__ == "__main__":
	
	order="nlo"
	#~ order="lo"
	
	hist = rt.TH2F("xsecs",";M1 (GeV);M3 (GeV);cross section (pb)",30,25,1525,31,975,2525)
	hist_pdf = rt.TH2F("pdf_uncertainty",";M1 (GeV);M3 (GeV);pdf uncertainty (%)",30,25,1525,31,975,2525)
	hist_alpha = rt.TH2F("alphaS_uncertainty",";M1 (GeV);M3 (GeV);alpha_s uncertainty (%)",30,25,1525,31,975,2525)
	hist_diffMeanMed = rt.TH2F("difference_Mean_Median",";M1 (GeV);M3 (GeV);diff Mean/Median (%)",30,25,1525,31,975,2525)
	hist_alphaUp = rt.TH2F("alphaUp",";M1 (GeV);M3 (GeV);alphaUp (pb)",30,25,1525,31,975,2525)
	hist_alphaDown = rt.TH2F("alphaDown",";M1 (GeV);M3 (GeV);alphaDown (pb)",30,25,1525,31,975,2525)
	hist_alphaDiff = rt.TH2F("alphaDiff",";M1 (GeV);M3 (GeV);diff alphaUp alphaDown (pb)",30,25,1525,31,975,2525)
	hist_scaleUp = rt.TH2F("scaleUp_uncertainty",";M1 (GeV);M3 (GeV);scaleUp uncertainty (%)",30,25,1525,31,975,2525)
	hist_scaleDown = rt.TH2F("scaleDown_uncertainty",";M1 (GeV);M3 (GeV);scaleDown uncertainty (%)",30,25,1525,31,975,2525)
	hist_totalUnc = rt.TH2F("total_uncertainty",";M1 (GeV);M3 (GeV);total uncertainty (%)",30,25,1525,31,975,2525)
	hist_diffKN1 = rt.TH2F("diffKN1",";M1 (GeV);M3 (GeV);diff k (%)",30,25,1525,31,975,2525)
	hist_diffKN2 = rt.TH2F("diffKN2",";M1 (GeV);M3 (GeV);diff k (%)",30,25,1525,31,975,2525)
	hist_KN2C1 = rt.TH2F("KN2C1",";M1 (GeV);M2 (GeV);k",27,175,1525,31,975,2525)
	hist_KC1C1 = rt.TH2F("KC1C1",";M1 (GeV);M2 (GeV);k",27,175,1525,31,975,2525)
	gr_massPlane = rt.TGraph2D()
	gr_chMass = rt.TGraph()
	gr_massPlane2 = rt.TGraph2D()
	
	if order=="nlo": 
		out="output/xsec_GGM_M1_M3_PDF4LHC.txt"
		pos=15
	else:
		out="output/xsec_GGM_M1_M3_PDF4LHC_LO.txt"
		pos=14
	
	f=open(out,"w")
	f.write("#M1 M3 xsec[pb] uncertainty[pb]\n")
	
	for point in os.listdir("./input/M1_M3/"):
		if point=="missingPoints_M1M3.py": continue
		m=getMasses(point)
		#~ if m[0]!=200 or m[1]!=1500: continue
		result=xsecUnc_alphaS(m[0],m[1])
		hist.Fill(m[0],m[1],result[0])
		hist_pdf.Fill(m[0],m[1],result[1]/result[0]*100)
		hist_alpha.Fill(m[0],m[1],result[2]/result[0]*100)
		hist_diffMeanMed.Fill(m[0],m[1],abs(result[0]-result[3])/result[0]*100)
		hist_scaleUp.Fill(m[0],m[1],result[6]/result[0]*100)
		hist_scaleDown.Fill(m[0],m[1],result[7]/result[0]*100)
		hist_alphaUp.Fill(m[0],m[1],result[4]-result[0])
		hist_alphaDown.Fill(m[0],m[1],result[5]-result[0])
		hist_alphaDiff.Fill(m[0],m[1],(result[4]-result[5])/result[0]*100)
		hist_diffKN1.Fill(m[0],m[1],abs(result[11][12]-result[11][6])/((result[11][12]+result[11][6])/2.)*100)
		hist_diffKN2.Fill(m[0],m[1],abs(result[11][12]-result[11][7])/((result[11][12]+result[11][7])/2.)*100)
		hist_KN2C1.Fill(m[0],m[1],result[11][9])
		hist_KC1C1.Fill(m[0],m[1],result[11][12])
		
		scale_temp=max(result[6],result[7])
		tot_unc=np.sqrt(result[1]**2+result[2]**2+scale_temp**2)
		hist_totalUnc.Fill(m[0],m[1],tot_unc/result[0]*100)
		gr_chMass.SetPoint(gr_chMass.GetN(),result[8],result[0])
		gr_massPlane.SetPoint(gr_massPlane.GetN(),abs(result[9]),abs(result[10]),result[0])
		gr_massPlane2.SetPoint(gr_massPlane2.GetN(),abs(result[9]),abs(result[8]),result[0])
		
		f.write(str(m[0])+" "+str(m[1])+" "+str(result[0])+" "+str(tot_unc)+"\n")
		
		print m[0],m[1],result[11][12],result[11][6],result[11][7]
	
	if order=="nlo": outRoot="output/xsecs_M1_M3.root"
	else: outRoot="output/xsecs_M1_M3_LO.root"
	
	rfile = rt.TFile(outRoot,"UPDATE")	
	hist.Write("xsecs",rt.TObject.kOverwrite)
	hist_pdf.Write("pdf_uncertainty",rt.TObject.kOverwrite)
	hist_alpha.Write("alphaS_uncertainty",rt.TObject.kOverwrite)
	hist_diffMeanMed.Write("difference_Mean_Median",rt.TObject.kOverwrite)
	hist_alphaUp.Write("alphaUp",rt.TObject.kOverwrite)
	hist_alphaDown.Write("alphaDown",rt.TObject.kOverwrite)
	hist_alphaDiff.Write("alphaDiff",rt.TObject.kOverwrite)
	hist_scaleUp.Write("scaleUp_uncertainty",rt.TObject.kOverwrite)
	hist_scaleDown.Write("scaleDown_uncertainty",rt.TObject.kOverwrite)
	hist_totalUnc.Write("total_uncertainty",rt.TObject.kOverwrite)
	hist_diffKN1.Write("hist_diffKN1",rt.TObject.kOverwrite)
	hist_diffKN2.Write("hist_diffKN2",rt.TObject.kOverwrite)
	hist_KN2C1.Write("hist_KN2C1",rt.TObject.kOverwrite)
	hist_KC1C1.Write("hist_KC1C1",rt.TObject.kOverwrite)
	gr_chMass.Write("chMassVSxsec",rt.TObject.kOverwrite)
	gr_massPlane.Write("N1N2_MassPlane",rt.TObject.kOverwrite)
	gr_massPlane2.Write("N1C1_MassPlane",rt.TObject.kOverwrite)
	gr_massPlane2.GetHistogram().Write("N1C1_MassPlane_hist",rt.TObject.kOverwrite)

	f.close

	#~ rfile.Close()

	#~ print xsecUnc_alphaS(400,500)
