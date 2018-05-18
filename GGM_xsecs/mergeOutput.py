#script to read in the prospino output and derive mean xsec and uncertainties
# Use CMSSW810 to be able to import matplotlib

import numpy as np
import matplotlib.pyplot as plt
import ROOT as rt
import os,re

def readOutput(M1,M2,rep):
	fileName="input/out_M1_"+str(M1)+"M2_"+str(M2)+".slha/prospino_out_M1_"+str(M1)+"M2_"+str(M2)+".slha_"+str(rep)+".dat"
	data=np.loadtxt(fileName,usecols=(1,2,5,15),dtype={"names":("n1","n2","mu","xsec"),"formats":("i","i","d","f")},comments="i1")
	return data
	
def xsecUnc_prospinoScale(M1,M2):
	xsecs=[]
	xsecs_05=[]
	xsecs_20=[]
	for rep in xrange(101):
		data=readOutput(M1,M2,rep)
		if len(data)==0: continue
		tempxsec=0
		tempxsec_05=0
		tempxsec_20=0
		for process in data:
			if process[2]==1.0 : tempxsec+=process[3]
			if process[2]==0.5 : tempxsec_05+=process[3]
			if process[2]==2.0 : tempxsec_20+=process[3]
		xsecs.append(tempxsec)
		xsecs_05.append(tempxsec_05)
		xsecs_20.append(tempxsec_20)
	xsecs=np.array(xsecs)
	xsecs_05=np.array(xsecs_05)
	xsecs_20=np.array(xsecs_20)
	plt.figure()
	plt.hist(xsecs)
	plt.savefig("test.pdf")
	return (xsecs.mean(), xsecs.std(ddof=1)), xsecs_20.mean(), xsecs_05.mean(), (xsecs_20.mean()-xsecs_05.mean())/2.0

def xsecUnc_alphaS(M1,M2):
	xsecs=[]
	xsec_up=0
	xsec_down=0
	for rep in xrange(103):
		data=readOutput(M1,M2,rep)
		if len(data)==0: continue
		tempxsec=0
		if rep<101:
			for process in data:
				if process[2]==1.0 : tempxsec+=process[3]
			xsecs.append(tempxsec)
		elif rep==101:
			for process in data:
				if process[2]==1.0 : xsec_up+=process[3]
		elif rep==102:
			for process in data:
				if process[2]==1.0 : xsec_down+=process[3]
			
	xsecs=np.array(xsecs)
	
	rootfile = rt.TFile("output/xsecs_distributions_M1_M2.root","UPDATE")
	histo = rt.TH1F(str(M1)+"_"+str(M2),";cross section (fb);nReplicas",10,xsecs.min(),xsecs.max())
	for value in xsecs:
		histo.Fill(value)
	histo.Write(str(M1)+"_"+str(M2),rt.TObject.kOverwrite)
	rootfile.Close()

	return (xsecs.mean(), xsecs.std(ddof=1), np.abs(xsec_up-xsec_down)/2.0)

def getMasses(folder):
	pattern="out_M1_(.*)M2_(.*).slha"
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
	
	hist = rt.TH2F("xsecs",";M1 (GeV);M2 (GeV);cross section (pb)",30,25,1525,30,25,1525)
	hist_pdf = rt.TH2F("pdf_uncertainty",";M1 (GeV);M2 (GeV);pdf uncertainty (%)",30,25,1525,30,25,1525)
	hist_scale = rt.TH2F("scale_uncertainty",";M1 (GeV);M2 (GeV);scale uncertainty (%)",30,25,1525,30,25,1525)
	
	for point in os.listdir("./input/"):
		m=getMasses(point)
		if m[1]==150 or m[1]==200 or m[0]<200: continue
		result=xsecUnc_alphaS(m[0],m[1])
		hist.Fill(m[0],m[1],result[0])
		hist_pdf.Fill(m[0],m[1],result[1]/result[0]*100)
		hist_scale.Fill(m[0],m[1],result[2]/result[0]*100)
		print m[0],m[1],result
	
	rfile = rt.TFile("output/xsecs_M1_M2.root","UPDATE")	
	hist.Write("xsecs",rt.TObject.kOverwrite)
	hist_pdf.Write("pdf_uncertainty",rt.TObject.kOverwrite)
	hist_scale.Write("scale_uncertainty",rt.TObject.kOverwrite)
	#~ rfile.Close()



#~ print xsecUnc_alphaS(1000,1500)
#~ print xsecUnc(1200,700)
