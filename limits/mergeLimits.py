#Script to create limit histograms for the combination

import os, re, subprocess as sp, math
from distutils import spawn
import ROOT as rt
import argparse
import sys


def getMasses(point):
	if sScan=="T5Wg" : pattern="datacard_T5Wg_(.*)_(.*)"
	elif sScan=="T5Wg_thirds" : pattern="datacard_T5Wg_(.*)_(.*)"
	elif sScan=="GGM_M1_M2" : pattern="datacard_GGM_M1_M2_(.*)_(.*)"
	elif sScan=="GGM1_NLSP500" : pattern="datacard_GGM_M1_M2_(.*)_(.*)"
	elif sScan=="GGM_M1_M3" : pattern="datacard_GGM_M1_M3_(.*)_(.*)"
	elif sScan=="TChiNg_BR" : pattern="datacard_TChiNg_BR_(.*)_(.*)"
	elif sScan=="CharginoBR" : pattern="datacard_CharginoBR_(.*)_(.*)"
	elif sScan=="CharginoBR_C1C1" : pattern="datacard_CharginoBR_(.*)_(.*)"
	elif sScan=="CharginoBRstrong" : pattern="datacard_CharginoBRstrong1950_(.*)_(.*)"
	elif sScan=="T6Wg" : pattern="datacard_T6Wg_(.*)_(.*)"
	elif sScan=="T6gg" : pattern="datacard_T6gg_(.*)_(.*)"
	elif sScan=="T5gg" : pattern="datacard_T5gg_(.*)_(.*)"
	elif sScan=="TChiWg": pattern="datacard_TChiWG_(.*)"
	#~ elif sScan=="TChiNg": pattern="datacard_TChiNG_(.*)"
	elif sScan=="TChiNg": pattern="datacard_TChiNg_(.*)_(.*)"
	elif sScan.find("CharginoBRstrongN")!=-1 : pattern="datacard_CharginoBRstrongN(.*)_(.*)_(.*)"
	elif sScan.find("CharginoBRstrongG")!=-1 : pattern="datacard_CharginoBRstrongG(.*)_(.*)_(.*)"
	if sScan=="T5Wg" and selection=="lepton" and selection!="lepton_final": pattern="counting_t5Wg_(.*)_(.*)"
	if sScan=="T5Wg" and (selection=="htg" or selection=="htgNN" or selection=="htgFinal") : pattern="(.*)_(.*)"
	if sScan=="T5Wg_thirds" and (selection=="htg" or selection=="htgNN" or selection=="htgFinal") : pattern="(.*)_(.*)"
	if sScan=="T6gg" and (selection=="htg" or selection=="htgNN" or selection=="htgFinal") : pattern="(.*)_(.*)"
	if sScan=="T5gg" and (selection=="htg" or selection=="htgNN" or selection=="htgFinal") : pattern="(.*)_(.*)"
	if sScan=="T6Wg" and (selection=="htg" or selection=="htgNN" or selection=="htgFinal") : pattern="(.*)_(.*)"
	if sScan=="GGM_M1_M2" and (selection=="htg" or selection=="htgNN" or selection=="htgFinal") : pattern="(.*)_(.*)"
	if sScan=="GGM_M1_M3" and (selection=="htg" or selection=="htgNN" or selection=="htgFinal") : pattern="(.*)_(.*)"
	if sScan=="GGM_M1_M3" and selection=="lepton" : pattern="counting_GMSB_(.*)_(.*)"
	#~ if sScan=="GGM_M1_M3" and selection=="diphoton_final" : pattern="counting_GGM_(.*)_(.*)"
	m=re.search(pattern,point)
	masses=[]
	print point
	if m and (len(m.groups())==2) and strongBR==False:
		for s in m.groups():
			masses.append(int(s))
	elif m and (len(m.groups())==3) and strongBR:
		masses.append(int(m.groups()[1]))
		masses.append(int(m.groups()[2]))
		masses.append(int(m.groups()[0]))
	elif m and (sScan=="TChiNg" or sScan=="TChiWg") and (len(m.groups())==1):
		masses=[int(m.groups()[0]),0]
	else:
		print "don't know what this is:",point
		exit(-1)
	return masses

	
def translateCombineOutput(output):
	def getLim(line):
		return float(line.split("<")[1])
	limit={}
	limit["obs"] = 0
	limit["exp"] = 0
	limit["exp+1"] = 0
	limit["exp-1"] = 0
	limit["exp+2"] = 0
	limit["exp-2"] = 0
	for line in output:
		if "Observed Li" in line: limit["obs"] = getLim(line)
		if "Expected 50" in line: limit["exp"] = getLim(line)
		if "Expected 84" in line: limit["exp+1"] = getLim(line)
		if "Expected 16" in line: limit["exp-1"] = getLim(line)
		if "Expected 97" in line: limit["exp+2"] = getLim(line)
		if "Expected  2" in line: limit["exp-2"] = getLim(line)
	return limit

def readInOutput(datacard,scan,selection):
	my_file = "temp/"+selection+"/"+scan+"/"+selection+"_combineOut_"+datacard.split("/")[2]
	output=""
	if os.path.isfile(my_file):
		f=open(my_file,"r")
		output=f.readlines()
	else:
		missingCards.append(datacard.split("/")[-1])
	return output


def mergeLimits():
	print "using",os.environ["CMSSW_VERSION"]
	print "combine tool",spawn.find_executable("combine")
	print

	if sScan=="T5Wg": xsecFile="input/xsec_gluglu.csv"
	elif sScan=="GGM_M1_M2": xsecFile="input/xsec_GGM_M1_M2.txt"
	elif sScan=="GGM_M1_M3": xsecFile="input/xsec_GGM_M1_M3.txt"
	elif sScan=="TChiNg_BR": xsecFile="input/xsec_comb_wino.csv"
	elif sScan=="TChiNg": xsecFile="input/xsec_comb_wino.csv"
	elif sScan=="CharginoBR": xsecFile="input/xsec_N2C1_wino.csv"
	elif sScan=="CharginoBR_C1C1": xsecFile="input/xsec_C1C1_wino.csv"
	elif sScan=="TChiWg": xsecFile="input/xsec_N2C1_wino.csv"
	elif strongBR: xsecFile="input/xsec_gluglu.csv"
	elif sScan=="T6Wg": xsecFile="input/xsec_sqsq.csv"
	elif sScan=="T6gg": xsecFile="input/xsec_sqsq.csv"
	elif sScan=="T5gg": xsecFile="input/xsec_gluglu.csv"
	elif sScan=="T5Wg_thirds": xsecFile="input/xsec_gluglu.csv"
	elif sScan=="GGM1_NLSP500": xsecFile="input/xsec_GGM_M1_M2.txt"

	points=[]
	i=0
	m1_lim=0
	m2_lim=0
	if (sScan=="T5Wg" or sScan=="T5gg" or sScan=="T5Wg_thirds") and selection!="leptonVeto_test":
		m2_lim=1100
	elif sScan=="GGM_M1_M2":
		m1_lim=200
		m2_lim=200
	#~ elif sScan=="GGM_M1_M3":
		#~ m1_lim=150
	for p in os.listdir("./input/"+sScan+"_"+selection+""):
		p=p.split(".")[0]
		if sScan=="T5Wg" or sScan=="T5Wg_thirds" or sScan=="TChiNg_BR" or sScan=="CharginoBR" or sScan=="CharginoBR_C1C1" or sScan=="T5gg" or sScan=="CharginoBRstrong":
			m2,m1=getMasses(p)
		elif strongBR:
			m=getMasses(p)
			m2=m[0]
			m1=m[1]
			m_set=m[2]
		else:
			m1,m2=getMasses(p)
		if m2<m2_lim or m1<m1_lim: continue
		if sScan=="CharginoBR" and selection=="lepton" and m1==0: continue
		points.append(p)
		i+=1

	xsec={}
	xsec_err={}
	if sScan=="GGM_M1_M2" or sScan=="GGM_M1_M3" or sScan=="GGM1_NLSP500":
		with open(xsecFile) as f:
			for line in f:
				if line.startswith("#"): continue
				xsec[int(line.split()[0])*100000+int(line.split()[1])]=float(line.split()[2])
				if float(line.split()[3])!=0:
					xsec_err[int(line.split()[0])*100000+int(line.split()[1])]=float(line.split()[3])/float(line.split()[2])
				else:
					xsec_err[int(line.split()[0])*100000+int(line.split()[1])]=0
					
	elif sScan=="TChiNg_BR" or sScan=="CharginoBR" or sScan=="CharginoBR_C1C1":
		with open(xsecFile) as f:
			for line in f:
				if line.startswith("#"): continue
				xsec[int(line.split()[0])]=float(line.split()[1])
				xsec_err[int(line.split()[0])]=float(line.split()[2])/float(line.split()[1])
	else:
		with open(xsecFile) as f:
			for line in f:
				if line.startswith("#"): continue
				xsec[int(line.split()[0])]=float(line.split()[1])
				xsec_err[int(line.split()[0])]=float(line.split()[2])/100.
	
	#~ print xsec[0]
	#~ print xsec_err

	gr={}
	for lvl in ["obs","obs+1","obs-1","exp","exp+1","exp-1","exp+2","exp-2","obs_xs"]:
		gr[lvl]=rt.TGraph2D()
	if sScan=="T5Wg" or sScan=="T5Wg_thirds" or sScan=="T5gg": h_exp =rt.TH2F("","",18,0,2500,21,0,2150)
	elif sScan=="GGM_M1_M2": h_exp =rt.TH2F("","",27,175,1525,27,175,1525)
	#~ elif sScan=="GGM_M1_M2": h_exp =rt.TH2F("","",25,275,1525,25,275,1525)
	elif sScan=="GGM_M1_M3": h_exp =rt.TH2F("","",30,25,1525,31,975,2525)
	elif sScan=="TChiNg_BR": h_exp =rt.TH2F("","",51,-1,101,41,287.5,1312.5)
	elif sScan=="CharginoBR": h_exp =rt.TH2F("","",51,-1,101,41,287.5,1312.5)
	elif sScan=="CharginoBR_C1C1": h_exp =rt.TH2F("","",51,-1,101,40,287.5,1287.5)
	elif sScan=="CharginoBRstrong": h_exp =rt.TH2F("","",51,-1,101,21,0,2150)
	elif strongBR: h_exp =rt.TH2F("","",51,-1,101,21,0,2150)
	elif sScan=="T6Wg": h_exp =rt.TH2F("","",23,1000-25,2100+25,41,50,2100)
	elif sScan=="T6gg": h_exp =rt.TH2F("","",23,1000-25,2100+25,41,50,2100)
	elif sScan=="TChiWg": h_exp =rt.TH2F("","",41,300-12.5,1300+12.5,1,-1,1)
	elif sScan=="TChiNg": h_exp =rt.TH2F("","",41,300-12.5,1300+12.5,1,-1,1)
	elif sScan=="GGM1_NLSP500": h_exp =rt.TH2F("","",27,175,1525,1,-1,1)
	
	h_obs   =rt.TH2F(h_exp)
	h_obs_m1   =rt.TH2F(h_exp)
	h_obs_p1   =rt.TH2F(h_exp)
	h_exp_m1   =rt.TH2F(h_exp)
	h_exp_p1   =rt.TH2F(h_exp)
	h_exp_m2   =rt.TH2F(h_exp)
	h_exp_p2   =rt.TH2F(h_exp)
	h_exp_xs=rt.TH2F(h_exp)
	h_obs_xs=rt.TH2F(h_exp)
	#print points
	
	for i,point in enumerate(points):
		print point
		if sScan=="GGM_M1_M2" or sScan=="GGM_M1_M3" or sScan=="GGM1_NLSP500":
			m1,m2=getMasses(point)
			key=m1*100000+m2
		elif sScan=="CharginoBRstrong":
			m1,m2=getMasses(point)
			key=1950
		elif sScan.find("CharginoBRstrongN")!=-1:
			m=getMasses(point)
			m1=m[0]
			m2=m[1]
			key=m[0]
		elif sScan.find("CharginoBRstrongG")!=-1:
			m=getMasses(point)
			m1=m[0]
			m2=m[1]
			key=m[2]
		else :
			m2,m1=getMasses(point)
			key=m2
		#~ if sScan=="CharginoBR" or sScan=="CharginoBR_C1C1":
			#~ m1=100-m1
		#~ if strongBR:
			#~ m2=100-m2
		xs=xsec[key]
		datacard="input/"+sScan+"_"+selection+"/"+point+".txt"
		rLimits=translateCombineOutput(readInOutput(datacard,sScan,selection))
		e_xs_rel=xsec_err[key]
		rLimits["obs+1"]=rLimits["obs"]*(1+e_xs_rel)
		rLimits["obs-1"]=rLimits["obs"]*(1-e_xs_rel)
		if sScan=="GGM_M1_M2" or sScan=="GGM_M1_M3" or sScan=="TChiNg_BR" or sScan=="CharginoBR" or sScan=="CharginoBR_C1C1":
			x,y=m1,m2
		else:
			x,y=m2,m1
		for lvl in ["obs","obs+1","obs-1","exp","exp+1","exp-1","exp+2","exp-2"]:
			rLim=rLimits[lvl]
			gr[lvl].SetPoint(i,x,y,rLim)
		
		#Fix bug in scaled BR
		if sScan.find("BR")!=-1:
			if x==43: 
				x=42
		rLim=rLimits['exp']
		h_exp.SetBinContent(h_exp.FindBin(x,y),rLim)
		h_exp_xs.SetBinContent(h_exp.FindBin(x,y),rLim*xs)
		#~ print "exp. signal strength limit",rLim
		rLim=rLimits['obs']
		h_obs.SetBinContent(h_obs.FindBin(x,y),rLim)
		h_obs_xs.SetBinContent(h_obs.FindBin(x,y),rLim*xs)
		gr["obs_xs"].SetPoint(i,x,y,rLim*xs)
		#~ print "obs. signal strength limit",rLim
		
		rLim=rLimits['obs+1']
		h_obs_p1.SetBinContent(h_obs.FindBin(x,y),rLim)
		rLim=rLimits['obs-1']
		h_obs_m1.SetBinContent(h_obs.FindBin(x,y),rLim)
		rLim=rLimits['exp+1']
		h_exp_p1.SetBinContent(h_obs.FindBin(x,y),rLim)
		rLim=rLimits['exp-1']
		h_exp_m1.SetBinContent(h_obs.FindBin(x,y),rLim)
		rLim=rLimits['exp+2']
		h_exp_p2.SetBinContent(h_obs.FindBin(x,y),rLim)
		rLim=rLimits['exp-2']
		h_exp_m2.SetBinContent(h_obs.FindBin(x,y),rLim)
	
	f=rt.TFile(outdir+"limits_%s_"%sScan+selection+".root","update")
	for lvl in ["obs","obs+1","obs-1","exp","exp+1","exp-1","exp+2","exp-2","obs_xs"]:
		gr[lvl].Write("gr_"+lvl,rt.TObject.kOverwrite)

	h_exp.Write("h_exp",rt.TObject.kOverwrite)
	h_exp_xs.Write("h_exp_xs",rt.TObject.kOverwrite)
	h_obs.Write("h_obs",rt.TObject.kOverwrite)
	h_obs_xs.Write("h_obs_xs",rt.TObject.kOverwrite)
	
	h_obs_p1.Write("h_obs+1",rt.TObject.kOverwrite)
	h_obs_m1.Write("h_obs-1",rt.TObject.kOverwrite)
	h_exp_p1.Write("h_exp+1",rt.TObject.kOverwrite)
	h_exp_m1.Write("h_exp-1",rt.TObject.kOverwrite)
	h_exp_p2.Write("h_exp+2",rt.TObject.kOverwrite)
	h_exp_m2.Write("h_exp-2",rt.TObject.kOverwrite)
	
	f.Close()

def getContour(gr2,lvl=""):
	""" inspired by https://hypernews.cern.ch/HyperNews/CMS/get/susy/2122/2.html
	"""
	c=rt.TCanvas()
	gr2.Draw("tri1")
	c.Update()
	cont=gr2.GetContourList(1.)
	print cont
	if sScan=="GGM_M1_M2" and lvl=="exp" or (sScan=="GGM_M1_M2" and selection=="allCombined_final") or (sScan=="CharginoBR_C1C1" and (selection=="htgNN" or selection=="htgFinal") and lvl=="exp" and cont):
		i=1
		for ci in cont:
			if selection=="allCombined_final" and lvl.find("2")==-1:
					if lvl.find("exp")!=-1 and i==2:
						ci.Write("xtra_multicont_gr_"+lvl,rt.TObject.kOverwrite)
					if lvl.find("obs")!=-1 and i==3:
						ci.Write("xtra_multicont_gr_"+lvl,rt.TObject.kOverwrite)
			
			if lvl=="exp":
				ci.Write("multicont_gr_"+str(i),rt.TObject.kOverwrite)
			i+=1
		if lvl=="exp":	
			global numMulti
			numMulti=i
	
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
		
def getContourHS(hist,lvl=""):
	
	print hist
	for i in xrange(hist.GetNbinsX()+1):
		for j in xrange(hist.GetNbinsY()+1):
			if hist.GetBinContent(i,j)<1 :
				hist.SetBinContent(i,j,0)
			else :
				hist.SetBinContent(i,j,2)
	
	c=rt.TCanvas()
	c.cd()
	hist.SetContour(1)
	hist.SetContourLevel(0,1.)
	hist.Draw("CONT Z LIST")
	c.Update()
	conts=rt.gROOT.GetListOfSpecials().FindObject("contours")
	cont=conts.At(0)
	
	if sScan=="GGM_M1_M2" and lvl=="exp":
		i=1
		for ci in cont:
			ci.Write("multicont_gr_"+str(i),rt.TObject.kOverwrite)
			i+=1
		global numMulti
		numMulti=i
	
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

def getContourInter(gr):
	hsInter=rt.TH2F("","",150,gr.GetXaxis().GetXmin(),gr.GetXaxis().GetXmax(),150,gr.GetYaxis().GetXmin(),gr.GetYaxis().GetXmax())
	for i in range(1,hsInter.GetXaxis().GetNbins()+1):
		for j in range(1,hsInter.GetYaxis().GetNbins()+1):
			hsInter.SetBinContent(i,j,gr.Interpolate(hsInter.GetXaxis().GetBinCenter(i),hsInter.GetYaxis().GetBinCenter(j)))
	return getContourHS(hsInter),hsInter

def getContours(flag=""):
	f=rt.TFile(outdir+"limits_%s_"%sScan+selection+".root","update")
	for lvl in ["obs","obs+1","obs-1","exp","exp+1","exp-1","exp+2","exp-2"]:
		
		if flag=="hist":
			print "hist contouring"
			hs=f.Get("h_"+lvl)
			grC=getContourHS(hs,lvl)
		elif flag=="inter" :
			print "inter contouring"
			gr=f.Get("gr_"+lvl)
			grC,hsInter=getContourInter(gr)
			hsInter.Write("hs_inter"+lvl,rt.TObject.kOverwrite)
		else :
			gr=f.Get("gr_"+lvl)
			grC=getContour(gr,lvl)
		if grC:
			grC.Write("gr_"+lvl+"C",rt.TObject.kOverwrite)
		else:
			print "could not get contour for",lvl
	f.Close()
	
def redoHistogram():
	f=rt.TFile(outdir+"limits_%s_"%sScan+selection+".root","update")
	gr=f.Get("gr_obs_xs")
	h=gr.GetHistogram()
	h.Write("h_obs_xs_redone",rt.TObject.kOverwrite)
	f.Close()

def smoothContour_knut(gr, neighbors=5, sigma=0.5):
#~ def smoothContour_knut(gr, neighbors=5, sigma=0.3):
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
    if sScan=="T5Wg":
		if selection=="ST_HTG_lepton_final":
			n=n-5
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

def smoothContours(flag):
	f=rt.TFile(outdir+"limits_%s_"%sScan+selection+".root","update")
	gs=rt.TGraphSmooth()
	for lvl in ["obs","obs+1","obs-1","exp","exp+1","exp-1","exp+2","exp-2"]:
		gr=f.Get("gr_"+lvl+"C")
		#~ gr_sm=gs.SmoothLowess(gr)
		#~ gr_sm=gs.Approx(gr)
		#~ gr_sm=gs.SmoothSuper(gr)
		gr_sm=smoothContour_knut(gr)
		gr_sm.Write("gr_"+lvl+"C_sm",rt.TObject.kOverwrite)
	if flag=="multi":
		for num in range(1,numMulti):
			gr=f.Get("multicont_gr_"+str(num))
			print "multicont_gr_"+str(num)
			gr_sm=smoothContour_knut(gr)
			gr_sm.Write("multicont_gr_"+str(num)+"_sm",rt.TObject.kOverwrite)
		if sScan=="GGM_M1_M2" and selection=="allCombined_final":
			for lvl in ["obs","obs+1","obs-1","exp","exp+1","exp-1"]:
				gr=f.Get("xtra_multicont_gr_"+lvl)
				print "xtra_multicont_gr_"+lvl
				gr_sm=smoothContour_knut(gr)
				gr_sm.Write("xtra_multicont_gr_"+lvl+"_sm",rt.TObject.kOverwrite)
	f.Close()


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('scan', nargs='?', help="choose a signal scan")
	parser.add_argument('selection', nargs='?', help="choose as selection like leptonVeto, htgVeto etc.")
	parser.add_argument('outdir', nargs='?', default="output/", help="output or test")
	parser.add_argument('--contourtype', type=str, default="", help="hist or inter (inter uses finer binning than default) or mulit for multiple graphs")
	args = parser.parse_args()
	
	missingCards=[]
	
	outdir=args.outdir
	selection=args.selection
	sScan=args.scan
	contourtype=args.contourtype
	
	strongBR=False
	if sScan.find("CharginoBRstrong")!=-1: strongBR=True
	
	print contourtype

	mergeLimits()
	if sScan!="TChiNg" and sScan!="GGM1_NLSP500":
		getContours(contourtype)
		redoHistogram()
		smoothContours(contourtype)
	
	for card in missingCards:
		sys.stdout.write(card+" ")


