#Script to create significance histograms for the combination

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
		return float(line.split(":")[1])
	sig=0
	for line in output:
		if "Significance:" in line: sig = getLim(line)
	return sig

def readInOutput(datacard,scan,selection):
	my_file = "temp/Significance/"+selection+"/"+scan+"/"+selection+"_combineOut_"+datacard.split("/")[2]
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


	gr=rt.TGraph2D()
	if sScan=="T5Wg" or sScan=="T5Wg_thirds" or sScan=="T5gg": h_sig =rt.TH2F("","",18,0,2500,21,0,2150)
	elif sScan=="GGM_M1_M2": h_sig =rt.TH2F("","",27,175,1525,27,175,1525)
	#~ elif sScan=="GGM_M1_M2": h_sig =rt.TH2F("","",25,275,1525,25,275,1525)
	elif sScan=="GGM_M1_M3": h_sig =rt.TH2F("","",30,25,1525,31,975,2525)
	elif sScan=="TChiNg_BR": h_sig =rt.TH2F("","",51,-1,101,41,287.5,1312.5)
	elif sScan=="CharginoBR": h_sig =rt.TH2F("","",51,-1,101,41,287.5,1312.5)
	elif sScan=="CharginoBR_C1C1": h_sig =rt.TH2F("","",51,-1,101,40,287.5,1287.5)
	elif sScan=="CharginoBRstrong": h_sig =rt.TH2F("","",51,-1,101,21,0,2150)
	elif strongBR: h_sig =rt.TH2F("","",51,-1,101,21,0,2150)
	elif sScan=="T6Wg": h_sig =rt.TH2F("","",23,1000-25,2100+25,41,50,2100)
	elif sScan=="T6gg": h_sig =rt.TH2F("","",23,1000-25,2100+25,41,50,2100)
	elif sScan=="TChiWg": h_sig =rt.TH2F("","",41,300-12.5,1300+12.5,1,-1,1)
	elif sScan=="TChiNg": h_sig =rt.TH2F("","",41,300-12.5,1300+12.5,1,-1,1)
	elif sScan=="GGM1_NLSP500": h_sig =rt.TH2F("","",27,175,1525,1,-1,1)
	
	
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
		datacard="input/"+sScan+"_"+selection+"/"+point+".txt"
		sig=translateCombineOutput(readInOutput(datacard,sScan,selection))
		if sScan=="GGM_M1_M2" or sScan=="GGM_M1_M3" or sScan=="TChiNg_BR" or sScan=="CharginoBR" or sScan=="CharginoBR_C1C1":
			x,y=m1,m2
		else:
			x,y=m2,m1
		gr.SetPoint(i,x,y,sig)
		
		#Fix bug in scaled BR
		if sScan.find("BR")!=-1:
			if x==43: 
				x=42
		h_sig.SetBinContent(h_sig.FindBin(x,y),sig)
	
	f=rt.TFile(outdir+"significance_%s_"%sScan+selection+".root","update")
	gr.Write("gr_sig",rt.TObject.kOverwrite)

	h_sig.Write("h_sig",rt.TObject.kOverwrite)
	
	f.Close()
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('scan', nargs='?', help="choose a signal scan")
	parser.add_argument('selection', nargs='?', help="choose as selection like leptonVeto, htgVeto etc.")
	parser.add_argument('outdir', nargs='?', default="output/Significance/", help="output or test")
	args = parser.parse_args()
	
	missingCards=[]
	
	outdir=args.outdir
	selection=args.selection
	sScan=args.scan
	
	strongBR=False
	if sScan.find("CharginoBRstrong")!=-1: strongBR=True

	mergeLimits()
	
	for card in missingCards:
		sys.stdout.write(card+" ")
