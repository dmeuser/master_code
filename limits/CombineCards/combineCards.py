#Script to combine multiple datacards at once

import os
import sys
import multiprocessing
import tqdm
import argparse
from functools import partial
from MyDatacard import MyDatacard

def getMasses(path):
	masses=[]
	for f in os.listdir("./"+path):
		mass=[]
		splitted=f.split("_")
		for part in splitted:
			part=part.split(".")[0]
			if part.isdigit():
				mass.append(part)
				
		if len(mass)==2:
			masses.append(str(mass[0])+"_"+mass[1])
		else:
			print "Wrong syntax for masses"
	return masses

def combine(mass,path1,path2,output,ignore,corr):
	sScan=(path1.split("/")[1]).split("_")[0]
	if path1.find("GGM")!=-1:
		temp=(path1.split("/")[1]).split("_")
		sScan=temp[0]+"_"+temp[1]+"_"+temp[2]
		
	n1=path1+"/datacard_"+sScan+"_"+mass+".txt"
	
	if path2.find("knut")!=-1 and sScan=="T5Wg":
		n2=path2+"/Wg_"+mass+".txt"
	elif path2.find("lepton")!=-1:
		if sScan=="T5Wg":
			n2=path2+"/counting_t5Wg_"+mass+".txt"
		elif sScan=="CharginoBR":
			n2=path2+"/datacard_CharginoBR_"+mass+".txt"
	elif path2.find("htg_leptonVeto")!=-1 and sScan=="T5Wg":
		n2=path2+"/"+mass+".txt"
	elif path2.find("htgHigh")!=-1:
		n2=path2+"/"+mass+".txt"
	
	if len(corr)!=0:
		corrpairs=[]
		for i in range(len(corr)):
			if i%2!=0: continue
			corrpairs.append([corr[i],corr[i+1]])
		
		card1=MyDatacard(n1)
		card2=MyDatacard(n2)
		
		for uncNames in corrpairs:
			found1=card1.renameUncertainty(uncNames[0],uncNames[0]+"_CORR")
			found2=card2.renameUncertainty(uncNames[1],uncNames[0]+"_CORR")
			if found1==0 or found2==0: sys.exit(uncNames[0]+" and "+uncNames[1]+" not found")
		
		n1="temp/"+n1.split("/")[-1]
		n2="temp/"+n2.split("/")[-1]
		card1.write(n1)
		card2.write(n2)
		
	command="combineCards.py"
	for binName in ignore:
		command=command+" --xc="+binName
	
	if path2.find("knut")!=-1 or path2.find("htg")!=-1:
		command=command+" Photon_ST="+n1+" Photon_HTG="+n2+" >"+output+"/datacard_"+sScan+"_"+mass+".txt"
	elif path2.find("lepton")!=-1:
		command=command+" Photon_ST="+n1+" Photon_Lepton="+n2+" >"+output+"/datacard_"+sScan+"_"+mass+".txt"
	
	if os.system(command)!=0:
		sys.exit(sScan+"_"+mass+" not succesfull")
	
	#~ if len(corr)!=0:
		#~ os.system("rm temp/*")
	return 0
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('path1', nargs='?', help="choose folder for first datacard set (has to start with scan name)")
	parser.add_argument('path2', nargs='?', help="choose folder for second datacard set")
	parser.add_argument('outdir', nargs='?', help="choose folder for output")
	parser.add_argument('--ignoreBins', type=str, default=[""], nargs="+", help="Name of bins which are ignored")
	parser.add_argument('--correlate', type=str, default=[], nargs="+", help="Uncertainties which are assumend to be correlated")
	parser.add_argument('--single', action="store_true", help="Only a Cardset is combined")
	args = parser.parse_args()
	
	if not os.path.exists(args.outdir):
		os.makedirs(args.outdir)
	
	if args.path1.endswith("/"):
		args.path1=args.path1[:-1]
	if args.path2.endswith("/"):
		args.path2=args.path2[:-1]
	if args.outdir.endswith("/"):
		args.outdir=args.outdir[:-1]
	
	masses=[]
	masses1=getMasses(args.path1)
	masses2=getMasses(args.path2)
	
	for mass in masses1:
			if mass in masses2:
				masses.append(mass)
		
	#Clean sys.argv to avoid error due to MyDatcard		
	path1=args.path1
	path2=args.path2
	outdir=args.outdir
	ignoreBins=args.ignoreBins
	correlate=args.correlate
	single=args.single
	sys.argv=[]
	
	if single:
		for mass in masses:
			print mass
			combine(mass,path1,path2,outdir,ignoreBins,correlate)
			break
	else:
		p = multiprocessing.Pool()
		for _ in tqdm.tqdm(p.imap_unordered(partial(combine,path1=path1,path2=path2,output=outdir,ignore=ignoreBins,corr=correlate),masses), total=len(masses)):
			pass
