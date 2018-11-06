import argparse
import os, re, subprocess as sp, math

def runCombine(datacard,selection,scan,single=False,sig=False):
	if sig: combineOutput=sp.check_output(["combine","-M","Significance","--uncapped","1","--rMin","-3","-n",(datacard.split("/")[2])[:-4],datacard])
	else: combineOutput=sp.check_output(["combine","-M","Asymptotic","--rMax","3","-n",(datacard.split("/")[2])[:-4],datacard])
	if single:
		print combineOutput
	if sig: directory="temp/Significance/"+selection+"/"+scan+"/"
	else: directory="temp/"+selection+"/"+scan+"/"
	if not os.path.exists(directory):
		os.makedirs(directory)
	f=open(directory+selection+"_"+"combineOut_"+datacard.split("/")[2],"w")
	f.write(combineOutput)
	f.close
	if sig: sp.check_output(["rm","higgsCombine"+(datacard.split("/")[2])[:-4]+".Significance.mH120.root"])
	else: sp.check_output(["rm","higgsCombine"+(datacard.split("/")[2])[:-4]+".Asymptotic.mH120.root"])
	return 0

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('files', default="", nargs="?")
	parser.add_argument('scan', nargs='?', help="choose a signal scan")
	parser.add_argument('selection', nargs='?', help="choose as selection like leptonVeto, htgVeto etc.")
	parser.add_argument('--significance', action="store_true", help="Calculated significance instead of limit")
	args = parser.parse_args()
	
	toDo=args.files.split(",")
	for x in toDo[:-1]:
		runCombine(x,args.selection,args.scan,sig=args.significance)
