#Script to create datacards for different BRs

from MyDatacard import MyDatacard
import numpy as np
import argparse
import os
import sys
import pickle

def getMasses(path):
	mass=0
	splitted=path.split("_")
	for part in splitted:
		part=part.split(".")[0]
		if part.isdigit():
			mass=part
			break
	return int(mass)
	
def getXsec( mother_mass, pklfilename):
    info = (0,0)
    with open( pklfilename, 'rb') as f:
        data = pickle.load( f )
        if mother_mass in data:
            info = data[mother_mass]
        else:
            print "could not find mass {} in file {}".format(mother_mass, pklfilename)
    return info

def Scan_NeutralinoBr(selection):
		
		for f in os.listdir("input/ST/"+selection):
			path="input/ST/"+selection+"/"+f
			outputDC = MyDatacard(path)
			print f
			mass=getMasses(f)
			oldXsec=getXsec(mass,"xSec_SMS_N2C1_13TeV.pkl")[0]
			newXsec=getXsec(mass,"xSec_SMS_C1C1_13TeV.pkl")[0]
			ratio=newXsec/oldXsec
			
			
			for nBin in ["bin1","bin2","bin3","bin4"]:
				newYield = outputDC.exp[nBin]["sig"]*ratio
				
				outputDC.newSignalYield({nBin: newYield},"sig")
			
			directory="output/ST/"+selection
			if not os.path.exists(directory):
				os.makedirs(directory)
			
			outputDC.write(filename=directory+"/"+f)
			
def Scan_NeutralinoBrHTG(selection):
		
		for f in os.listdir("input/HTG/"+selection):
			path="input/HTG/"+selection+"/"+f
			outputDC = MyDatacard(path)
			print f
			mass=getMasses(f)
			oldXsec=getXsec(mass,"xSec_SMS_N2C1_13TeV.pkl")[0]
			newXsec=getXsec(mass,"xSec_SMS_C1C1_13TeV.pkl")[0]
			ratio=newXsec/oldXsec
			
			bins = ["binlowEMHT_24","binlowEMHT_25","binlowEMHT_26","binhighEMHT_24","binhighEMHT_25","binhighEMHT_26"]
			if selection=="highHTG" or selection=="DILEPcleanedHighHtgNN" or selection=="DILEPcleanedHighHtgFinal":
				bins = ["binhighEMHT_24","binhighEMHT_25","binhighEMHT_26"]
			
			for nBin in bins:
				newYield = outputDC.exp[nBin]["signal"]*ratio
				
				outputDC.newSignalYield({nBin: newYield},"signal")
			
			directory="output/HTG/"+selection
			if not os.path.exists(directory):
				os.makedirs(directory)
			
			outputDC.write(filename=directory+"/"+f)
			
def Scan_NeutralinoBrDiphoton(selection):
		
		for f in os.listdir("input/Diphoton"):
			path="input/Diphoton/"+f
			outputDC = MyDatacard(path)
			print f
			mass=getMasses(f)
			oldXsec=getXsec(mass,"xSec_SMS_N2C1_13TeV.pkl")[0]
			newXsec=getXsec(mass,"xSec_SMS_C1C1_13TeV.pkl")[0]
			ratio=newXsec/oldXsec
			
			bins = ["bin1","bin2","bin3","bin4","bin5","bin6"]

			for nBin in bins:
				newYield = outputDC.exp[nBin]["t5gg"]*ratio
				scale = outputDC.getUncertaintyGamma("mcStats_"+str(int(nBin[3])-1),nBin,"t5gg")*ratio
				
				outputDC.newSignalDiphoton({nBin: newYield},{"mcStats_"+str(int(nBin[3])-1): {nBin: scale}})
			
			directory="output/Diphoton"
			if not os.path.exists(directory):
				os.makedirs(directory)
			
			outputDC.write(filename=directory+"/"+f)
			
def Scan_NeutralinoBrLepton(selection):
		
		for f in os.listdir("input/Lepton"):
			path="input/Lepton/"+f
			outputDC = MyDatacard(path)
			print f
			mass=getMasses(f)
			oldXsec=getXsec(mass,"xSec_SMS_N2C1_13TeV.pkl")[0]
			newXsec=getXsec(mass,"xSec_SMS_C1C1_13TeV.pkl")[0]
			ratio=newXsec/oldXsec

			for n in range(1,37,1):
				nBin="bin"+str(n)
				newYield = outputDC.exp[nBin]["SUSY"]*ratio
				
				outputDC.newSignalYield({nBin: newYield},"SUSY")
			
			directory="output/Lepton"
			if not os.path.exists(directory):
				os.makedirs(directory)
			
			outputDC.write(filename=directory+"/"+f)
	
#~ selection="inclusiveNN"
#~ selection="inclusiveFinal"
#~ selection="exclusiveHighHtgNN"
selection="exclusiveHighHtgFinal"
#~ selection="originalNN"
#~ selection="originalFinal"
#~ selection="DILEPcleanedHighHtgNN"
#~ selection="DILEPcleanedHighHtgFinal"

Scan_NeutralinoBr(selection)
#~ Scan_NeutralinoBrHTG(selection)
#~ Scan_NeutralinoBrDiphoton(selection)
#~ Scan_NeutralinoBrLepton(selection)
