
from MyDatacard import MyDatacard
import numpy as np
import argparse
import os


def cancelSignal(datacard):

	DC = MyDatacard(datacard)
			
	for nBin in DC.bins:		
		DC.exp[nBin]["SUSY"]=0.0
			
	#~ print DC.systs[60][3]
	DC.write("out.txt")
			

cancelSignal("/home/home4/institut_1b/dmeuser/master_code/datacards/lepton/CharginoBR/datacard_CharginoBR_300_2.txt")

