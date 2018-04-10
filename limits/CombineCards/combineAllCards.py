#Script to combine multiple datacards at once

import os
import sys

def multi_combine(path1,path2,output,ignore=[]):
	masses=[]
	
	if path2.find("lepton")!=-1:
		masses1=[]
		masses2=[]
		for f in os.listdir("./"+path1):
			masses1.append(f[14:-4])	
		for f in os.listdir("./"+path2):
			masses2.append(f[14:-4])
		for mass in masses1:
			if mass in masses2:
				masses.append(mass)
				
	else:
		for f in os.listdir("./"+path1):
			masses.append(f[14:-4])
			
	for mass in masses:
		n1=path1+"/datacard_T5Wg_"+mass+".txt"
		
		if path2.find("knut")!=-1:
			n2=path2+"/Wg_"+mass+".txt"
		elif path2.find("lepton")!=-1:
			n2=path2+"/counting_t5Wg_"+mass+".txt"
			
		command="combineCards.py"
		for binName in ignore:
			command=command+" --xc="+binName
		
		if path2.find("knut")!=-1:
			command=command+" Photon_ST="+n1+" Photon_HTG="+n2+" >"+output+"/datacard_T5Wg_"+mass+".txt"
		elif path2.find("lepton")!=-1:
			command=command+" Photon_ST="+n1+" Photon_Lepton="+n2+" >"+output+"/datacard_T5Wg_"+mass+".txt"
		
		if os.system(command)==0:
			print "T5Wg_"+mass+" succesfull"
		else:
			sys.exit("T5Wg_"+mass+" not succesfull")
	return 0
	

#~ multi_combine("input/T5Wg_htgVeto","input/T5Wg_knut","output/T5Wg_htgVeto_KnutFull")
#~ multi_combine("input/T5Wg_htgHighVeto","input/T5Wg_knut","output/T5Wg_htgHighVeto_KnutHighHtg",["binlowEMHT_24","binlowEMHT_25","binlowEMHT_26"])
multi_combine("input/T5Wg_leptonVeto","input/T5Wg_lepton","output/T5Wg_leptonVeto_leptonFull")
