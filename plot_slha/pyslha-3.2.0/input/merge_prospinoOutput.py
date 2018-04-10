import numpy as np
import os,sys

def merge(mass,scan):
	first = open("Xsec_"+scan+"_"+mass+".txt","r")
	result = open("Xsec_"+scan+"_"+mass+"_merged.txt","w")
	
	i=0
	
	for line in first:
		i+=1
		failed = False
		if i<=2 : continue
		
		add = open(mass+"_addition.txt","r")
		for line2 in add:
			masses = line2[0:9]
			if line[0:9]==masses:
				failed = True
				break
		
		if failed:
			result.write(line2)
			print "okay"
		else:
			result.write(line)
	
	add.close()
	first.close()
	result.close()
		

merge("M1_M2","Scan1")
merge("M1_M3","Scan2")
