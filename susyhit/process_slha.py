#script to produce susyhit and hdecay output for slha files and replace h/H BR from susyhit by hdecay

import os,sys


def susyhit(name,scan):
	path="scan/"+scan+"/"+name
	command_cp="cp "+path+" ./"
	command_mv="mv "+name+" slhaspectrum.in"
	command_susyhit="./run"
	command_cp2="cp susyhit_slha.out output/"+scan+"/susyhit/"
	command_mv2="mv output/"+scan+"/susyhit/susyhit_slha.out output/"+scan+"/susyhit/"+name
	if os.system(command_cp)==0 and os.system(command_mv)==0 and os.system(command_susyhit)==0 and os.system(command_cp2)==0 and os.system(command_mv2)==0:
		print "susyhit successfull for "+name
	else:
		sys.exit("susyhit for "+name+" not succesfull")
	return 0

def hdecay(name,scan):
	path="scan/"+scan+"/"+name
	path_hdecay="../hdecay/"
	command_cp="cp "+path+" ../hdecay/"
	command_mv="mv "+path_hdecay+name+" "+path_hdecay+"slha.in"
	command_hdecay="cd "+path_hdecay+"; ./run 2> /dev/null"
	command_cp2="cp "+path_hdecay+"slha.out output/"+scan+"/hdecay/"
	command_mv2="mv output/"+scan+"/hdecay/slha.out output/"+scan+"/hdecay/"+name
	if os.system(command_cp)==0 and os.system(command_mv)==0 and os.system(command_hdecay)==0 and os.system(command_cp2)==0 and os.system(command_mv2)==0:
		print "hdecay successfull for "+name
	else:
			sys.exit("hdecay for "+name+" not succesfull")
	return 0

def merge(name,scan):
	f=open("output/"+scan+"/hdecay/"+name,"r")
	content=f.readlines()
	start_h=0
	stop_h=0
	for i,j in enumerate(content):
		if j.find("DECAY        25")!=-1:
			start_h=i-2
		if j.find("DECAY        36")!=-1:
			stop_h=i-2
	hdecay=content[start_h:stop_h]
	f.close()
	
	f=open("output/"+scan+"/susyhit/"+name,"r")
	content=f.readlines()
	start_s=0
	stop_s=0
	for i,j in enumerate(content):
		if j.find("DECAY        25")!=-1:
			start_s=i-2
		if j.find("DECAY        36")!=-1:
			stop_s=i-2
	del content[start_s:stop_s]
	content[start_s:start_s]=hdecay
	f.close()
	
	f=open("output/"+scan+"/merged/"+name,"w")
	output=""
	for i in xrange(len(content)):
		output+=content[i]
	f.write(output)
	f.close
	
for scan in ["Bench_M1_M2","Bench_M1_M3"]:
	for f in os.listdir("scan/"+scan):
		susyhit(f,scan)
		hdecay(f,scan)
		merge(f,scan)

