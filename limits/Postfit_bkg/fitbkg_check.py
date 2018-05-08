import ROOT as rt

f1=rt.TFile("fitDiagnostics_plain.root","read")
f2=rt.TFile("fitDiagnostics_VgGjet.root","read")
f3=rt.TFile("fitDiagnostics_allCORR.root","read")

post_plain=f1.Get("norm_fit_b")
pre_plain=f1.Get("norm_prefit")

post_VgGjet=f2.Get("norm_fit_b")
pre_VgGjet=f2.Get("norm_prefit")

post_all=f3.Get("norm_fit_b")
pre_all=f3.Get("norm_prefit")

norm_ST={}
norm_Lepton={}

for bkg in ["GJ","TTcomb","Vg","diboson","efake"]:
	norm_ST[bkg]=[[],[],[]]
	for i in range(1,5):
		temp1=post_plain.find("ch1_Photon_ST_bin"+str(i)+"/"+bkg)
		temp2=pre_plain.find("ch1_Photon_ST_bin"+str(i)+"/"+bkg)
		pull=(temp1.getVal()-temp2.getVal())/temp2.getError()
		norm_ST[bkg][0].append(pull)
		
		temp1=post_VgGjet.find("ch1_Photon_ST_bin"+str(i)+"/"+bkg)
		temp2=pre_VgGjet.find("ch1_Photon_ST_bin"+str(i)+"/"+bkg)
		pull=(temp1.getVal()-temp2.getVal())/temp2.getError()
		norm_ST[bkg][1].append(pull)
		
		temp1=post_all.find("ch1_Photon_ST_bin"+str(i)+"/"+bkg)
		temp2=pre_all.find("ch1_Photon_ST_bin"+str(i)+"/"+bkg)
		pull=(temp1.getVal()-temp2.getVal())/temp2.getError()
		norm_ST[bkg][2].append(pull)


for bkg in ["VGamma","elefakepho","jetfakepho","qcdfakelep","rare"]:
	norm_Lepton[bkg]=[[],[],[]]
	for i in range(1,37):
		temp1=post_plain.find("ch1_Photon_Lepton_bin"+str(i)+"/"+bkg)
		temp2=pre_plain.find("ch1_Photon_Lepton_bin"+str(i)+"/"+bkg)
		pull=(temp1.getVal()-temp2.getVal())/temp2.getError()
		norm_Lepton[bkg][0].append(pull)
		
		temp1=post_VgGjet.find("ch1_Photon_Lepton_bin"+str(i)+"/"+bkg)
		temp2=pre_VgGjet.find("ch1_Photon_Lepton_bin"+str(i)+"/"+bkg)
		pull=(temp1.getVal()-temp2.getVal())/temp2.getError()
		norm_Lepton[bkg][1].append(pull)
		
		temp1=post_all.find("ch1_Photon_Lepton_bin"+str(i)+"/"+bkg)
		temp2=pre_all.find("ch1_Photon_Lepton_bin"+str(i)+"/"+bkg)
		pull=(temp1.getVal()-temp2.getVal())/temp2.getError()
		norm_Lepton[bkg][2].append(pull)

corrName=["noCorr","VgGjet","allCorr"]
f=rt.TFile("hists_prefit_bkg.root","update")

for key in norm_ST:
	i=0
	for corr in norm_ST[key]:
		hist=rt.TH1F("",";binNo;pull",4,0,4)
		for j in range(4):
			hist.SetBinContent(j+1,corr[j])
		hist.Write("ST_"+key+"_"+corrName[i],rt.TObject.kOverwrite)
		i+=1
		
for key in norm_Lepton:
	i=0
	for corr in norm_Lepton[key]:
		hist=rt.TH1F("",";binNo;pull",36,0,36)
		for j in range(36):
			hist.SetBinContent(j+1,corr[j])
		hist.Write("Lepton_"+key+"_"+corrName[i],rt.TObject.kOverwrite)
		i+=1

f.Close()
		

