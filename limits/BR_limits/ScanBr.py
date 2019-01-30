#Script to create datacards for different BRs

from MyDatacard import MyDatacard
import numpy as np
import argparse
import os
import sys

def getMasses(path,gluino_mass,setGluinoMass):
	masses=[]
	for f in os.listdir("./"+path):
		mass=[]
		splitted=f.split("_")
		for part in splitted:
			part=part.split(".")[0]
			if part.isdigit():
				mass.append(part)
			
		if setGluinoMass:	
			if len(mass)==2 and int(mass[0])==gluino_mass:
				masses.append(mass[1])
		else:
			if len(mass)==2 and int(mass[1])==gluino_mass:
				masses.append(mass[0])
	return masses

def Scan_NeutralinoBr():

	for mass in range(300,1325,25):
		DC_gg = MyDatacard("input/ST/"+selection+"/TChiNg_gg/datacard_TChiNg_gg_"+str(mass)+".txt")
		DC_gz = MyDatacard("input/ST/"+selection+"/TChiNg_gz/datacard_TChiNg_gz_"+str(mass)+".txt")
		DC_zz = MyDatacard("input/ST/"+selection+"/TChiNg_zz/datacard_TChiNg_zz_"+str(mass)+".txt")
		outputDC = MyDatacard("input/ST/"+selection+"/TChiNg_gg/datacard_TChiNg_gg_"+str(mass)+".txt")
		
		print mass
		
		for x_BR in range(0,102,2):
			x_BR = x_BR/100.
			
			for nBin in ["bin1","bin2","bin3","bin4"]:
				newYield = x_BR**2*DC_gg.exp[nBin]["sig"]+2*x_BR*(1-x_BR)*DC_gz.exp[nBin]["sig"]+(1-x_BR)**2*DC_zz.exp[nBin]["sig"]
				
				#Statistic uncertainty
				statUnc = np.sqrt((x_BR**2*DC_gg.getStatUncertainty(nBin,"sig"))**2+(2*x_BR*(1-x_BR)*DC_gz.getStatUncertainty(nBin,"sig"))**2+((1-x_BR)**2*DC_zz.getStatUncertainty(nBin,"sig"))**2)
				if newYield==0:
					statUnc = 1.0
				else:
					statUnc = 1+statUnc/newYield
				
				#Systematic uncertainty
				syst_sig="syst_sig"
				if "Final" in selection: syst_sig="GenMet"
				ggUpSyst = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty(syst_sig,nBin,"sig")
				gzUpSyst = DC_gz.exp[nBin]["sig"]+DC_gz.getUncertainty(syst_sig,nBin,"sig")
				zzUpSyst = DC_zz.exp[nBin]["sig"]+DC_zz.getUncertainty(syst_sig,nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSyst+2*x_BR*(1-x_BR)*gzUpSyst+(1-x_BR)**2*zzUpSyst
				
				ggDownSyst = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty(syst_sig,nBin,"sig")
				gzDownSyst = DC_gz.exp[nBin]["sig"]-DC_gz.getUncertainty(syst_sig,nBin,"sig")
				zzDownSyst = DC_zz.exp[nBin]["sig"]-DC_zz.getUncertainty(syst_sig,nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSyst+2*x_BR*(1-x_BR)*gzDownSyst+(1-x_BR)**2*zzDownSyst
				
				if newYield==0:
					systUnc = 1.0
				else:
					systUnc = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				#Systematic ISR uncertainty	
				isr_sig="ISRsyst_sig"
				if "Final" in selection: isr_sig="ISR"
				ggUpSystISR = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty(isr_sig,nBin,"sig")
				gzUpSystISR = DC_gz.exp[nBin]["sig"]+DC_gz.getUncertainty(isr_sig,nBin,"sig")
				zzUpSystISR = DC_zz.exp[nBin]["sig"]+DC_zz.getUncertainty(isr_sig,nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSystISR+2*x_BR*(1-x_BR)*gzUpSystISR+(1-x_BR)**2*zzUpSystISR
				
				ggDownSystISR = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty(isr_sig,nBin,"sig")
				gzDownSystISR = DC_gz.exp[nBin]["sig"]-DC_gz.getUncertainty(isr_sig,nBin,"sig")
				zzDownSystISR = DC_zz.exp[nBin]["sig"]-DC_zz.getUncertainty(isr_sig,nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSystISR+2*x_BR*(1-x_BR)*gzDownSystISR+(1-x_BR)**2*zzDownSystISR
				
				if newYield==0:
					systUncISR = 1.0
				else:
					systUncISR = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				if "Final" not in selection:
					outputDC.newSignal({nBin: newYield},{"statB"+nBin.split("n")[1]+"sig": {nBin: round(statUnc,2)},
					isr_sig: {nBin: systUncISR}, syst_sig: {nBin: systUnc}})
					continue
				
				#Systematic PU uncertainty	
				PU_sig="PU"
				ggUpSystPU = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty(PU_sig,nBin,"sig")
				gzUpSystPU = DC_gz.exp[nBin]["sig"]+DC_gz.getUncertainty(PU_sig,nBin,"sig")
				zzUpSystPU = DC_zz.exp[nBin]["sig"]+DC_zz.getUncertainty(PU_sig,nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSystPU+2*x_BR*(1-x_BR)*gzUpSystPU+(1-x_BR)**2*zzUpSystPU
				
				ggDownSystPU = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty(PU_sig,nBin,"sig")
				gzDownSystPU = DC_gz.exp[nBin]["sig"]-DC_gz.getUncertainty(PU_sig,nBin,"sig")
				zzDownSystPU = DC_zz.exp[nBin]["sig"]-DC_zz.getUncertainty(PU_sig,nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSystPU+2*x_BR*(1-x_BR)*gzDownSystPU+(1-x_BR)**2*zzDownSystPU
				
				if newYield==0:
					systUncPU = 1.0
				else:
					systUncPU = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				#Systematic JES uncertainty	
				JES_sig="JES"
				ggUpSystJES = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty(JES_sig,nBin,"sig")
				gzUpSystJES = DC_gz.exp[nBin]["sig"]+DC_gz.getUncertainty(JES_sig,nBin,"sig")
				zzUpSystJES = DC_zz.exp[nBin]["sig"]+DC_zz.getUncertainty(JES_sig,nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSystJES+2*x_BR*(1-x_BR)*gzUpSystJES+(1-x_BR)**2*zzUpSystJES
				
				ggDownSystJES = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty(JES_sig,nBin,"sig")
				gzDownSystJES = DC_gz.exp[nBin]["sig"]-DC_gz.getUncertainty(JES_sig,nBin,"sig")
				zzDownSystJES = DC_zz.exp[nBin]["sig"]-DC_zz.getUncertainty(JES_sig,nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSystJES+2*x_BR*(1-x_BR)*gzDownSystJES+(1-x_BR)**2*zzDownSystJES
				
				if newYield==0:
					systUncJES = 1.0
				else:
					systUncJES = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				#Systematic PhotonSF uncertainty	
				PhotonSF_sig="PhotonSF"
				ggUpSystPhotonSF = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty(PhotonSF_sig,nBin,"sig")
				gzUpSystPhotonSF = DC_gz.exp[nBin]["sig"]+DC_gz.getUncertainty(PhotonSF_sig,nBin,"sig")
				zzUpSystPhotonSF = DC_zz.exp[nBin]["sig"]+DC_zz.getUncertainty(PhotonSF_sig,nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSystPhotonSF+2*x_BR*(1-x_BR)*gzUpSystPhotonSF+(1-x_BR)**2*zzUpSystPhotonSF
				
				ggDownSystPhotonSF = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty(PhotonSF_sig,nBin,"sig")
				gzDownSystPhotonSF = DC_gz.exp[nBin]["sig"]-DC_gz.getUncertainty(PhotonSF_sig,nBin,"sig")
				zzDownSystPhotonSF = DC_zz.exp[nBin]["sig"]-DC_zz.getUncertainty(PhotonSF_sig,nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSystPhotonSF+2*x_BR*(1-x_BR)*gzDownSystPhotonSF+(1-x_BR)**2*zzDownSystPhotonSF
				
				if newYield==0:
					systUncPhotonSF = 1.0
				else:
					systUncPhotonSF = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				outputDC.newSignal({nBin: newYield},{"statB"+nBin.split("n")[1]+"sig": {nBin: round(statUnc,2)},
				isr_sig: {nBin: systUncISR}, syst_sig: {nBin: systUnc}, PU_sig: {nBin: systUncPU}, JES_sig: {nBin: systUncJES}, PhotonSF_sig: {nBin: systUncPhotonSF}})
			
			directory="output/ST/NeutralinoBR/"+selection
			if not os.path.exists(directory):
				os.makedirs(directory)
			
			outputDC.write(filename=directory+"/datacard_TChiNg_BR_"+str(mass)+"_"+"%i"%(x_BR*100)+".txt")

def Scan_CharginoBr():

	for mass in range(300,1325,25):
		print mass
		DC_gg = MyDatacard("input/ST/"+selection+"/TChiNg_gg_C1N2/datacard_TChiNg_gg_C1N2_"+str(mass)+".txt")
		DC_wg = MyDatacard("input/ST/"+selection+"/TChiWg/datacard_TChiWG_"+str(mass)+".txt")
		outputDC = MyDatacard("input/ST/"+selection+"/TChiNg_gg_C1N2/datacard_TChiNg_gg_C1N2_"+str(mass)+".txt")
		
		for x_BR in range(0,102,2):
			x_BR = x_BR/100.
			
			for nBin in ["bin1","bin2","bin3","bin4"]:
				newYield = x_BR**2*DC_gg.exp[nBin]["sig"]+2*x_BR*(1-x_BR)*DC_wg.exp[nBin]["sig"]
				
				#Statistic uncertainty
				statUnc = np.sqrt((x_BR**2*DC_gg.getStatUncertainty(nBin,"sig"))**2+(2*x_BR*(1-x_BR)*DC_wg.getStatUncertainty(nBin,"sig"))**2)
				if newYield==0:
					statUnc = 1.0
				else:
					statUnc = 1+statUnc/newYield
				
				#Systematic uncertainty
				syst_sig="syst_sig"
				if "Final" in selection: syst_sig="GenMet"
				ggUpSyst = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty(syst_sig,nBin,"sig")
				wgUpSyst = DC_wg.exp[nBin]["sig"]+DC_wg.getUncertainty(syst_sig,nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSyst+2*x_BR*(1-x_BR)*wgUpSyst
				
				ggDownSyst = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty(syst_sig,nBin,"sig")
				wgDownSyst = DC_wg.exp[nBin]["sig"]-DC_wg.getUncertainty(syst_sig,nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSyst+2*x_BR*(1-x_BR)*wgDownSyst
				
				if newYield==0:
					systUnc = 1.0
				else:
					systUnc = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				#Systematic ISR uncertainty	
				isr_sig="ISRsyst_sig"
				if "Final" in selection: isr_sig="ISR"
				ggUpSystISR = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty(isr_sig,nBin,"sig")
				wgUpSystISR = DC_wg.exp[nBin]["sig"]+DC_wg.getUncertainty(isr_sig,nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSystISR+2*x_BR*(1-x_BR)*wgUpSystISR
				
				ggDownSystISR = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty(isr_sig,nBin,"sig")
				wgDownSystISR = DC_wg.exp[nBin]["sig"]-DC_wg.getUncertainty(isr_sig,nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSystISR+2*x_BR*(1-x_BR)*wgDownSystISR
				
				if newYield==0:
					systUncISR = 1.0
				else:
					systUncISR = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				if "Final" not in selection:
					outputDC.newSignal({nBin: newYield},{"statB"+nBin.split("n")[1]+"sig": {nBin: round(statUnc,2)},
					isr_sig: {nBin: systUncISR}, syst_sig: {nBin: systUnc}})
					continue
				
				#Systematic PU uncertainty	
				PU_sig="PU"
				ggUpSystPU = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty(PU_sig,nBin,"sig")
				wgUpSystPU = DC_wg.exp[nBin]["sig"]+DC_wg.getUncertainty(PU_sig,nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSystPU+2*x_BR*(1-x_BR)*wgUpSystPU
				
				ggDownSystPU = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty(PU_sig,nBin,"sig")
				wgDownSystPU = DC_wg.exp[nBin]["sig"]-DC_wg.getUncertainty(PU_sig,nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSystPU+2*x_BR*(1-x_BR)*wgDownSystPU
				
				if newYield==0:
					systUncPU = 1.0
				else:
					systUncPU = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				#Systematic JES uncertainty	
				JES_sig="JES"
				ggUpSystJES = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty(JES_sig,nBin,"sig")
				wgUpSystJES = DC_wg.exp[nBin]["sig"]+DC_wg.getUncertainty(JES_sig,nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSystJES+2*x_BR*(1-x_BR)*wgUpSystJES
				
				ggDownSystJES = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty(JES_sig,nBin,"sig")
				wgDownSystJES = DC_wg.exp[nBin]["sig"]-DC_wg.getUncertainty(JES_sig,nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSystJES+2*x_BR*(1-x_BR)*wgDownSystJES
				
				if newYield==0:
					systUncJES = 1.0
				else:
					systUncJES = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				#Systematic PhotonSF uncertainty	
				PhotonSF_sig="PhotonSF"
				ggUpSystPhotonSF = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty(PhotonSF_sig,nBin,"sig")
				wgUpSystPhotonSF = DC_wg.exp[nBin]["sig"]+DC_wg.getUncertainty(PhotonSF_sig,nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSystPhotonSF+2*x_BR*(1-x_BR)*wgUpSystPhotonSF
				
				ggDownSystPhotonSF = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty(PhotonSF_sig,nBin,"sig")
				wgDownSystPhotonSF = DC_wg.exp[nBin]["sig"]-DC_wg.getUncertainty(PhotonSF_sig,nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSystPhotonSF+2*x_BR*(1-x_BR)*wgDownSystPhotonSF
				
				if newYield==0:
					systUncPhotonSF = 1.0
				else:
					systUncPhotonSF = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				outputDC.newSignal({nBin: newYield},{"statB"+nBin.split("n")[1]+"sig": {nBin: round(statUnc,2)},
				isr_sig: {nBin: systUncISR}, syst_sig: {nBin: systUnc}, PU_sig: {nBin: systUncPU}, JES_sig: {nBin: systUncJES}, PhotonSF_sig: {nBin: systUncPhotonSF}})
			
			directory="output/ST/CharginoBR/"+selection
			if not os.path.exists(directory):
				os.makedirs(directory)
				
			if x_BR==0.58:
				outputDC.write(filename=directory+"/datacard_CharginoBR_"+str(mass)+"_58"+".txt")
			else:
				outputDC.write(filename=directory+"/datacard_CharginoBR_"+str(mass)+"_"+"%i"%(x_BR*100)+".txt")
	print "Remember to rescale to C1_C1 xsec!"
			
def Scan_CharginoBr_strong(gluino_mass,neutralino_mass):
	
	setGluinomass=True
	setmass=gluino_mass
	if neutralino_mass!=0:
		setGluinomass=False
		setmass=neutralino_mass
	
	masses=getMasses("input/ST/"+selection+"/T5gg/",setmass,setGluinomass)

	for mass in masses:
		if setGluinomass:
			DC_gg = MyDatacard("input/ST/"+selection+"/T5gg/datacard_T5gg_"+str(setmass)+"_"+str(mass)+".txt")
			DC_wg = MyDatacard("input/ST/"+selection+"/T5Wg/datacard_T5Wg_"+str(setmass)+"_"+str(mass)+".txt")
			outputDC = MyDatacard("input/ST/"+selection+"/T5gg/datacard_T5gg_"+str(setmass)+"_"+str(mass)+".txt")
		else:
			DC_gg = MyDatacard("input/ST/"+selection+"/T5gg/datacard_T5gg_"+str(mass)+"_"+str(setmass)+".txt")
			DC_wg = MyDatacard("input/ST/"+selection+"/T5Wg/datacard_T5Wg_"+str(mass)+"_"+str(setmass)+".txt")
			outputDC = MyDatacard("input/ST/"+selection+"/T5gg/datacard_T5gg_"+str(mass)+"_"+str(setmass)+".txt")
		
		print mass
		
		for x_BR in range(0,102,2):
			x_BR = x_BR/100.
			
			for nBin in ["bin1","bin2","bin3","bin4"]:
				newYield = x_BR**2*DC_gg.exp[nBin]["sig"]+2*x_BR*(1-x_BR)*DC_wg.exp[nBin]["sig"]
				
				#Statistic uncertainty
				statUnc = np.sqrt((x_BR**2*DC_gg.getStatUncertainty(nBin,"sig"))**2+(2*x_BR*(1-x_BR)*DC_wg.getStatUncertainty(nBin,"sig"))**2)
				if newYield==0:
					statUnc = 1.0
				else:
					statUnc = 1+statUnc/newYield
				#Systematic uncertainty
				syst_sig="syst_sig"
				if "Final" in selection: syst_sig="GenMet"
				ggUpSyst = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty(syst_sig,nBin,"sig")
				wgUpSyst = DC_wg.exp[nBin]["sig"]+DC_wg.getUncertainty(syst_sig,nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSyst+2*x_BR*(1-x_BR)*wgUpSyst
				
				ggDownSyst = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty(syst_sig,nBin,"sig")
				wgDownSyst = DC_wg.exp[nBin]["sig"]-DC_wg.getUncertainty(syst_sig,nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSyst+2*x_BR*(1-x_BR)*wgDownSyst
				
				if newYield==0:
					systUnc = 1.0
				else:
					systUnc = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				#Systematic ISR uncertainty	
				isr_sig="ISRsyst_sig"
				if "Final" in selection: isr_sig="ISR"
				ggUpSystISR = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty(isr_sig,nBin,"sig")
				wgUpSystISR = DC_wg.exp[nBin]["sig"]+DC_wg.getUncertainty(isr_sig,nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSystISR+2*x_BR*(1-x_BR)*wgUpSystISR
				
				ggDownSystISR = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty(isr_sig,nBin,"sig")
				wgDownSystISR = DC_wg.exp[nBin]["sig"]-DC_wg.getUncertainty(isr_sig,nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSystISR+2*x_BR*(1-x_BR)*wgDownSystISR
				
				if newYield==0:
					systUncISR = 1.0
				else:
					systUncISR = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				if "Final" not in selection:
					outputDC.newSignal({nBin: newYield},{"statB"+nBin.split("n")[1]+"sig": {nBin: round(statUnc,2)},
					isr_sig: {nBin: systUncISR}, syst_sig: {nBin: systUnc}})
					continue
				
				#Systematic PU uncertainty	
				PU_sig="PU"
				ggUpSystPU = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty(PU_sig,nBin,"sig")
				wgUpSystPU = DC_wg.exp[nBin]["sig"]+DC_wg.getUncertainty(PU_sig,nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSystPU+2*x_BR*(1-x_BR)*wgUpSystPU
				
				ggDownSystPU = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty(PU_sig,nBin,"sig")
				wgDownSystPU = DC_wg.exp[nBin]["sig"]-DC_wg.getUncertainty(PU_sig,nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSystPU+2*x_BR*(1-x_BR)*wgDownSystPU
				
				if newYield==0:
					systUncPU = 1.0
				else:
					systUncPU = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				#Systematic JES uncertainty	
				JES_sig="JES"
				ggUpSystJES = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty(JES_sig,nBin,"sig")
				wgUpSystJES = DC_wg.exp[nBin]["sig"]+DC_wg.getUncertainty(JES_sig,nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSystJES+2*x_BR*(1-x_BR)*wgUpSystJES
				
				ggDownSystJES = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty(JES_sig,nBin,"sig")
				wgDownSystJES = DC_wg.exp[nBin]["sig"]-DC_wg.getUncertainty(JES_sig,nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSystJES+2*x_BR*(1-x_BR)*wgDownSystJES
				
				if newYield==0:
					systUncJES = 1.0
				else:
					systUncJES = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				#Systematic PhotonSF uncertainty	
				PhotonSF_sig="PhotonSF"
				ggUpSystPhotonSF = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty(PhotonSF_sig,nBin,"sig")
				wgUpSystPhotonSF = DC_wg.exp[nBin]["sig"]+DC_wg.getUncertainty(PhotonSF_sig,nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSystPhotonSF+2*x_BR*(1-x_BR)*wgUpSystPhotonSF
				
				ggDownSystPhotonSF = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty(PhotonSF_sig,nBin,"sig")
				wgDownSystPhotonSF = DC_wg.exp[nBin]["sig"]-DC_wg.getUncertainty(PhotonSF_sig,nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSystPhotonSF+2*x_BR*(1-x_BR)*wgDownSystPhotonSF
				
				if newYield==0:
					systUncPhotonSF = 1.0
				else:
					systUncPhotonSF = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				outputDC.newSignal({nBin: newYield},{"statB"+nBin.split("n")[1]+"sig": {nBin: round(statUnc,2)},
				isr_sig: {nBin: systUncISR}, syst_sig: {nBin: systUnc}, PU_sig: {nBin: systUncPU}, JES_sig: {nBin: systUncJES}, PhotonSF_sig: {nBin: systUncPhotonSF}})
			
			if setGluinomass:
				directory="output/ST/CharginoBRstrongG"+str(setmass)+"/"+selection
				fname=directory+"/datacard_CharginoBRstrongG"+str(setmass)+"_"+str(mass)+"_"+"%i"%(x_BR*100)+".txt"
			else:
				directory="output/ST/CharginoBRstrongN"+str(setmass)+"/"+selection
				fname=directory+"/datacard_CharginoBRstrongN"+str(setmass)+"_"+str(mass)+"_"+"%i"%(x_BR*100)+".txt"
			if not os.path.exists(directory):
				os.makedirs(directory)
				
			outputDC.write(filename=fname)
			
def Scan_NeutralinoBr_HTG():

	for mass in range(300,1325,25):
		DC_gg = MyDatacard("input/HTG/"+selection+"/TChiNg_gg/"+str(mass)+"_0.txt")
		DC_gz = MyDatacard("input/HTG/"+selection+"/TChiNg_gz/"+str(mass)+"_0.txt")
		DC_zz = MyDatacard("input/HTG/"+selection+"/TChiNg_zz/"+str(mass)+"_0.txt")
		outputDC = MyDatacard("input/HTG/"+selection+"/TChiNg_gg/"+str(mass)+"_0.txt")
		print mass
		
		for x_BR in range(0,102,2):
			x_BR = x_BR/100.
			
			bins = ["binlowEMHT_24","binlowEMHT_25","binlowEMHT_26","binhighEMHT_24","binhighEMHT_25","binhighEMHT_26"]
			if selection=="highHTG" or selection=="DILEPcleanedHighHtgNN" or selection=="DILEPcleanedHighHtgFinal" or selection=="DILEPcleanedHighHtgFinalPre":
				bins = ["binhighEMHT_24","binhighEMHT_25","binhighEMHT_26"]
			
			for nBin in bins:
				newYield = x_BR**2*DC_gg.exp[nBin]["signal"]+2*x_BR*(1-x_BR)*DC_gz.exp[nBin]["signal"]+(1-x_BR)**2*DC_zz.exp[nBin]["signal"]
				
				
				#Statistic uncertainty
				statUnc = np.sqrt((x_BR**2*DC_gg.getUncertainty("signalStat_"+nBin,nBin,"signal"))**2+(2*x_BR*(1-x_BR)*DC_gz.getUncertainty("signalStat_"+nBin,nBin,"signal"))**2+((1-x_BR)**2*DC_zz.getUncertainty("signalStat_"+nBin,nBin,"signal"))**2)
				if newYield==0:
					statUnc = 1.0
				else:
					statUnc = 1+statUnc/newYield
				
				#Systematic uncertainty genMet
				genMet_sig="genMet"
				if "Final" in selection: genMet_sig="GenMet"
				ggUpSyst = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty(genMet_sig,nBin,"signal")
				gzUpSyst = DC_gz.exp[nBin]["signal"]+DC_gz.getUncertainty(genMet_sig,nBin,"signal")
				zzUpSyst = DC_zz.exp[nBin]["signal"]+DC_zz.getUncertainty(genMet_sig,nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSyst+2*x_BR*(1-x_BR)*gzUpSyst+(1-x_BR)**2*zzUpSyst
				
				ggDownSyst = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty(genMet_sig,nBin,"signal")
				gzDownSyst = DC_gz.exp[nBin]["signal"]-DC_gz.getUncertainty(genMet_sig,nBin,"signal")
				zzDownSyst = DC_zz.exp[nBin]["signal"]-DC_zz.getUncertainty(genMet_sig,nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSyst+2*x_BR*(1-x_BR)*gzDownSyst+(1-x_BR)**2*zzDownSyst
				
				if newYield==0:
					systUnc = 1.0
				else:
					systUnc = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				#Systematic uncertainty	jes
				jes_sig="jes"
				if "Final" in selection: jes_sig="JES"
				ggUpSystJES = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty(jes_sig,nBin,"signal")
				gzUpSystJES = DC_gz.exp[nBin]["signal"]+DC_gz.getUncertainty(jes_sig,nBin,"signal")
				zzUpSystJES = DC_zz.exp[nBin]["signal"]+DC_zz.getUncertainty(jes_sig,nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSystJES+2*x_BR*(1-x_BR)*gzUpSystJES+(1-x_BR)**2*zzUpSystJES
				
				ggDownSystJES = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty(jes_sig,nBin,"signal")
				gzDownSystJES = DC_gz.exp[nBin]["signal"]-DC_gz.getUncertainty(jes_sig,nBin,"signal")
				zzDownSystJES = DC_zz.exp[nBin]["signal"]-DC_zz.getUncertainty(jes_sig,nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSystJES+2*x_BR*(1-x_BR)*gzDownSystJES+(1-x_BR)**2*zzDownSystJES
				
				if newYield==0:
					systUncJES = 1.0
				else:
					systUncJES = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				#Systematic uncertainty	pu
				pu_sig="pu"
				if "Final" in selection: pu_sig="PU"
				ggUpSystPU = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty(pu_sig,nBin,"signal")
				gzUpSystPU = DC_gz.exp[nBin]["signal"]+DC_gz.getUncertainty(pu_sig,nBin,"signal")
				zzUpSystPU = DC_zz.exp[nBin]["signal"]+DC_zz.getUncertainty(pu_sig,nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSystPU+2*x_BR*(1-x_BR)*gzUpSystPU+(1-x_BR)**2*zzUpSystPU
				
				ggDownSystPU = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty(pu_sig,nBin,"signal")
				gzDownSystPU = DC_gz.exp[nBin]["signal"]-DC_gz.getUncertainty(pu_sig,nBin,"signal")
				zzDownSystPU = DC_zz.exp[nBin]["signal"]-DC_zz.getUncertainty(pu_sig,nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSystPU+2*x_BR*(1-x_BR)*gzDownSystPU+(1-x_BR)**2*zzDownSystPU
				
				if newYield==0:
					systUncPU = 1.0
				else:
					systUncPU = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				#Systematic uncertainty	scale
				ggUpSystSCALE = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty("scale",nBin,"signal")
				gzUpSystSCALE = DC_gz.exp[nBin]["signal"]+DC_gz.getUncertainty("scale",nBin,"signal")
				zzUpSystSCALE = DC_zz.exp[nBin]["signal"]+DC_zz.getUncertainty("scale",nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSystSCALE+2*x_BR*(1-x_BR)*gzUpSystSCALE+(1-x_BR)**2*zzUpSystSCALE
				
				ggDownSystSCALE = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty("scale",nBin,"signal")
				gzDownSystSCALE = DC_gz.exp[nBin]["signal"]-DC_gz.getUncertainty("scale",nBin,"signal")
				zzDownSystSCALE = DC_zz.exp[nBin]["signal"]-DC_zz.getUncertainty("scale",nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSystSCALE+2*x_BR*(1-x_BR)*gzDownSystSCALE+(1-x_BR)**2*zzDownSystSCALE
				
				if newYield==0:
					systUncSCALE = 1.0
				else:
					systUncSCALE = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				if "Final" not in selection:
					outputDC.newSignalHTG({nBin: newYield},{"signalStat_"+nBin: {nBin: round(statUnc,2)},
					"jes": {nBin: systUncJES}, "genMet": {nBin: systUnc}, "pu": {nBin: systUncPU}, "scale": {nBin: systUncSCALE}})
					continue
				
				#Systematic uncertainty	isr
				isr_sig="ISR"
				ggUpSystISR = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty(isr_sig,nBin,"signal")
				gzUpSystISR = DC_gz.exp[nBin]["signal"]+DC_gz.getUncertainty(isr_sig,nBin,"signal")
				zzUpSystISR = DC_zz.exp[nBin]["signal"]+DC_zz.getUncertainty(isr_sig,nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSystISR+2*x_BR*(1-x_BR)*gzUpSystISR+(1-x_BR)**2*zzUpSystISR
				
				ggDownSystISR = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty(isr_sig,nBin,"signal")
				gzDownSystISR = DC_gz.exp[nBin]["signal"]-DC_gz.getUncertainty(isr_sig,nBin,"signal")
				zzDownSystISR = DC_zz.exp[nBin]["signal"]-DC_zz.getUncertainty(isr_sig,nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSystISR+2*x_BR*(1-x_BR)*gzDownSystISR+(1-x_BR)**2*zzDownSystISR
				
				if newYield==0:
					systUncISR = 1.0
				else:
					systUncISR = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				#Systematic uncertainty	isr
				PhotonSF_sig="PhotonSF"
				ggUpSystPhotonSF = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty(PhotonSF_sig,nBin,"signal")
				gzUpSystPhotonSF = DC_gz.exp[nBin]["signal"]+DC_gz.getUncertainty(PhotonSF_sig,nBin,"signal")
				zzUpSystPhotonSF = DC_zz.exp[nBin]["signal"]+DC_zz.getUncertainty(PhotonSF_sig,nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSystPhotonSF+2*x_BR*(1-x_BR)*gzUpSystPhotonSF+(1-x_BR)**2*zzUpSystPhotonSF
				
				ggDownSystPhotonSF = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty(PhotonSF_sig,nBin,"signal")
				gzDownSystPhotonSF = DC_gz.exp[nBin]["signal"]-DC_gz.getUncertainty(PhotonSF_sig,nBin,"signal")
				zzDownSystPhotonSF = DC_zz.exp[nBin]["signal"]-DC_zz.getUncertainty(PhotonSF_sig,nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSystPhotonSF+2*x_BR*(1-x_BR)*gzDownSystPhotonSF+(1-x_BR)**2*zzDownSystPhotonSF
				
				if newYield==0:
					systUncPhotonSF = 1.0
				else:
					systUncPhotonSF = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				outputDC.newSignalHTG({nBin: newYield},{"signalStat_"+nBin: {nBin: round(statUnc,2)},
					jes_sig: {nBin: systUncJES}, genMet_sig: {nBin: systUnc}, pu_sig: {nBin: systUncPU}, "scale": {nBin: systUncSCALE}, isr_sig: {nBin: systUncISR}, PhotonSF_sig: {nBin: systUncPhotonSF}})
			
			directory="output/HTG/NeutralinoBR/"+selection
			if not os.path.exists(directory):
				os.makedirs(directory)
			
			outputDC.write(filename=directory+"/datacard_TChiNg_BR_"+str(mass)+"_"+"%i"%(x_BR*100)+".txt")

def Scan_CharginoBr_HTG():

	for mass in range(300,1325,25):
		DC_gg = MyDatacard("input/HTG/"+selection+"/TChiNg_gg_C1N2/"+str(mass)+"_0.txt")
		DC_wg = MyDatacard("input/HTG/"+selection+"/TChiWg/"+str(mass)+"_0.txt")
		outputDC = MyDatacard("input/HTG/"+selection+"/TChiNg_gg_C1N2/"+str(mass)+"_0.txt")
		print mass
		
		for x_BR in range(0,102,2):
			x_BR = x_BR/100.
			
			bins = ["binlowEMHT_24","binlowEMHT_25","binlowEMHT_26","binhighEMHT_24","binhighEMHT_25","binhighEMHT_26"]
			if selection=="highHTG" or selection=="DILEPcleanedHighHtgNN" or selection=="DILEPcleanedHighHtgFinal" or selection=="DILEPcleanedHighHtgFinalPre":
				bins = ["binhighEMHT_24","binhighEMHT_25","binhighEMHT_26"]
			
			for nBin in bins:
				newYield = x_BR**2*DC_gg.exp[nBin]["signal"]+2*x_BR*(1-x_BR)*DC_wg.exp[nBin]["signal"]
				
				
				#Statistic uncertainty
				statUnc = np.sqrt((x_BR**2*DC_gg.getUncertainty("signalStat_"+nBin,nBin,"signal"))**2+(2*x_BR*(1-x_BR)*DC_wg.getUncertainty("signalStat_"+nBin,nBin,"signal"))**2)
				if newYield==0:
					statUnc = 1.0
				else:
					statUnc = 1+statUnc/newYield
				
				#Systematic uncertainty genMet
				genMet_sig="genMet"
				if "Final" in selection: genMet_sig="GenMet"
				ggUpSyst = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty(genMet_sig,nBin,"signal")
				wgUpSyst = DC_wg.exp[nBin]["signal"]+DC_wg.getUncertainty(genMet_sig,nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSyst+2*x_BR*(1-x_BR)*wgUpSyst
				
				ggDownSyst = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty(genMet_sig,nBin,"signal")
				wgDownSyst = DC_wg.exp[nBin]["signal"]-DC_wg.getUncertainty(genMet_sig,nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSyst+2*x_BR*(1-x_BR)*wgDownSyst
				
				if newYield==0:
					systUnc = 1.0
				else:
					systUnc = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				#Systematic uncertainty	jes
				jes_sig="jes"
				if "Final" in selection: jes_sig="JES"
				ggUpSystJES = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty(jes_sig,nBin,"signal")
				wgUpSystJES = DC_wg.exp[nBin]["signal"]+DC_wg.getUncertainty(jes_sig,nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSystJES+2*x_BR*(1-x_BR)*wgUpSystJES
				
				ggDownSystJES = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty(jes_sig,nBin,"signal")
				wgDownSystJES = DC_wg.exp[nBin]["signal"]-DC_wg.getUncertainty(jes_sig,nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSystJES+2*x_BR*(1-x_BR)*wgDownSystJES
				
				if newYield==0:
					systUncJES = 1.0
				else:
					systUncJES = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				#Systematic uncertainty	pu
				pu_sig="pu"
				if "Final" in selection: pu_sig="PU"
				ggUpSystPU = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty(pu_sig,nBin,"signal")
				wgUpSystPU = DC_wg.exp[nBin]["signal"]+DC_wg.getUncertainty(pu_sig,nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSystPU+2*x_BR*(1-x_BR)*wgUpSystPU
				
				ggDownSystPU = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty(pu_sig,nBin,"signal")
				wgDownSystPU = DC_wg.exp[nBin]["signal"]-DC_wg.getUncertainty(pu_sig,nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSystPU+2*x_BR*(1-x_BR)*wgDownSystPU
				
				if newYield==0:
					systUncPU = 1.0
				else:
					systUncPU = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				#Systematic uncertainty	scale
				ggUpSystSCALE = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty("scale",nBin,"signal")
				wgUpSystSCALE = DC_wg.exp[nBin]["signal"]+DC_wg.getUncertainty("scale",nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSystSCALE+2*x_BR*(1-x_BR)*wgUpSystSCALE
				
				ggDownSystSCALE = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty("scale",nBin,"signal")
				wgDownSystSCALE = DC_wg.exp[nBin]["signal"]-DC_wg.getUncertainty("scale",nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSystSCALE+2*x_BR*(1-x_BR)*wgDownSystSCALE
				
				if newYield==0:
					systUncSCALE = 1.0
				else:
					systUncSCALE = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				if "Final" not in selection:
					outputDC.newSignalHTG({nBin: newYield},{"signalStat_"+nBin: {nBin: round(statUnc,2)},
					"jes": {nBin: systUncJES}, "genMet": {nBin: systUnc}, "pu": {nBin: systUncPU}, "scale": {nBin: systUncSCALE}})
					continue
				
				#Systematic uncertainty	isr
				isr_sig="ISR"
				ggUpSystISR = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty(isr_sig,nBin,"signal")
				wgUpSystISR = DC_wg.exp[nBin]["signal"]+DC_wg.getUncertainty(isr_sig,nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSystISR+2*x_BR*(1-x_BR)*wgUpSystISR
				
				ggDownSystISR = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty(isr_sig,nBin,"signal")
				wgDownSystISR = DC_wg.exp[nBin]["signal"]-DC_wg.getUncertainty(isr_sig,nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSystISR+2*x_BR*(1-x_BR)*wgDownSystISR
				
				if newYield==0:
					systUncISR = 1.0
				else:
					systUncISR = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				#Systematic uncertainty	isr
				PhotonSF_sig="PhotonSF"
				ggUpSystPhotonSF = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty(PhotonSF_sig,nBin,"signal")
				wgUpSystPhotonSF = DC_wg.exp[nBin]["signal"]+DC_wg.getUncertainty(PhotonSF_sig,nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSystPhotonSF+2*x_BR*(1-x_BR)*wgUpSystPhotonSF
				
				ggDownSystPhotonSF = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty(PhotonSF_sig,nBin,"signal")
				wgDownSystPhotonSF = DC_wg.exp[nBin]["signal"]-DC_wg.getUncertainty(PhotonSF_sig,nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSystPhotonSF+2*x_BR*(1-x_BR)*wgDownSystPhotonSF
				
				if newYield==0:
					systUncPhotonSF = 1.0
				else:
					systUncPhotonSF = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				outputDC.newSignalHTG({nBin: newYield},{"signalStat_"+nBin: {nBin: round(statUnc,2)},
					jes_sig: {nBin: systUncJES}, genMet_sig: {nBin: systUnc}, pu_sig: {nBin: systUncPU}, "scale": {nBin: systUncSCALE}, isr_sig: {nBin: systUncISR}, PhotonSF_sig: {nBin: systUncPhotonSF}})
			
			directory="output/HTG/CharginoBR/"+selection
			if not os.path.exists(directory):
				os.makedirs(directory)
			
			if x_BR==0.58:
				outputDC.write(filename=directory+"/datacard_CharginoBR_"+str(mass)+"_58"+".txt")
			else:
				outputDC.write(filename=directory+"/datacard_CharginoBR_"+str(mass)+"_"+"%i"%(x_BR*100)+".txt")
	print "Remember to rescale to C1_C1 xsec!"

def Scan_CharginoBr_strong_HTG(gluino_mass,neutralino_mass):
	
	setGluinomass=True
	setmass=gluino_mass
	if neutralino_mass!=0:
		setGluinomass=False
		setmass=neutralino_mass
	
	masses=getMasses("input/HTG/"+selection+"/T5gg/",setmass,setGluinomass)
	
	for mass in masses:
		if setGluinomass:
			DC_gg = MyDatacard("input/HTG/"+selection+"/T5gg/"+str(setmass)+"_"+str(mass)+".txt")
			DC_wg = MyDatacard("input/HTG/"+selection+"/T5Wg/"+str(setmass)+"_"+str(mass)+".txt")
			outputDC = MyDatacard("input/HTG/"+selection+"/T5gg/"+str(setmass)+"_"+str(mass)+".txt")
		else:
			DC_gg = MyDatacard("input/HTG/"+selection+"/T5gg/"+str(mass)+"_"+str(setmass)+".txt")
			DC_wg = MyDatacard("input/HTG/"+selection+"/T5Wg/"+str(mass)+"_"+str(setmass)+".txt")
			outputDC = MyDatacard("input/HTG/"+selection+"/T5gg/"+str(mass)+"_"+str(setmass)+".txt")
		print mass
		
		for x_BR in range(0,102,2):
			x_BR = x_BR/100.
			
			bins = ["binlowEMHT_24","binlowEMHT_25","binlowEMHT_26","binhighEMHT_24","binhighEMHT_25","binhighEMHT_26"]
			if selection=="highHTG" or selection=="LEPcleanedHighHtgNN" or selection=="LEPcleanedHighHtgFinal" or selection=="LEPcleanedHighHtgFinalPre":
				bins = ["binhighEMHT_24","binhighEMHT_25","binhighEMHT_26"]
			
			for nBin in bins:
				newYield = x_BR**2*DC_gg.exp[nBin]["signal"]+2*x_BR*(1-x_BR)*DC_wg.exp[nBin]["signal"]
				
				
				#Statistic uncertainty
				statUnc = np.sqrt((x_BR**2*DC_gg.getUncertainty("signalStat_"+nBin,nBin,"signal"))**2+(2*x_BR*(1-x_BR)*DC_wg.getUncertainty("signalStat_"+nBin,nBin,"signal"))**2)
				if newYield==0:
					statUnc = 1.0
				else:
					statUnc = 1+statUnc/newYield
				
				#Systematic uncertainty genMet
				genMet_sig="genMet"
				if "Final" in selection: genMet_sig="GenMet"
				ggUpSyst = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty(genMet_sig,nBin,"signal")
				wgUpSyst = DC_wg.exp[nBin]["signal"]+DC_wg.getUncertainty(genMet_sig,nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSyst+2*x_BR*(1-x_BR)*wgUpSyst
				
				ggDownSyst = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty(genMet_sig,nBin,"signal")
				wgDownSyst = DC_wg.exp[nBin]["signal"]-DC_wg.getUncertainty(genMet_sig,nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSyst+2*x_BR*(1-x_BR)*wgDownSyst
				
				if newYield==0:
					systUnc = 1.0
				else:
					systUnc = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				#Systematic uncertainty	jes
				jes_sig="jes"
				if "Final" in selection: jes_sig="JES"
				ggUpSystJES = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty(jes_sig,nBin,"signal")
				wgUpSystJES = DC_wg.exp[nBin]["signal"]+DC_wg.getUncertainty(jes_sig,nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSystJES+2*x_BR*(1-x_BR)*wgUpSystJES
				
				ggDownSystJES = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty(jes_sig,nBin,"signal")
				wgDownSystJES = DC_wg.exp[nBin]["signal"]-DC_wg.getUncertainty(jes_sig,nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSystJES+2*x_BR*(1-x_BR)*wgDownSystJES
				
				if newYield==0:
					systUncJES = 1.0
				else:
					systUncJES = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				#Systematic uncertainty	pu
				pu_sig="pu"
				if "Final" in selection: pu_sig="PU"
				ggUpSystPU = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty(pu_sig,nBin,"signal")
				wgUpSystPU = DC_wg.exp[nBin]["signal"]+DC_wg.getUncertainty(pu_sig,nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSystPU+2*x_BR*(1-x_BR)*wgUpSystPU
				
				ggDownSystPU = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty(pu_sig,nBin,"signal")
				wgDownSystPU = DC_wg.exp[nBin]["signal"]-DC_wg.getUncertainty(pu_sig,nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSystPU+2*x_BR*(1-x_BR)*wgDownSystPU
				
				if newYield==0:
					systUncPU = 1.0
				else:
					systUncPU = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				#Systematic uncertainty	scale
				ggUpSystSCALE = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty("scale",nBin,"signal")
				wgUpSystSCALE = DC_wg.exp[nBin]["signal"]+DC_wg.getUncertainty("scale",nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSystSCALE+2*x_BR*(1-x_BR)*wgUpSystSCALE
				
				ggDownSystSCALE = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty("scale",nBin,"signal")
				wgDownSystSCALE = DC_wg.exp[nBin]["signal"]-DC_wg.getUncertainty("scale",nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSystSCALE+2*x_BR*(1-x_BR)*wgDownSystSCALE
				
				if newYield==0:
					systUncSCALE = 1.0
				else:
					systUncSCALE = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				if "Final" not in selection:
					outputDC.newSignalHTG({nBin: newYield},{"signalStat_"+nBin: {nBin: round(statUnc,2)},
					"jes": {nBin: systUncJES}, "genMet": {nBin: systUnc}, "pu": {nBin: systUncPU}, "scale": {nBin: systUncSCALE}})
					continue
				
				#Systematic uncertainty	isr
				isr_sig="ISR"
				ggUpSystISR = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty(isr_sig,nBin,"signal")
				wgUpSystISR = DC_wg.exp[nBin]["signal"]+DC_wg.getUncertainty(isr_sig,nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSystISR+2*x_BR*(1-x_BR)*wgUpSystISR
				
				ggDownSystISR = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty(isr_sig,nBin,"signal")
				wgDownSystISR = DC_wg.exp[nBin]["signal"]-DC_wg.getUncertainty(isr_sig,nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSystISR+2*x_BR*(1-x_BR)*wgDownSystISR
				
				if newYield==0:
					systUncISR = 1.0
				else:
					systUncISR = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				#Systematic uncertainty	isr
				PhotonSF_sig="PhotonSF"
				ggUpSystPhotonSF = DC_gg.exp[nBin]["signal"]+DC_gg.getUncertainty(PhotonSF_sig,nBin,"signal")
				wgUpSystPhotonSF = DC_wg.exp[nBin]["signal"]+DC_wg.getUncertainty(PhotonSF_sig,nBin,"signal")
				yieldUpSyst = x_BR**2*ggUpSystPhotonSF+2*x_BR*(1-x_BR)*wgUpSystPhotonSF
				
				ggDownSystPhotonSF = DC_gg.exp[nBin]["signal"]-DC_gg.getUncertainty(PhotonSF_sig,nBin,"signal")
				wgDownSystPhotonSF = DC_wg.exp[nBin]["signal"]-DC_wg.getUncertainty(PhotonSF_sig,nBin,"signal")
				yieldDownSyst = x_BR**2*ggDownSystPhotonSF+2*x_BR*(1-x_BR)*wgDownSystPhotonSF
				
				if newYield==0:
					systUncPhotonSF = 1.0
				else:
					systUncPhotonSF = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				outputDC.newSignalHTG({nBin: newYield},{"signalStat_"+nBin: {nBin: round(statUnc,2)},
					jes_sig: {nBin: systUncJES}, genMet_sig: {nBin: systUnc}, pu_sig: {nBin: systUncPU}, "scale": {nBin: systUncSCALE}, isr_sig: {nBin: systUncISR}, PhotonSF_sig: {nBin: systUncPhotonSF}})
			
			if setGluinomass:
				directory="output/HTG/CharginoBRstrongG"+str(setmass)+"/"+selection
				fname=directory+"/datacard_CharginoBRstrongG"+str(setmass)+"_"+str(mass)+"_"+"%i"%(x_BR*100)+".txt"
			else:
				directory="output/HTG/CharginoBRstrongN"+str(setmass)+"/"+selection
				fname=directory+"/datacard_CharginoBRstrongN"+str(setmass)+"_"+str(mass)+"_"+"%i"%(x_BR*100)+".txt"
			if not os.path.exists(directory):
				os.makedirs(directory)
			
			outputDC.write(filename=fname)
			
def Scan_CharginoBr_strong_diphoton(gluino_mass,neutralino_mass):
	
	setGluinomass=True
	setmass=gluino_mass
	if neutralino_mass!=0:
		setGluinomass=False
		setmass=neutralino_mass
	
	masses=getMasses("input/Diphoton/"+selection+"/T5gg/",setmass,setGluinomass)
	
	for mass in masses:
		if setGluinomass:
			DC_gg = MyDatacard("input/Diphoton/"+selection+"/T5gg/datacard_T5gg_"+str(setmass)+"_"+str(mass)+".txt")
			outputDC = MyDatacard("input/Diphoton/"+selection+"/T5gg/datacard_T5gg_"+str(setmass)+"_"+str(mass)+".txt")
		else:
			DC_gg = MyDatacard("input/Diphoton/"+selection+"/T5gg/datacard_T5gg_"+str(setmass)+"_"+str(mass)+".txt")
			outputDC = MyDatacard("input/Diphoton/"+selection+"/T5gg/datacard_T5gg_"+str(setmass)+"_"+str(mass)+".txt")
		print mass
		
		for x_BR in range(0,102,2):
			x_BR = x_BR/100.
			
			bins = ["bin1","bin2","bin3","bin4","bin5","bin6"]
			
			for nBin in bins:
				newYield = x_BR**2*DC_gg.exp[nBin]["t5gg"]
				
				
				#Statistic uncertainty
				statUnc = x_BR**2*DC_gg.getUncertaintyGamma("mcStats_"+str(int(nBin[3])-1),nBin,"t5gg")
				if newYield==0:
					statUnc = 1.0
				
				#Systematic uncertainty photon scale factor
				ggUpSyst = DC_gg.exp[nBin]["t5gg"]+DC_gg.getUncertainty("phoSf",nBin,"t5gg")
				yieldUpSyst = x_BR**2*ggUpSyst
				
				ggDownSyst = DC_gg.exp[nBin]["t5gg"]-DC_gg.getUncertainty("phoSf",nBin,"t5gg")
				yieldDownSyst = x_BR**2*ggDownSyst
				
				if newYield==0:
					systUnc = 1.0
				else:
					systUnc = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				#Systematic uncertainty	jes
				ggUpSystJES = DC_gg.exp[nBin]["t5gg"]+DC_gg.getUncertainty("jes",nBin,"t5gg")
				yieldUpSyst = x_BR**2*ggUpSystJES
				
				ggDownSystJES = DC_gg.exp[nBin]["t5gg"]-DC_gg.getUncertainty("jes",nBin,"t5gg")
				yieldDownSyst = x_BR**2*ggDownSystJES
				
				if newYield==0:
					systUncJES = 1.0
				else:
					systUncJES = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				
				outputDC.newSignalDiphoton({nBin: newYield},{"mcStats_"+str(int(nBin[3])-1): {nBin: statUnc}, "phoSf": {nBin: systUnc},
				"jes": {nBin: systUncJES}})
			
			if setGluinomass:
				directory="output/Diphoton/CharginoBRstrongG"+str(setmass)+"/"+selection
				fname=directory+"/datacard_CharginoBRstrongG"+str(setmass)+"_"+str(mass)+"_"+"%i"%(x_BR*100)+".txt"
			else:
				directory="output/Diphoton/CharginoBRstrongN"+str(setmass)+"/"+selection
				fname=directory+"/datacard_CharginoBRstrongN"+str(setmass)+"_"+str(mass)+"_"+"%i"%(x_BR*100)+".txt"
			if not os.path.exists(directory):
				os.makedirs(directory)
			
			outputDC.write(filename=fname)
			
def Scan_CharginoBr_strong_lepton(gluino_mass,neutralino_mass):
	
	setGluinomass=True
	setmass=gluino_mass
	if neutralino_mass!=0:
		setGluinomass=False
		setmass=neutralino_mass
	
	masses=getMasses("input/Lepton/"+selection+"/T5gg/",setmass,setGluinomass)

	for mass in masses:
		if setGluinomass:
			DC_gg = MyDatacard("input/Lepton/"+selection+"/T5gg/datacard_CharginoBRstrong_GG_"+str(setmass)+"_"+str(mass)+".txt")
			DC_wg = MyDatacard("input/Lepton/"+selection+"/T5Wg/datacard_CharginoBRstrong_WG_"+str(setmass)+"_"+str(mass)+".txt")
			outputDC = MyDatacard("input/Lepton/"+selection+"/T5gg/datacard_CharginoBRstrong_GG_"+str(setmass)+"_"+str(mass)+".txt")
		else:
			DC_gg = MyDatacard("input/Lepton/"+selection+"/T5gg/datacard_CharginoBRstrong_GG_"+str(mass)+"_"+str(setmass)+".txt")
			DC_wg = MyDatacard("input/Lepton/"+selection+"/T5Wg/datacard_CharginoBRstrong_WG_"+str(mass)+"_"+str(setmass)+".txt")
			outputDC = MyDatacard("input/Lepton/"+selection+"/T5gg/datacard_CharginoBRstrong_GG_"+str(mass)+"_"+str(setmass)+".txt")
		
		print mass
		
		for x_BR in range(0,102,2):
			x_BR = x_BR/100.
			
			for n in range(1,37,1):
				nBin="bin"+str(n)
				newYield = x_BR**2*DC_gg.exp[nBin]["SUSY"]+2*x_BR*(1-x_BR)*DC_wg.exp[nBin]["SUSY"]
				
				#Statistic uncertainty
				statUnc = np.sqrt((x_BR**2*DC_gg.getUncertainty("SUSY_stat"+str(n),nBin,"SUSY"))**2+(2*x_BR*(1-x_BR)*DC_wg.getUncertainty("SUSY_stat"+str(n),nBin,"SUSY"))**2)
				if newYield==0:
					statUnc = 1.0
				else:
					statUnc = 1+statUnc/newYield
				
				#Systematic jes uncertainty
				ggUpSyst = DC_gg.exp[nBin]["SUSY"]+DC_gg.getUncertainty("JES",nBin,"SUSY")
				wgUpSyst = DC_wg.exp[nBin]["SUSY"]+DC_wg.getUncertainty("JES",nBin,"SUSY")
				yieldUpSyst = x_BR**2*ggUpSyst+2*x_BR*(1-x_BR)*wgUpSyst
				
				ggDownSyst = DC_gg.exp[nBin]["SUSY"]-DC_gg.getUncertainty("JES",nBin,"SUSY")
				wgDownSyst = DC_wg.exp[nBin]["SUSY"]-DC_wg.getUncertainty("JES",nBin,"SUSY")
				yieldDownSyst = x_BR**2*ggDownSyst+2*x_BR*(1-x_BR)*wgDownSyst
				
				if newYield==0:
					systUnc = 1.0
				else:
					systUnc = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				#Systematic ISR uncertainty	
				ggUpSystISR = DC_gg.exp[nBin]["SUSY"]+DC_gg.getUncertainty("ISR",nBin,"SUSY")
				wgUpSystISR = DC_wg.exp[nBin]["SUSY"]+DC_wg.getUncertainty("ISR",nBin,"SUSY")
				yieldUpSyst = x_BR**2*ggUpSystISR+2*x_BR*(1-x_BR)*wgUpSystISR
				
				ggDownSystISR = DC_gg.exp[nBin]["SUSY"]-DC_gg.getUncertainty("ISR",nBin,"SUSY")
				wgDownSystISR = DC_wg.exp[nBin]["SUSY"]-DC_wg.getUncertainty("ISR",nBin,"SUSY")
				yieldDownSyst = x_BR**2*ggDownSystISR+2*x_BR*(1-x_BR)*wgDownSystISR
				
				if newYield==0:
					systUncISR = 1.0
				else:
					systUncISR = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
					
				#Systematic ElectronSF uncertainty	
				ggUpSystESF = DC_gg.exp[nBin]["SUSY"]+DC_gg.getUncertainty("ElectronSF",nBin,"SUSY")
				wgUpSystESF = DC_wg.exp[nBin]["SUSY"]+DC_wg.getUncertainty("ElectronSF",nBin,"SUSY")
				yieldUpSyst = x_BR**2*ggUpSystESF+2*x_BR*(1-x_BR)*wgUpSystESF
				
				ggDownSystESF = DC_gg.exp[nBin]["SUSY"]-DC_gg.getUncertainty("ElectronSF",nBin,"SUSY")
				wgDownSystESF = DC_wg.exp[nBin]["SUSY"]-DC_wg.getUncertainty("ElectronSF",nBin,"SUSY")
				yieldDownSyst = x_BR**2*ggDownSystESF+2*x_BR*(1-x_BR)*wgDownSystESF
				
				if newYield==0:
					systUncESF = 1.0
				else:
					systUncESF = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				#Systematic PhotonSF uncertainty	
				ggUpSystPSF = DC_gg.exp[nBin]["SUSY"]+DC_gg.getUncertainty("PhotonSF",nBin,"SUSY")
				wgUpSystPSF = DC_wg.exp[nBin]["SUSY"]+DC_wg.getUncertainty("PhotonSF",nBin,"SUSY")
				yieldUpSyst = x_BR**2*ggUpSystPSF+2*x_BR*(1-x_BR)*wgUpSystPSF
				
				ggDownSystPSF = DC_gg.exp[nBin]["SUSY"]-DC_gg.getUncertainty("PhotonSF",nBin,"SUSY")
				wgDownSystPSF = DC_wg.exp[nBin]["SUSY"]-DC_wg.getUncertainty("PhotonSF",nBin,"SUSY")
				yieldDownSyst = x_BR**2*ggDownSystPSF+2*x_BR*(1-x_BR)*wgDownSystPSF
				
				if newYield==0:
					systUncPSF = 1.0
				else:
					systUncPSF = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				#Systematic GenMet uncertainty	
				ggUpSystGenMet = DC_gg.exp[nBin]["SUSY"]+DC_gg.getUncertainty("GenMet",nBin,"SUSY")
				wgUpSystGenMet = DC_wg.exp[nBin]["SUSY"]+DC_wg.getUncertainty("GenMet",nBin,"SUSY")
				yieldUpSyst = x_BR**2*ggUpSystGenMet+2*x_BR*(1-x_BR)*wgUpSystGenMet
				
				ggDownSystGenMet = DC_gg.exp[nBin]["SUSY"]-DC_gg.getUncertainty("GenMet",nBin,"SUSY")
				wgDownSystGenMet = DC_wg.exp[nBin]["SUSY"]-DC_wg.getUncertainty("GenMet",nBin,"SUSY")
				yieldDownSyst = x_BR**2*ggDownSystGenMet+2*x_BR*(1-x_BR)*wgDownSystGenMet
				
				if newYield==0:
					systUncGenMet = 1.0
				else:
					systUncGenMet = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				#Systematic PDFSCALE uncertainty	
				ggUpSystPDFSCALE = DC_gg.exp[nBin]["SUSY"]+DC_gg.getUncertainty("PDFSCALE",nBin,"SUSY")
				wgUpSystPDFSCALE = DC_wg.exp[nBin]["SUSY"]+DC_wg.getUncertainty("PDFSCALE",nBin,"SUSY")
				yieldUpSyst = x_BR**2*ggUpSystPDFSCALE+2*x_BR*(1-x_BR)*wgUpSystPDFSCALE
				
				ggDownSystPDFSCALE = DC_gg.exp[nBin]["SUSY"]-DC_gg.getUncertainty("PDFSCALE",nBin,"SUSY")
				wgDownSystPDFSCALE = DC_wg.exp[nBin]["SUSY"]-DC_wg.getUncertainty("PDFSCALE",nBin,"SUSY")
				yieldDownSyst = x_BR**2*ggDownSystPDFSCALE+2*x_BR*(1-x_BR)*wgDownSystPDFSCALE
				
				if newYield==0:
					systUncPDFSCALE = 1.0
				else:
					systUncPDFSCALE = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				#Systematic PU uncertainty	
				ggUpSystPU = DC_gg.exp[nBin]["SUSY"]+DC_gg.getUncertainty("PU",nBin,"SUSY")
				wgUpSystPU = DC_wg.exp[nBin]["SUSY"]+DC_wg.getUncertainty("PU",nBin,"SUSY")
				yieldUpSyst = x_BR**2*ggUpSystPU+2*x_BR*(1-x_BR)*wgUpSystPU
				
				ggDownSystPU = DC_gg.exp[nBin]["SUSY"]-DC_gg.getUncertainty("PU",nBin,"SUSY")
				wgDownSystPU = DC_wg.exp[nBin]["SUSY"]-DC_wg.getUncertainty("PU",nBin,"SUSY")
				yieldDownSyst = x_BR**2*ggDownSystPU+2*x_BR*(1-x_BR)*wgDownSystPU
				
				if newYield==0:
					systUncPU = 1.0
				else:
					systUncPU = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				outputDC.newSignalLepton({nBin: newYield},{"SUSY_stat"+str(n): {nBin: round(statUnc,2)},
					"JES": {nBin: systUnc}, "ElectronSF": {nBin: systUncESF}, "PhotonSF": {nBin: systUncPSF}, "ISR": {nBin: systUncISR}
					, "GenMet": {nBin: systUncGenMet}, "PDFSCALE": {nBin: systUncPDFSCALE}, "PU": {nBin: systUncPU}})
			
			if setGluinomass:
				directory="output/Lepton/CharginoBRstrongG"+str(setmass)+"/"+selection
				fname=directory+"/datacard_CharginoBRstrongG"+str(setmass)+"_"+str(mass)+"_"+"%i"%(x_BR*100)+".txt"
			else:
				directory="output/Lepton/CharginoBRstrongN"+str(setmass)+"/"+selection
				fname=directory+"/datacard_CharginoBRstrongN"+str(setmass)+"_"+str(mass)+"_"+"%i"%(x_BR*100)+".txt"
			if not os.path.exists(directory):
				os.makedirs(directory)
				
			outputDC.write(filename=fname)
			
#Run different scanning types:

parser = argparse.ArgumentParser()
parser.add_argument('analysis', nargs='?', help="ST or HTG or Diphoton")
parser.add_argument('scan', nargs='?', help="NeutralinoBR or CharginoBR or CharginoBRstrong")
parser.add_argument('selection', nargs='?', help="choose as selection like leptonVeto, htgVeto etc.")
parser.add_argument('--gluinomass', type=int, default=0, help="For CharginoBRstrong choose gluino mass or neutralino mass")
parser.add_argument('--neutralinomass', type=int, default=0, help="For CharginoBRstrong choose gluino mass or neutralino mass")
args = parser.parse_args()

selection = args.selection
sys.argv=[]


if (args.analysis=="ST"):
	if args.scan=="NeutralinoBR":
		Scan_NeutralinoBr()
	elif args.scan=="CharginoBR":
		Scan_CharginoBr()
	elif args.scan=="CharginoBRstrong":
		if args.neutralinomass==0 and args.gluinomass==0:
			print "Choose neutralino or gluinomass"
		else:
			Scan_CharginoBr_strong(args.gluinomass,args.neutralinomass)
	else:
		print "Unknown scan"

elif (args.analysis=="HTG"):
	if args.scan=="NeutralinoBR":
		Scan_NeutralinoBr_HTG()
	elif args.scan=="CharginoBR":
		Scan_CharginoBr_HTG()
	elif args.scan=="CharginoBRstrong":
		if args.neutralinomass==0 and args.gluinomass==0:
			print "Choose neutralino or gluinomass"
		else:
			Scan_CharginoBr_strong_HTG(args.gluinomass,args.neutralinomass)
	else:
		print "Unknown scan"
		
elif (args.analysis=="Diphoton"):
	#~ if args.scan=="NeutralinoBR":
		#~ Scan_NeutralinoBr_HTG()
	#~ elif args.scan=="CharginoBR":
		#~ Scan_CharginoBr_HTG()
	if args.scan=="CharginoBRstrong":
		if args.neutralinomass==0 and args.gluinomass==0:
			print "Choose neutralino or gluinomass"
		else:
			Scan_CharginoBr_strong_diphoton(args.gluinomass,args.neutralinomass)
	else:
		print "Unknown scan"

elif (args.analysis=="Lepton"):
	#~ if args.scan=="NeutralinoBR":
		#~ Scan_NeutralinoBr_HTG()
	#~ elif args.scan=="CharginoBR":
		#~ Scan_CharginoBr_HTG()
	if args.scan=="CharginoBRstrong":
		if args.neutralinomass==0 and args.gluinomass==0:
			print "Choose neutralino or gluinomass"
		else:
			Scan_CharginoBr_strong_lepton(args.gluinomass,args.neutralinomass)
	else:
		print "Unknown scan"

else: print "Unknown analysis"

