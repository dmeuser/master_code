#Script to create datacards for different BRs

from MyDatacard import MyDatacard
import numpy as np
import argparse
import os


def Scan_NeutralinoBr():

	for mass in range(300,1325,25):
		DC_gg = MyDatacard("input/"+selection+"/TChiNg_gg/datacard_TChiNg_gg_"+str(mass)+".txt")
		DC_gz = MyDatacard("input/"+selection+"/TChiNg_gz/datacard_TChiNg_gz_"+str(mass)+".txt")
		DC_zz = MyDatacard("input/"+selection+"/TChiNg_zz/datacard_TChiNg_zz_"+str(mass)+".txt")
		outputDC = MyDatacard("input/"+selection+"/TChiNg_gg/datacard_TChiNg_gg_"+str(mass)+".txt")
		
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
				ggUpSyst = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty("syst_sig",nBin,"sig")
				gzUpSyst = DC_gz.exp[nBin]["sig"]+DC_gz.getUncertainty("syst_sig",nBin,"sig")
				zzUpSyst = DC_zz.exp[nBin]["sig"]+DC_zz.getUncertainty("syst_sig",nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSyst+2*x_BR*(1-x_BR)*gzUpSyst+(1-x_BR)**2*zzUpSyst
				
				ggDownSyst = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty("syst_sig",nBin,"sig")
				gzDownSyst = DC_gz.exp[nBin]["sig"]-DC_gz.getUncertainty("syst_sig",nBin,"sig")
				zzDownSyst = DC_zz.exp[nBin]["sig"]-DC_zz.getUncertainty("syst_sig",nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSyst+2*x_BR*(1-x_BR)*gzDownSyst+(1-x_BR)**2*zzDownSyst
				
				if newYield==0:
					systUnc = 1.0
				else:
					systUnc = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				#Systematic ISR uncertainty	
				ggUpSystISR = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty("ISRsyst_sig",nBin,"sig")
				gzUpSystISR = DC_gz.exp[nBin]["sig"]+DC_gz.getUncertainty("ISRsyst_sig",nBin,"sig")
				zzUpSystISR = DC_zz.exp[nBin]["sig"]+DC_zz.getUncertainty("ISRsyst_sig",nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSystISR+2*x_BR*(1-x_BR)*gzUpSystISR+(1-x_BR)**2*zzUpSystISR
				
				ggDownSystISR = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty("ISRsyst_sig",nBin,"sig")
				gzDownSystISR = DC_gz.exp[nBin]["sig"]-DC_gz.getUncertainty("ISRsyst_sig",nBin,"sig")
				zzDownSystISR = DC_zz.exp[nBin]["sig"]-DC_zz.getUncertainty("ISRsyst_sig",nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSystISR+2*x_BR*(1-x_BR)*gzDownSystISR+(1-x_BR)**2*zzDownSystISR
				
				if newYield==0:
					systUncISR = 1.0
				else:
					systUncISR = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				outputDC.newSignal({nBin: newYield},{"statB"+nBin.split("n")[1]+"sig": {nBin: round(statUnc,2)},
				"ISRsyst_sig": {nBin: systUncISR}, "syst_sig": {nBin: systUnc}})
			
			directory="output/NeutralinoBR/"+selection
			if not os.path.exists(directory):
				os.makedirs(directory)
			
			outputDC.write(filename=directory+"/datacard_TChiNg_BR_"+str(mass)+"_"+"%i"%(x_BR*100)+".txt")

def Scan_CharginoBr():

	for mass in range(300,1325,25):
		DC_gg = MyDatacard("input/"+selection+"/TChiNg_gg_C1N2/datacard_TChiNg_gg_C1N2_"+str(mass)+".txt")
		DC_wg = MyDatacard("input/"+selection+"/TChiWg/datacard_TChiWG_"+str(mass)+".txt")
		outputDC = MyDatacard("input/"+selection+"/TChiNg_gg_C1N2/datacard_TChiNg_gg_C1N2_"+str(mass)+".txt")
		
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
				ggUpSyst = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty("syst_sig",nBin,"sig")
				wgUpSyst = DC_wg.exp[nBin]["sig"]+DC_wg.getUncertainty("syst_sig",nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSyst+2*x_BR*(1-x_BR)*wgUpSyst
				
				ggDownSyst = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty("syst_sig",nBin,"sig")
				wgDownSyst = DC_wg.exp[nBin]["sig"]-DC_wg.getUncertainty("syst_sig",nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSyst+2*x_BR*(1-x_BR)*wgDownSyst
				
				if newYield==0:
					systUnc = 1.0
				else:
					systUnc = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				#Systematic ISR uncertainty	
				ggUpSystISR = DC_gg.exp[nBin]["sig"]+DC_gg.getUncertainty("ISRsyst_sig",nBin,"sig")
				wgUpSystISR = DC_wg.exp[nBin]["sig"]+DC_wg.getUncertainty("ISRsyst_sig",nBin,"sig")
				yieldUpSyst = x_BR**2*ggUpSystISR+2*x_BR*(1-x_BR)*wgUpSystISR
				
				ggDownSystISR = DC_gg.exp[nBin]["sig"]-DC_gg.getUncertainty("ISRsyst_sig",nBin,"sig")
				wgDownSystISR = DC_wg.exp[nBin]["sig"]-DC_wg.getUncertainty("ISRsyst_sig",nBin,"sig")
				yieldDownSyst = x_BR**2*ggDownSystISR+2*x_BR*(1-x_BR)*wgDownSystISR
				
				if newYield==0:
					systUncISR = 1.0
				else:
					systUncISR = 1+(yieldUpSyst-yieldDownSyst)/(2.*newYield)
				
				outputDC.newSignal({nBin: newYield},{"statB"+nBin.split("n")[1]+"sig": {nBin: round(statUnc,2)},
				"ISRsyst_sig": {nBin: systUncISR}, "syst_sig": {nBin: systUnc}})
			
			directory="output/CharginoBR/"+selection
			if not os.path.exists(directory):
				os.makedirs(directory)
				
			outputDC.write(filename=directory+"/datacard_CharginoBR_"+str(mass)+"_"+"%i"%(x_BR*100)+".txt")
			

#Run different scanning types:

parser = argparse.ArgumentParser()
parser.add_argument('scan', nargs='?', help="NeutralinoBR or CharginoBR")
parser.add_argument('selection', nargs='?', help="choose as selection like leptonVeto, htgVeto etc.")
args = parser.parse_args()

selection = args.selection

if args.scan=="NeutralinoBR":
	Scan_NeutralinoBr()
elif args.scan=="CharginoBR":
	Scan_CharginoBr()
else:
	print "Unknown scan"

