import ROOT as rt

outdir="output_old/"
selection="inclusiv"
sScan="GGM_M1_M2"

f=rt.TFile(outdir+"limits_%s_"%sScan+selection+".root","read")

gr=f.Get("gr_obs")

out=rt.TFile("test/contour.root","update")
count=gr.GetContourList(1.)

#~ hist=f.Get("h_obs")
hist=f.Get("h_exp")
hist.Write("hist",rt.TObject.kOverwrite)
hist.SetContour(1)
hist.SetContourLevel(0,1.)

c=rt.TCanvas()
c.cd()
hist.Draw("CONT1")
histSM=rt.TH2F(hist)
histSM.Smooth()
histSM.SetLineColor(rt.kRed)
histSM.Draw("CONT2 same")
c.Write("histCont",rt.TObject.kOverwrite)



for i in xrange(hist.GetNbinsX()+1):
	for j in xrange(hist.GetNbinsY()+1):
		if hist.GetBinContent(i,j)<1 :
			hist.SetBinContent(i,j,0)
		else :
			hist.SetBinContent(i,j,2)

c2=rt.TCanvas()
c2.cd()
hist.SetLineColor(rt.kGreen)
hist.Draw("colz")
hist.Draw("CONT2 same")
c2.Write("contBinaer",rt.TObject.kOverwrite)

c.cd()
hist.Draw("CONT Z LIST same")
c.Update()
conts=rt.gROOT.GetListOfSpecials().FindObject("contours")
contLevel=conts.At(0)

c.cd()
Nmax=0
i=0
for contGraph in contLevel:
	contGraph.Write("gr%i"%i,rt.TObject.kOverwrite)
	i+=1
	if contGraph.GetN()>Nmax:
		contGraphFinal=contGraph
		Nmax=contGraph.GetN()

hist.Draw("colz")
contGraphFinal.SetLineColor(rt.kRed)
contGraphFinal.Draw("same")
contGraphFinal.Write("contBinOnly",rt.TObject.kOverwrite)

#~ oldCont=f.Get("gr_obsC")
oldCont=f.Get("gr_expC")
oldCont.SetLineColor(rt.kBlue)
oldCont.Draw("same")
c.Write("contDiff",rt.TObject.kOverwrite)
	
out.Close()
f.Close()

