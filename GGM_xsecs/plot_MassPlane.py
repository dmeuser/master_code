import ROOT as rt

rFile=rt.TFile("output/xsecs_M1_M2.root","read")

for temp in ["N1N2_MassPlane","N1C1_MassPlane"]:

	graph=rFile.Get(temp)

	cd=rt.TCanvas()
	cd.SetLogz()
	graph.GetHistogram().SetMinimum(0.0001)
	graph.Draw("hist colz")
	cd.SaveAs("output/"+temp+".pdf")
