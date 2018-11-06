from array import *

commonZmin = 4e-4
commonZmax = 2e-2

squarkZmin = commonZmin
squarkZmax = commonZmax
gluinoZmin = commonZmin
gluinoZmax = commonZmax

class sms():

    def __init__(self, modelname):
        if modelname.find("T1tttt") != -1: self.T1tttt()
        if modelname.find("T5ttttDM175") != -1: self.T5ttttDM175()
        if modelname.find("T1bbbb") != -1: self.T1bbbb()
        if modelname.find("T1qqqq") != -1: self.T1qqqq()
        if modelname.find("T5gg") != -1: self.T5gg()
        if modelname.find("T5Wg") != -1: self.T5Wg()
        if modelname.find("T6gg") != -1: self.T6gg()
        if modelname.find("T6Wg") != -1: self.T6Wg()
        if modelname.find("GGM1") != -1: self.GGM1()
        if modelname.find("GGM2") != -1: self.GGM2()
        if modelname.find("TChiNg") != -1: self.TChiNg_BR()
        if modelname.find("CharginoBR") != -1: self.CharginoBR()
        if modelname.find("CharginoBRstrong") != -1: self.CharginoBRstrong()

    def T6gg(self):
        # model name
        self.modelname = "T6gg"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma#tilde{G}"%(lsp_s,lsp_s)
        self.label2= "";
        # scan range to plot
        self.Xmin = 1100.
        self.Xmax = 2050.
        self.Ymin = 0.
        self.Ymax = 2700.
        self.Zmin = squarkZmin
        self.Zmax = squarkZmax
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV)"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} (GeV)"
        # diagonal lines
        self.diagOn = True

    def T6Wg(self):
        # model name
        self.modelname = "T6Wg"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}"%(lsp_s,lsp_s)
        self.label2= "";
        # scan range to plot
        self.Xmin = 1100.
        self.Xmax = 2050.
        self.Ymin = 0.
        self.Ymax = 2700.
        self.Zmin = squarkZmin
        self.Zmax = squarkZmax
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV)"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} (GeV)"
        # diagonal lines
        self.diagOn = True

    def T5gg(self):
        # model name
        self.modelname = "T5gg"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma#tilde{G}"%(lsp_s,lsp_s)
        self.label2= "";
        # scan range to plot
        self.Xmin = 1400.
        self.Xmax = 2500.
        self.Ymin = 0.
        self.Ymax = 3300.
        self.Zmin = gluinoZmin
        self.Zmax = gluinoZmax
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV)"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} (GeV)"
        # diagonal lines
        self.diagOn = True

    def T5Wg(self):
        # model name
        self.modelname = "T5Wg"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}"%(lsp_s,lsp_s)
        self.label2= "";
        # scan range to plot
        self.Xmin = 1400.
        self.Xmax = 2500.
        self.Ymin = 0.
        self.Ymax = 3300.
        self.Zmin = gluinoZmin
        self.Zmax = gluinoZmax
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV)"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} (GeV)"
        # diagonal lines
        self.diagOn = True

    #~ def GGM(self):
        #~ # model name
        #~ self.modelname = "GGM"
        #~ # decay chain
        #~ self.label= "GGM"
        #~ self.label2= "";
        #~ # scan range to plot
        #~ self.Xmin = 205-12.5
        #~ self.Xmax = 840+12.5
        #~ self.Ymin = 390-12.5
        #~ # self.Ymax = 1015+12.5
        #~ self.Ymax = 900
        #~ self.Zmin = 0.001
        #~ self.Zmax = 2.
        #~ # produce sparticle
        #~ self.sParticle = "m_{#tilde{B}} [GeV]"
        #~ # LSP
        #~ self.LSP = "m_{#tilde{W}} [GeV]"
        #~ # diagonal lines
        #~ self.diagOn = True

    def T1tttt(self):
        # model name
        self.modelname = "T1tttt"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow t #bar{t} "+lsp_s;
        self.label2= "";
        # scan range to plot
        self.Xmin = 600.
        self.Xmax = 2000.
        self.Ymin = 0.
        self.Ymax = 1800.
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} [GeV]"
        # turn off diagonal lines
        self.diagOn = False

    def T5ttttDM175(self):
        # model name
        self.modelname = "T5ttttDM175"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{g} #tilde{g},  #tilde{g} #rightarrow #tilde{t}_{1} t,  #tilde{t}_{1} #rightarrow #bar{t} "+lsp_s;
        self.label2= "m_{#tilde{t}_{1}} - m_{#tilde{#chi}^{0}_{1}} = 175 GeV";
        # scan range to plot
        self.Xmin = 600.
        self.Xmax = 1700.
        self.Ymin = 0.
        self.Ymax = 1800.
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} [GeV]"
        # turn off diagonal lines
        self.diagOn = False

    def T1bbbb(self):
        # model name
        self.modelname = "T1bbbb"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow b #bar{b} "+lsp_s;
        self.label2= "";
        # plot boundary. The top 1/4 of the y axis is taken by the legend
        self.Xmin = 600.
        self.Xmax = 1950.
        self.Ymin = 0.
        self.Ymax = 1800.
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} [GeV]"
        # turn off diagonal lines
        self.diagOn = False

    def T1qqqq(self):
        # model name
        self.modelname = "T1qqqq"
        # decay chain
        self.label= "pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow q #bar{q} #tilde{#chi}^{0}_{1}";
        self.label2= "";
        # plot boundary. The top 1/4 of the y axis is taken by the legend
        self.Xmin = 600.
        self.Xmax = 1950.
        self.Ymin = 0.
        self.Ymax = 1600.
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} [GeV]"
        # turn off diagonal lines
        self.diagOn = False

    def TChiNg_BR(self):
        # model name
        self.modelname = "TChiNg_BR"
        # decay chain
        self.label= "pp #rightarrow #lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"+" #lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} / "+"#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"+" #lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}, "+"#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow "+"#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}+soft, "+"#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow (#gamma,Z)#tilde{G}"
        self.label2= "";
        # plot boundary. The top 1/4 of the y axis is taken by the legend
        self.Xmin = 0.
        self.Xmax = 100.
        self.Ymin = 300.
        self.Ymax = 1600.
        self.Zmin = 0.01
        self.Zmax = 1000.
        # produce sparticle
        self.sParticle = "BF(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #gamma + #tilde{G}) (%)"
        # LSP
        self.LSP = "m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)"
        # turn off diagonal lines
        self.diagOn = False
        
    def CharginoBR(self):
        # model name
        self.modelname = "CharginoBR"
        # decay chain
        self.label= "pp #rightarrow #lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"+" #lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}, "+"#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow "+"#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}+soft/"+"W^{#pm}+#tilde{G}, "+"#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #gamma + #tilde{G}"
        self.label2= "";
        # plot boundary. The top 1/4 of the y axis is taken by the legend
        self.Xmin = 0.
        self.Xmax = 100.
        self.Ymin = 300.
        self.Ymax = 1600.
        self.Zmin = 0.01
        self.Zmax = 1000.
        # produce sparticle
        self.sParticle = "BF(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #tilde{#chi}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}(#gamma#tilde{G}) + soft) (%)";
        # LSP
        self.LSP = "m#kern[0.1]{_{#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}} (GeV)"
        # turn off diagonal lines
        self.diagOn = False

    def GGM1(self):
        # model name
        self.modelname = "GGM1"
        # decay chain
        self.label= "GGM M1M2"
        self.label2= "";
        # scan range to plot
        self.Xmin = 275
        self.Xmax = 1500
        self.Ymin = 275
        self.Ymax = 1800
        #~ self.Ymax = 900
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "#it{M}_{1} (GeV)"
        # LSP
        self.LSP = "#it{M}_{2} (GeV)"
        # diagonal lines
        self.diagOn = False
    
    def GGM2(self):
        # model name
        self.modelname = "GGM2"
        # decay chain
        self.label= "GGM M1M3"
        self.label2= "";
        # scan range to plot
        self.Xmin = 50
        self.Xmax = 1500
        self.Ymin = 1000
        self.Ymax = 2900
        #~ self.Ymax = 900
        self.Zmin = 0.001
        self.Zmax = 0.01
        # produce sparticle
        self.sParticle = "#it{M}_{1} (GeV)"
        # LSP
        self.LSP = "#it{M}_{3} (GeV)"
        # diagonal lines
        self.diagOn = False
       
    def CharginoBRstrong(self):
        # model name
        self.modelname = "CharginoBRstrong"
        # decay chain
        self.label= "pp #rightarrow #lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"+" #lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}, "+"#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow "+"#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}+soft/"+"W^{#pm}+#tilde{G}, "+"#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #gamma + #tilde{G}"
        self.label2= "";
        # plot boundary. The top 1/4 of the y axis is taken by the legend
        self.Xmin = 0.
        self.Xmax = 100.
        self.Ymin = 0.
        self.Ymax = 2400.
        self.Zmin = 0.000001
        self.Zmax = 10.
        # produce sparticle
        self.sParticle = "BF(#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}} #rightarrow #tilde{#chi}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}(#gamma#tilde{G}) + soft) (%)";
        # LSP
        self.LSP = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV)"
        # turn off diagonal lines
        self.diagOn = False
