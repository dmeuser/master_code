#!/usr/bin/env python2

import re
import subprocess
import ROOT
import optparse
import os
from Datacard import Datacard

class MyDatacard(Datacard):
    def __init__(self, dc=""):
        if not dc:
            self.bins = []
            self.obs = {}
            #self.processes = ['signal', 'gqcd', 'ele', 'zg', 'wg', 'ttg']
            self.processes = []
            self.signals = ['signal']
            #self.isSignal = {'wg': False, 'signal': True, 'zg': False, 'gqcd': False, 'ele': False, 'ttg': False}
            self.isSignal = {}
            self.keyline = []
            self.exp = {}
            #self.systs = [(x, False, 'lnN', [], {}) for x in "lumi", "jec", "pdf", "gqcdSyst", "eleSyst", "wgSyst", "zgSyst", "ttgSyst"]
            self.systs = []
            self.shapeMap = {}
            self.hasShape = False
            self.flatParamNuisances = {}
            self.rateParams = {}
            self.rateParamsOrder = []
        elif isinstance(dc, Datacard):
            self.bins = dc.bins
            self.obs = dc.obs
            self.processes = dc.processes
            self.signals = dc.signals
            self.isSignal = dc.isSignal
            self.keyline = dc.keyline
            self.exp = dc.exp
            self.systs = dc.systs
            self.shapeMap = dc.shapeMap
            self.hasShape = dc.hasShape
            self.flatParamNuisances = dc.flatParamNuisances
            self.rateParams = dc.rateParams
            self.rateParamsOrder = dc.rateParamsOrder
        elif dc.endswith(".txt"):
            import DatacardParser
            options, b = DatacardParser.addDatacardParserOptions(optparse.OptionParser())
            mydc = MyDatacard(DatacardParser.parseCard(file(dc), options))
            self.__init__(mydc)
        else:
            print "Do not know how to initialize MyDatacard with", dc


    def __str__(self):
        print "Class MyDatacard"
        print "bins               ", self.bins
        print "obs                ", self.obs
        print "processes          ", self.processes
        print "signals            ", self.signals
        print "isSignal           ", self.isSignal
        print "keyline            ", self.keyline
        print "exp                ", self.exp
        for s in self.systs:
            print "systs              ", s
        #print "shapeMap           ", self.shapeMap
        #print "hasShape           ", self.hasShape
        #print "flatParamNuicances ", self.flatParamNuisances
        #print "rateParams         ", self.rateParams
        #print "rateParamsOrder    ", self.rateParamsOrder
        return ""

    def _getProcessNumbers(self):
        counter = 1
        processNumbers = {}
        for a, b in self.isSignal.iteritems():
            if b: processNumbers[a] = 0
            else:
                processNumbers[a] = counter
                counter += 1
        return processNumbers


    def write(self, filename=""):
        
        
        maxInfoLen = max( [len(line[0])+7 for line in self.systs])
        out = ""
        out += "\nimax "+str(len(self.bins))
        out += "\njmax *"
        out += "\nkmax *"
        out += "\n\nbin         " + ("{:>40}"*len(self.bins)).format(*self.bins)
        out += "\nobservation " + ("{:>40}"*len(self.bins)).format(*[str(int(self.obs[x])) for x in self.bins])
        out += "\n\n"

        # create table for syst uncerts
        binNames, processNames, processNumbers = zip(*self.keyline)
        table = []
        table.append(["bin", ""] + list(binNames))
        table.append(["process", ""] + list(processNames))
        processNumbers = self._getProcessNumbers()
        table.append(["process", ""] + [str(processNumbers[x]) for x in processNames])
        table.append(["rate", ""] + [str(round(self.exp[bN][processNames[i]],3)) for i, bN in enumerate(binNames)])
        for line in self.systs:
            relUncerts = [line[4][bN][processNames[i]] for i, bN in enumerate(binNames)]
            templine = line[2]
            if len(line[3])==1:
                templine = line[2]+" "+str(line[3][0])
            table.append([line[0], templine] + ["-" if x==1 or x==0 else str(round(x,4)) for x in relUncerts])
        # format lengts of strings
        columnWidths = [max([len(i) for i in line])+1 for line in zip(*table)]
        for irow, row in enumerate(table):
            for icol, col in enumerate(row):
                table[irow][icol] = "{{:>{}}}".format(columnWidths[icol]+1).format(col)
        # append table to output
        for row in table: out += ''.join(row) + "\n"

        if filename:
            with open(filename, "wb") as f:
                f.write(out)
                #~ print "Writing to file:", filename
        else:
            print out

    def addBin(self, name, obs, bkgRates, bkgUncertainties):
        if self.processes:
            if self.processes != bkgRates.keys(): print "ERROR: Old processes", self.processes, " New processes:", bkgRates.keys()
        else:
            self.processes = bkgRates.keys()
            self.isSignal = dict([(r, r=="signal") for r in self.processes])
        self.bins.append(name)
        self.obs[name] = obs
        self.keyline.extend([(name, process, process=="signal") for process in bkgRates.keys()]) # TODD check order, take ordered dict???
        self.exp[name] = bkgRates

        systDict = dict([(l[0],l) for l in self.systs])
        for source, newUncerts in bkgUncertainties.iteritems():
            for p in self.processes:
                if p not in newUncerts:
                    newUncerts[p] = 0
            if source not in systDict:
                systDict[source] = (source, False, "lnN", [], dict([(b,dict([(r,0) for r in self.processes])) for b in self.bins]))
            systDict[source][4][name] = newUncerts
        for source, line in systDict.iteritems():
            for bin in self.bins:
                if bin not in line[4]:
                    systDict[source][4][bin] = dict([(r,0) for r in self.processes])
        self.systs = sorted(systDict.values())

    def newSignal(self, exp, unc):
        for bName, newRate in exp.iteritems():
            self.exp[bName]["sig"] = newRate
        systDict = dict([(l[0],l) for l in self.systs])
        for uncName, uncertaintyDict in unc.iteritems():
            for binName, u in uncertaintyDict.iteritems():
                systDict[uncName][4][binName]["sig"] = u
        self.systs = sorted(systDict.values())

    def newSignalHTG(self, exp, unc):
        for bName, newRate in exp.iteritems():
            self.exp[bName]["signal"] = newRate
        systDict = dict([(l[0],l) for l in self.systs])
        for uncName, uncertaintyDict in unc.iteritems():
            for binName, u in uncertaintyDict.iteritems():
                systDict[uncName][4][binName]["signal"] = u
        self.systs = sorted(systDict.values())
    
    def newSignalDiphoton(self, exp, unc):
        for bName, newRate in exp.iteritems():
            self.exp[bName]["t5gg"] = newRate
        systDict = dict([(l[0],l) for l in self.systs])
        #~ print systDict
        for uncName, uncertaintyDict in unc.iteritems():
            for binName, u in uncertaintyDict.iteritems():
                systDict[uncName][4][binName]["t5gg"] = u
        self.systs = sorted(systDict.values())
    
    def newSignalLepton(self, exp, unc):
        for bName, newRate in exp.iteritems():
            self.exp[bName]["SUSY"] = newRate
        systDict = dict([(l[0],l) for l in self.systs])
        #~ print systDict
        for uncName, uncertaintyDict in unc.iteritems():
            for binName, u in uncertaintyDict.iteritems():
                systDict[uncName][4][binName]["SUSY"] = u
        self.systs = sorted(systDict.values())

    def limit(self):
        self.write("/tmp/tmpDataCard.txt")
        return infosFromDatacard("/tmp/tmpDataCard.txt")

    def limitFast(self):
        infos = { "obs":0, "exp":0, "exp1up":0, "exp1dn":0, "exp2up":0, "exp2dn":0 }

        for bin in self.bins:
            obs = self.obs[bin]
            bkg = sum([b for a,b in self.exp[bin].iteritems() if a is not "sig"])
            signal = self.exp[bin]["sig"]
            if not signal: continue
            err = ROOT.TMath.Sqrt(bkg)
            r = signal/abs(err - abs(obs - bkg))
            infos["obs"] = max(infos["obs"], r)
            r_exp = signal/err
            if r_exp > infos["exp"]:
                infos["exp"] = r_exp
                infos["exp1up"] = (signal+err)/err
                infos["exp2up"] = (signal+2*err)/err
                infos["exp1dn"] = (signal-err)/err
                infos["exp2dn"] = (signal-2*err)/err
        return infos

    def setExpection(self):
        for binName, expDict in self.exp.iteritems():
            totRate = 0
            for process, rate in expDict.iteritems():
                totRate += 0 if process == "signal" else rate
            self.obs[binName] = round(totRate)
    
    def getStatUncertainty(self,binName,processName):
        systDict = dict([(l[0],l) for l in self.systs])
        temp = 0
        uncName = "statB"+binName.split("n")[1]+processName
        if uncName in systDict.keys():
            temp = systDict[uncName][4][binName][processName]
            temp = (temp-1)*self.exp[binName][processName]
        return temp
    
    def getUncertainty(self,uncName,binName,processName):
        systDict = dict([(l[0],l) for l in self.systs])
        temp = 0
        if uncName in systDict.keys():
            temp = systDict[uncName][4][binName][processName]
            temp = (temp-1)*self.exp[binName][processName]
        
        return temp
    def getUncertaintyGamma(self,uncName,binName,processName):
        systDict = dict([(l[0],l) for l in self.systs])
        temp = 0
        if uncName in systDict.keys():
            temp = systDict[uncName][4][binName][processName]
        return temp
        
    def renameUncertainty(self,uncNameOld,uncNameNew):
        found=False
        for line in self.systs:
            if line[0]==uncNameOld:
                self.systs[self.systs.index(line)][0]=uncNameNew
                found=True
                break
        return found

