#!/usr/bin/env python2
import suppressor
with suppressor.suppress_stdout_stderr():
    import ROOT
import ROOT
import argparse
import os
import subprocess

def run(infile=""):
    extName = infile.replace("_nTuple", "_ext_nTuple")
    if os.path.isfile(extName):
        print "Add file", extName
        subprocess.call(["./CombinationHistogramProducer", infile, extName])
    else:
        subprocess.call(["./CombinationHistogramProducer", infile])

def runExt(infile="", selector="HistogramProducer.cc"):
    # wrapper
    run(infile, selector, True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file', default="", nargs="?")
    parser.add_argument('--ext', action='store_true')
    parser.add_argument('--signal', action='store_true')

    args = parser.parse_args()
    signalScans = ["SMS-T5Wg_nTuple.root", "SMS-T6Wg_nTuple.root", "SMS-T5Wg_mGo2150To2500_nTuple.root", "SMS-T6Wg_mSq1850To2150_nTuple.root", "SMS-TChiWG_nTuple.root", "SMS-TChiNG_nTuple.root"]
    if os.path.basename(args.file) in signalScans or args.signal:
        run(args.file)
    else:
        run(args.file)
