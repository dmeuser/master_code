import argparse
import os
import runCombine
import multiprocessing
import subprocess
import sys
import time
import tqdm
from functools import partial

parser = argparse.ArgumentParser()
parser.add_argument('scan', nargs='?', help="choose a signal scan")
parser.add_argument('selection', nargs='?', help="choose as selection like leptonVeto, htgVeto etc.")
parser.add_argument('--condor', action="store_true", help="Limits are calculated on condor")
parser.add_argument('--jobsize', type=int, default=10, help="Number of Datacards in one condor job (default=10)")
args = parser.parse_args()

dir="input/"+args.scan+"_"+args.selection

toProcess=[]

for point in os.listdir(dir):
    toProcess.append(dir+"/"+point)


if args.condor:
    subprocess.call(["rm","condor/submits/*"])
    subprocess.call(["rm","condor/logs/*"])
    toProcessSeq=[toProcess[i:i+args.jobsize] for i  in range(0, len(toProcess), args.jobsize)]
    for part in toProcessSeq:
        toProcessString=""
        for x in part:
            toProcessString+=x+","
        with open("condor/submits/submitCondor_"+part[0].split("/")[2],"w") as f:
            f.write("""
Universe   = vanilla
Executable = runCombine.sh
Arguments  = {0} {1} {2}
Log        = condor/logs/{3}.log
Output     = condor/logs/{3}.out
Error      = condor/logs/{3}.error
Queue
""".format(toProcessString, args.scan, args.selection, args.selection+"_"+(part[0].split("/")[2])))
        subprocess.call(["condor_submit", "condor/submits/submitCondor_"+part[0].split("/")[2]])

else: # local processing
    p = multiprocessing.Pool()
    
    for _ in tqdm.tqdm(p.imap_unordered(partial(runCombine.runCombine,selection=args.selection,scan=args.scan),toProcess), total=len(toProcess)):
        pass
