#!/usr/bin/python
#!/apps/python/3.4.3/bin/python3


# Toshiyuki Gogami
# Nov 2, 2018

import sys
import time, os.path
from subprocess import call
#import concurrent.futures
from logging import StreamHandler, Formatter, INFO, getLogger
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures.process import ProcessPoolExecutor 
import numpy as np

nworkers=5
#nworkers=15
#runfile = "all_replay.dat"
#runfile = "h2_replay.dat"
#runfile = "h22_replay.dat"
#runfile = "He3_replay.dat"
#runfile = "T2_replay.dat"
runfile = "dummy_replay.dat"

thisfile = "gogogo.py"

def nude_start(command):
    time.sleep(1.0)
    call(command,shell=True)

def main():
    comlist = []
    #inputfile = open("h2_2.dat","r")
    inputfile = open(runfile,"r")
    lines = inputfile.readlines()
    for line in lines:
        data = line.split()
        com = "./fullReplay_gogogo " + data[0]
        comlist.append(com)
    with ProcessPoolExecutor(max_workers=nworkers) as executor:
        executor.map(nude_start,comlist)

    #com3 = "./merge_rootfiles" + runfile
    #call(com3, shell=True)


stime = time.time()
main()
print("\n Jobs were done in %.0f sec \n" % float(time.time()-stime))

