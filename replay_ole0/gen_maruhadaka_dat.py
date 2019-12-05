#!/usr/bin/python

# This file will create dat file which is used for Rootfiles/maruhadaka.cc
# 1. Replay data files using gogogo.py in which you specify *_replay.dat
# 2. Run this script (*_replay.dat --> *.dat)
# 3. Copy *.dat to Rootfiles/ to be used for Rootfiles/maruhadaka.cc 

# --- Toshiyuki Gogami ----  #
# --- December 5, 2019 ----  #

import sys
import time, os.path
from subprocess import call
#import concurrent.futures
from logging import StreamHandler, Formatter, INFO, getLogger
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures.process import ProcessPoolExecutor 
import numpy as np

target = "h2"
#target = "h22"
#target = "T2"
#target = "He3"

runfile = target+"_replay.dat" 
outfile = target+".dat"

if target=="h2":
    hypid = 1
elif target=="h22":
    hypid = 2
elif target=="T2": 
    hypid = 3
elif target=="He3": 
    hypid = 4
else:
    hypid = 27

def main():
    inputfile = open(runfile,"r")
    lines = inputfile.readlines()
    f = open(outfile,"w")
    for line in lines:
        data = line.split()
        com = "ls -l Rootfiles/nnL/*" + data[0] + "*.root > temptemp"
        call(com, shell=True)
        
        inputfile2 = open("temptemp","r")
        lines2 = inputfile2.readlines()
        fnum = 0
        for line2 in lines2:
            data2 = line2.rsplit(' ',1)
            length_data2 = len(data2)
            #print(data2[length_data2-1][-15:])
            #f.write(data2[length_data2-1])
            #print(data2[length_data2-1])
            fnum = fnum +1
        print(data[0],fnum)
        
        fnum_this = 0
        for i in range(fnum):
            content = data[0]+" "+str(fnum_this)+" "+str(hypid)+"\n"
            f.write(content)
            fnum_this = fnum_this + 1

    com2 = "rm -f temptemp"
    call (com2, shell=True)
        

stime = time.time()
main()
print("\n Jobs were done in %.0f sec \n" % float(time.time()-stime))
