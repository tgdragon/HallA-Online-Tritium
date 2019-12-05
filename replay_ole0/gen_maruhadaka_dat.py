#!/usr/bin/python

# Toshiyuki Gogami
# December 5, 2019

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
        for line2 in lines2:
            data2 = line2.rsplit(' ',1)
            length_data2 = len(data2)
            #print(data2[length_data2-1][-15:])
            f.write(data2[length_data2-1])
            print(data2[length_data2-1])
    com2 = "rm -f temptemp"
    call (com2, shell=True)
        

stime = time.time()
main()
print("\n Jobs were done in %.0f sec \n" % float(time.time()-stime))
