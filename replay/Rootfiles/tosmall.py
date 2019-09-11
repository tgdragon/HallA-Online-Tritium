#!/apps/python/3.4.3/bin/python3
#!/usr/bin/python

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
#runfile = "multiF.dat"
runfile = "h2.dat"
#runfile = "h22.dat"
#runfile = "He3.dat"
#runfile = "T2.dat"
#runfile = "ss.dat"
#runfile = "Al_h22.dat"
#runfile = "Al_h2.dat"
#runfile = "Al_He3.dat"
#runfile = "Al_T2.dat"
#runfile = "ssR.dat"

thisfile = "nude_coin.py"

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
        com = "./nude_coin " + data[0]+ " " +data[1]
        #com = "./maruhadaka " + data[0]+ " " +data[1] + " " +data[2] #stricter cut
        #com = "./maruhadakaR " + data[0]+ " " +data[1] + " " +data[2] #stricter cut
        #com = "root -l -q \"nude3.cc(" + data[0]+ "," +data[1] + "," + str(trigflag)
        #com2 = com + ")\";"
        #call(com,shell=True)
        #comlist.append(com2)
        comlist.append(com)
    with ProcessPoolExecutor(max_workers=nworkers) as executor:
        executor.map(nude_start,comlist)

    #com3 = "./merge_rootfiles" + runfile
    #call(com3, shell=True)


stime = time.time()
main()
print("\n Jobs were done in %.0f sec \n" % float(time.time()-stime))

