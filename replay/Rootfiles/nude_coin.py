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

#trigflag = 5 # (coin)
#trigflag = 4 # (RHRS)
#trigflag = 1 # (LHRS)
nworkers=5
#runfile = "h2_1.dat"
#runfile = "h2_2.dat"
#runfile = "h2_3.dat"
#runfile = "h2_4.dat"
runfile = "temp.dat"
#runfile = "h2_5.dat"
#runfile = "h2_6.dat"
#runfile = "multiF.dat"
#runfile = "h22_1.dat"
#runfile = "h22_2.dat"
#runfile = "h22_3.dat"
#runfile = "h2_2_temp.dat"
#runfile = "ss.dat"
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
        #com = "root -l -q \"nude3.cc(" + data[0]+ "," +data[1] + "," + str(trigflag)
        #com2 = com + ")\";"
        #call(com,shell=True)
        #comlist.append(com2)
        comlist.append(com)
    with ProcessPoolExecutor(max_workers=nworkers) as executor:
        executor.map(nude_start,comlist)


stime = time.time()
main()
print("\n Jobs were done in %.0f sec \n" % float(time.time()-stime))

