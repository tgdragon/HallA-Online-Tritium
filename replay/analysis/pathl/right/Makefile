#Makefile

ROOTLIBS = -L/u/apps/root/6.10.02/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic 
#ROOTLIBS = -L/home/gogami/root/root-5.34.36-install/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic
#ROOTINCLUDE = `root-config --cflags`
#ROOTINCLUDE = /home/gogami/root/root-5.34.36-install/include/root
#GLIBS = $(ROOTLIBS) $(ROOTGLIBS) -lm -lCore -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lGui -pthread -lm -ldl -rdynamic -lMinuit
GLIBS = $(ROOTLIBS) $(ROOTGLIBS) -lMinuit
CC = g++
#CC = /home/gogami/gcc/gcc-6.3.0-install/bin/g++
#CCFLAGS = -O  -D__osf__ -D__alpha -I$(ROOTINCLUDE)
#CCFLAGS = -I$(ROOTINCLUDE) -I./
CCFLAGS = -pthread -std=c++11 -m64 -I/u/apps/root/6.10.02/root/include
#LD = g++
#LD = /home/gogami/gcc/gcc-6.3.0-install/bin/g++
LDFLAGS = -g
#TARGET1 = getMatrix 
#OBJS1 =  getMatrix.o calcMatrix.o readvalue.o 
TARGET1 = pathl_opt
OBJS1 =  pathl_opt.o calcRecon.o

all: ${TARGET1} 

${TARGET1}: ${OBJS1}
	${CC} ${CCFLAGS} -o $@ ${OBJS1} ${GLIBS}

pathl_opt.o: pathl_opt.cc pathl_opt.hh
	${CC} ${CCFLAGS} -c  pathl_opt.cc ${GLIBS}

#getMatrix.o: getMatrix.C getMatrix.h 
#	${CC} ${CCFLAGS} -c  getMatrix.C ${GLIBS}
#
calcRecon.o: calcRecon.C pathl_opt.hh
	${CC} ${CCFLAGS} -c  calcRecon.C ${GLIBS}

#readvalue.o: readvalue.C 
#	${CC} ${CCFLAGS} -c  readvalue.C ${GLIBS}
#

#.cc.o:
#	${CC} ${CCFLAGS} -c $<

clean:
	rm -rf *.o core ${TARGET1} *~

