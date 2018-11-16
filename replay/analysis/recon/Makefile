#Makefile

#ROOTLIBS = `root-config --libs`
#ROOTINCLUDE = `root-config --cflags`
#ROOTLIBS = -L/home/gogami/root/root-5.34.36-install/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic
ROOTLIBS = -L/u/apps/root/6.10.02/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic
#ROOTINCLUDE = /home/gogami/root/root-5.34.36-install/include/root
#ROOTLIBS = `-L/home/gogami/root/root-5.34.36-install/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic`
#ROOTINCLUDE = `-m64 -I/home/gogami/root/root-5.34.36-install/include/root`
#GLIBS = $(ROOTLIBS) $(ROOTGLIBS) -lm -lCore -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lGui -pthread -lm -ldl -rdynamic -lMinuit
#GLIBS = $(ROOTLIBS) $(ROOTGLIBS) -lm -lCore -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lGui -pthread -lm -ldl -rdynamic -lMinuit
#GLIBS = $(ROOTLIBS) $(ROOTGLIBS) -lMinuit
GLIBS = $(ROOTLIBS) -lMinuit
CC = g++
#CC = /home/gogami/gcc/gcc-6.3.0-install/bin/g++
#CCFLAGS = -O  -D__osf__ -D__alpha -I$(ROOTINCLUDE)
#CCFLAGS = -I$(ROOTINCLUDE) -I./
#CCFLAGS = -pthread -m64 -I$(ROOTINCLUDE)
CCFLAGS =-pthread -std=c++11 -m64 -I/u/apps/root/6.10.02/root/include
#LD = g++
#LDFLAGS = -g
TARGET1 = recon 
OBJS1 =  recon.o calcRecon.o readvalue.o 

all: ${TARGET1} 

${TARGET1}: ${OBJS1}
	${CC} ${CCFLAGS} -o $@ ${OBJS1} ${GLIBS}

readvalue.o: readvalue.C 
	${CC} ${CCFLAGS} -c  readvalue.C ${GLIBS}

recon.o: recon.C recon.h
	${CC} ${CCFLAGS} -c  recon.C ${GLIBS}

calcRecon.o: calcRecon.C recon.h
	${CC} ${CCFLAGS} -c  calcRecon.C ${GLIBS}

#.cc.o:
#	${CC} ${CCFLAGS} -c $<

clean:
	rm -rf *.o core ${TARGET1}

