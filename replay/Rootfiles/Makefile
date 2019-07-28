#Makefile


#ROOTLIBS = -L/u/apps/root/6.10.02/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic
ROOTLIBS    = `root-config --libs`
ROOTINCLUDE = `root-config --cflags`
GLIBS = $(ROOTLIBS) -lMinuit
CC = g++
#CCFLAGS = -pthread -std=c++11 -m64 -I/u/apps/root/6.10.02/root/include
CCFLAGS = $(ROOTINCLUDE)
#LD = g++
#LD = /home/gogami/gcc/gcc-6.3.0-install/bin/g++
LDFLAGS = -g
TARGET1 = nude_coin
OBJS1 = nude_coin.o

all: ${TARGET1} 

${TARGET1}: ${OBJS1}
	${CC} ${CCFLAGS} -o $@ ${OBJS1} ${GLIBS}
${OBJS1}: nude_coin.cc
	${CC} ${CCFLAGS} -c  nude_coin.cc ${GLIBS}



#.cc.o:
#	${CC} ${CCFLAGS} -c $<

clean:
	rm -rf *.o core ${TARGET1}

