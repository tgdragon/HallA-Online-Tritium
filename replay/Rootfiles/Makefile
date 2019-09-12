#Makefile


#ROOTLIBS = -L/u/apps/root/6.10.02/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic
ROOTLIBS    = `root-config --libs`
ROOTINCLUDE = `root-config --cflags`
GLIBS = $(ROOTLIBS) -lMinuit
CC = g++
#CC = icpc
#CCFLAGS = -pthread -std=c++11 -m64 -I/u/apps/root/6.10.02/root/include
CCFLAGS = $(ROOTINCLUDE)
#LD = g++
#LD = /home/gogami/gcc/gcc-6.3.0-install/bin/g++
LDFLAGS = -g

#TARGET1 = nude_coin
#OBJS1 = nude_coin.o

TARGET2 = merge_rootfiles
OBJS2 = merge_rootfiles.o

TARGET3 = maruhadaka
OBJS3 = maruhadaka.o

TARGET4 = maruhadakaR
OBJS4 = maruhadakaR.o

TARGET5 = maruhadakaL
OBJS5 = maruhadakaL.o



#all: ${TARGET1} ${TARGET2} ${TARGET3} ${TARGET4} ${TARGET5}
all: ${TARGET2} ${TARGET3} ${TARGET4} ${TARGET5}

#${TARGET1}: ${OBJS1}
#	${CC} ${CCFLAGS} -o $@ ${OBJS1} ${GLIBS}
#${OBJS1}: nude_coin.cc
#	${CC} ${CCFLAGS} -c  nude_coin.cc ${GLIBS}

${TARGET2}: ${OBJS2}
	${CC} ${CCFLAGS} -o $@ ${OBJS2} ${GLIBS}
${OBJS2}: merge_rootfiles.cc
	${CC} ${CCFLAGS} -c  merge_rootfiles.cc ${GLIBS}

${TARGET3}: ${OBJS3}
	${CC} ${CCFLAGS} -o $@ ${OBJS3} ${GLIBS}
${OBJS3}: maruhadaka.cc
	${CC} ${CCFLAGS} -c  maruhadaka.cc ${GLIBS}

${TARGET4}: ${OBJS4}
	${CC} ${CCFLAGS} -o $@ ${OBJS4} ${GLIBS}
${OBJS4}: maruhadakaR.cc
	${CC} ${CCFLAGS} -c  maruhadakaR.cc ${GLIBS}

${TARGET5}: ${OBJS5}
	${CC} ${CCFLAGS} -o $@ ${OBJS5} ${GLIBS}
${OBJS5}: maruhadakaL.cc
	${CC} ${CCFLAGS} -c  maruhadakaL.cc ${GLIBS}




#.cc.o:
#	${CC} ${CCFLAGS} -c $<

clean:
	rm -rf *.o core ${TARGET1} ${TARGET2}

