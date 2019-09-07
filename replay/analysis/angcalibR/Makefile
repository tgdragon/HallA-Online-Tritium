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

TARGET1 = angcalibR
OBJS1 = angcalibR.o
#TARGET2 = check
#OBJS2 = check.o

#all: ${TARGET1}  ${TARGET2} 
all: ${TARGET1}

${TARGET1}: ${OBJS1}
	${CC} ${CCFLAGS} -o $@ ${OBJS1} ${GLIBS}
${OBJS1}: angcalibR.cc
	${CC} ${CCFLAGS} -c  angcalibR.cc ${GLIBS}

#${TARGET2}: ${OBJS2}
#	${CC} ${CCFLAGS} -o $@ ${OBJS2} ${GLIBS}
#${OBJS2}: check.cc
#	${CC} ${CCFLAGS} -c  check.cc ${GLIBS}



#.cc.o:
#	${CC} ${CCFLAGS} -c $<

clean:
	rm -rf *.o core ${TARGET1} ${TARGET2}

