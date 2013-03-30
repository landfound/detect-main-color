all: maincolor 
debug: maincolord

OBJS=maincolor.o
OBJSd=maincolord.o
CC= g++
CFLAGS=-Wall -O3 -g
LIBS=-lcv -lcxcore -lhighgui -lcvaux -L/opt/opencv/lib
INCLUDE=-I/opt/opencv/include/opencv

maincolor: $(OBJS)
	$(CC) -o maincolor $(OBJS) $(CFLAGS) $(LIBS) $(INCLUDE) 

maincolor.o: maincolor.cpp
	$(CC) $(CFALGS) $(INCLUDE) -c maincolor.cpp -o maincolor.o

maincolord: $(OBJSd)
	$(CC) -o maincolord $(OBJSd) $(CFLAGS) $(LIBS) $(INCLUDE) 

maincolord.o: maincolor.cpp
	$(CC) $(CFALGS) $(INCLUDE) -D_DEBUG -c maincolor.cpp -o maincolord.o
