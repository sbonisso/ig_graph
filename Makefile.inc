# put common definitions in here
CC	= g++
PRJCFLAGS	= -g
LD	= g++
LDFLAGS	= 
AR	= ar
ARFLAGS	=
RANLIB	= ranlib
RM	= rm -f
ECHO	= echo
CFLAGS	= -c -Wall -std=c++11 -I. -O3
DBGFLAGS=-c -Wall -std=c++11 -I. -g -DDEBUG
SHELL	= /bin/sh
.SILENT :