#
# Makefile
# uqcskenn, 2013-11-22 00:01
#
CC = gcc
CFLAGS = -g -O2 -Wall -pthread
HTS = /usr/local/lib#$(HOME)/local/htslib
LIBS = -lz -lpthread -lm -lcfu
EXECUTABLE = bamParser

SOURCES = example.c bamParser.c pairedLink.c

$(EXECUTABLE): $(SOURCES)
	$(CC) $(CFLAGS) -L$(HTS) -I$(HTS) -o $@ $^ $(HTS)/libhts.a $(LIBS)

all: $(EXECUTABLE)


# vim:ft=make
#

