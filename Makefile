# Makefile for mcmc
# v1 DP 30 Nov 2018

#For the GNU C compiler
CC=gcc 
CFLAGS=-w

#Libraries
LIBSGEN=-lm

#Executables
EXEC=mcmc

#all rule
all: $(EXEC)

readdata.o: readdata.c
	$(CC) $(CFLAGS) -c readdata.c $(LIBSGEN)

twister.o: twister.c
	$(CC) $(CFLAGS) -c twister.c $(LIBSGEN)

chain.o: chain.c twister.o
	$(CC) $(CFLAGS) -c chain.c twister.o $(LIBSGEN)

mcmc: mcmc.c chain.o readdata.o twister.o
	$(CC) $(CFLAGS) mcmc.c chain.o readdata.o twister.o -o mcmc  $(LIBSGEN)

clean:
	rm -f *.o *.trace *~

.PHONY : all clean
