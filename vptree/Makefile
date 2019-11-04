SHELL := /bin/bash



CC       = gcc-7
FLAGS    = -O3 -pthread -Wall
CFLAGS = -std=gnu99  -O3  -lpthread  -fopenmp -Wall


INCLUDES =	-I inc

all:	src/main.c src/vptree_sequential.c src/vptree_pthreads.c
	$(CC) $(CFLAGS) $(INCLUDES) -c src/vptree_sequential.c -o  src/vptree_sequential.o -lm
	ar rcs lib/vptree_sequential.a src/vptree_sequential.o
	rm src/vptree_sequential.o
	$(CC) $(CFLAGS) $(INCLUDES) -o main_sequential src/main.c lib/vptree_sequential.a -lm
	
	$(CC) $(CFLAGS) $(INCLUDES) -c src/vptree_pthreads.c -o  src/vptree_pthreads.o -lm
	ar rcs lib/vptree_pthreads.a src/vptree_pthreads.o
	rm src/vptree_pthreads.o
	$(CC) $(CFLAGS) $(INCLUDES) -o main_pthreads src/main.c lib/vptree_pthreads.a -lm
	
	
lib:	src/vptree_sequential.c src/vptree_pthreads.c
	$(CC) $(CFLAGS) $(INCLUDES) -c src/vptree_sequential.c -o  src/vptree_sequential.o -lm
	ar rcs lib/vptree_sequential.a src/vptree_sequential.o
	rm src/vptree_sequential.o
	

	$(CC) $(CFLAGS) $(INCLUDES) -c src/vptree_pthreads.c -o  src/vptree_pthreads.o -lm
	ar rcs lib/vptree_pthreads.a src/vptree_pthreads.o
	rm src/vptree_pthreads.o

