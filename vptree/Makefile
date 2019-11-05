SHELL := /bin/bash

# define the C/C++ compiler to use,default here is clang
CC = gcc

# define compile-time flags
CFLAGS = -std=gnu99  -O3  -lpthread  -fopenmp -Wall


INCLUDES =	-I inc

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib specify
#   their path using -Lpath, something like:


all:	src/main.c src/vptree_sequential.c src/vptree_pthreads.c src/vptree_openmp.c src/vptree_cilk.c
	$(CC) $(CFLAGS) $(INCLUDES) -c src/vptree_sequential.c -o  src/vptree_sequential.o -lm
	ar rcs lib/vptree_sequential.a src/vptree_sequential.o
	rm src/vptree_sequential.o
	$(CC) $(CFLAGS) $(INCLUDES) -o main_sequential src/main.c lib/vptree_sequential.a -lm
	$(CC) $(CFLAGS) $(INCLUDES) -c src/vptree_pthreads.c -o  src/vptree_pthreads.o -lm
	ar rcs lib/vptree_pthreads.a src/vptree_pthreads.o
	rm src/vptree_pthreads.o
	$(CC) $(CFLAGS) $(INCLUDES) -o main_pthreads src/main.c lib/vptree_pthreads.a -lm
	$(CC) $(CFLAGS) $(INCLUDES) -c src/vptree_openmp.c -o  src/vptree_openmp.o -lm
	ar rcs lib/vptree_openmp.a src/vptree_openmp.o
	rm src/vptree_openmp.o
	$(CC) $(CFLAGS) $(INCLUDES) -o main_openmp src/main.c lib/vptree_openmp.a -lm
	
	#$(CC) $(CFLAGS) $(INCLUDES) -c src/vptree_cilk.c -o  src/vptree_cilk.o -lm
	#ar rcs lib/vptree_cilk.a src/vptree_cilk.o
	#rm src/vptree_cilk.o
	#$(CC) $(CFLAGS) $(INCLUDES) -o main_cilk src/main.c lib/vptree_cilk.a -lm

lib:	src/vptree_sequential.c src/vptree_pthreads.c src/vptree_openmp.c src/vptree_cilk.c
	$(CC) $(CFLAGS) $(INCLUDES) -c src/vptree_sequential.c -o  src/vptree_sequential.o -lm
	ar rcs lib/vptree_sequential.a src/vptree_sequential.o
	rm src/vptree_sequential.o
	$(CC) $(CFLAGS) $(INCLUDES) -c src/vptree_pthreads.c -o  src/vptree_pthreads.o -lm
	ar rcs lib/vptree_pthreads.a src/vptree_pthreads.o
	rm src/vptree_pthreads.o
	$(CC) $(CFLAGS) $(INCLUDES) -c src/vptree_openmp.c -o  src/vptree_openmp.o -lm
	ar rcs lib/vptree_openmp.a src/vptree_openmp.o
	rm src/vptree_openmp.o
	$(CC) $(CFLAGS) $(INCLUDES) -c src/vptree_cilk.c -o  src/vptree_cilk.o -lm
	ar rcs lib/vptree_cilk.a src/vptree_cilk.o
	rm src/vptree_cilk.o
