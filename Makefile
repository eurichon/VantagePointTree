# define the shell to bash
SHELL := /bin/bash

# define the C/C++ compiler to use,default here is clang
CC = gcc-7


test_sequential:	
	#tar -xvzf code.tar.gz
	cd vptree; make lib; cd ..
	cd vptree; cp lib/*.a inc/vptree.h ../; cd ..
	$(CC) tester.c vptree_sequential.a -o $@ -lm
	./test_sequential

test_pthreads:
	#tar -xvzf code.tar.gz
	cd vptree; make lib; cd ..
	cd vptree; cp lib/*.a inc/vptree.h ../; cd ..
	$(CC) tester.c vptree_pthreads.a -o $@ -lm
	./test_pthreads

test_openmp:
	tar -xvzf code.tar.gz
	cd vptree; make lib; cd ..
	cd vptree; cp lib/*.a inc/vptree.h ../; cd ..
	$(CC) tester.c vptree_openmp.a -o $@ -lm
	./test_openmp

test_cilk:
	tar -xvzf code.tar.gz
	cd vptree; make lib; cd ..
	cd vptree; cp lib/*.a inc/vptree.h ../; cd ..
	$(CC) tester.c vptree_cilk.a -o $@ -lm
	./test_cilk
	
clean: 
	rm vptree_seq.a
	rm vptree_sequential.a
	rm test_sequential
	