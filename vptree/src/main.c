#include <stdlib.h>
#include <stdio.h>
#include "vptree.h"
#include <time.h>

#define SET_SIZE 10000
#define POINT_DIM 2
#define MIN -1000.0
#define MAX 1000.0

double getRandom(float min, float max);


int main(int argc, char *argv[]) {
	int n, d;
	
	if(argc != 3){
		printf("Wrong Number of arguments %i!\nSettign default...\n",argc);
		n = 100000;
		d = 2;
	}else{
		n = atoi(argv[1]);
		d = atoi(argv[2]);
	}
	
	srand((int)time(NULL));
	
	unsigned i;
	

	
		
	//================array initialization================
	double *data = (double *)malloc(n*d*sizeof(double));
	for(i = 0;i < n*d;++i){
		data[i] = getRandom(MIN,MAX);
	}
	if(data == NULL){
		printf("Not enought memory for creation of data.Exiting....\n");
		exit(-1);
	}
	//printf("The dataset has been initialized with %i points of %i dimensions\n", n, d);

	struct timespec start, finish;
	double elapsed;
	
	clock_gettime(CLOCK_MONOTONIC, &start);
		
	vptree *T = buildvp(data,n,d);
	
	
	clock_gettime(CLOCK_MONOTONIC, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("Time is %lf\n", elapsed);
	

	free(data);
	
	//getchar();
	return 0;
}



double getRandom(float min, float max) {
	float random = ((float)rand()) / (float)RAND_MAX;
	float diff = max - min;
	float r = random * diff;
	return (min + r);
}


