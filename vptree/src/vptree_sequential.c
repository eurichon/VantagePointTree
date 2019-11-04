#include "vptree.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

void createRecursivly(vptree *T, int index, double **arrayset, double *distances, int start, int end);
void createTree(vptree *T, double **arraySet, double *distances);
double partition(double **arraySet, double *distances, int low, int high);
double kthSmallest(double **arraySet, double *distances, int left, int right, int k);
void getDistanceSquared(double *p1, double *p2, double *result,int dim);
void swap2(double *array, int i, int j);
void swap(double **n1, double **n2);


/* Trasforms the data in efficient form and proceeds to create the Tree*/
vptree *buildvp(double *X, int n, int d){
	unsigned i,j;
	double *ptr;
	
	
	int length = n * sizeof(double *) + n * d * sizeof(double);
	double **dataset = (double **)malloc(length);
	if(dataset == NULL){
		printf("Not enought memory for creation of dataset table.Exiting....\n");
		exit(-1);
	}
	
	//printf("Dataset is: ");
	ptr = (double *)(dataset + n);
	for (i = 0; i < n; i++) {
		dataset[i] = (double *)(ptr + d * i);
		
		for (j = 0;j < d; j++) {
			dataset[i][j] = X[i * d + j];
			//printf("%lf, ",dataset[i][j]);
		}
	}
	//printf("\n\n");
	
	//printf("Auxilary dataset initialized\n\n");
	
	struct timespec start, finish;
	double elapsed;
	
	clock_gettime(CLOCK_MONOTONIC, &start);
	unsigned max_nodes = (int)pow(2, ((int)log2(n) + 1)) - 1;
	vptree *T = (vptree *)malloc(sizeof(vptree));
	T->tree = (double **)calloc(max_nodes * (d+3),sizeof(double));
	if(T->tree == NULL){
		printf("Not enought memory for creation of the tree.Exiting....\n");
		exit(-1);
	}
	ptr = (double *)(T->tree + max_nodes);
	
	for(i = 0; i < max_nodes; ++i){
		T->tree[i] = (ptr + (d + 2) * i); 
		T->tree[i][d + 1] = -1.0; 
	}
	
	
	
	T->index = 0;
	T->d = d;
	T->n = n;
	
	
	//create tree
	double *distances = (double *)malloc(max_nodes * sizeof(double)); //auxilary table to hold the distances
	if(distances == NULL){
		printf("Not enought memory for the creation of distances.Exiting....\n");
		exit(-1);
	}
	
	
	createTree(T, dataset ,distances);
	
	
	
	/*
	printf("\nTree structure is\n");
	for(int i = 0;i < max_nodes; ++i){
		printf("Vantage point: %lf with median: %lf index %i\n",T->tree[i][0],T->tree[i][1],(int)T->tree[i][2]);
	}*/
	
	
	clock_gettime(CLOCK_MONOTONIC, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("%i %lf\n", n,elapsed *4 / 5);
	
	
	free(distances);
	free(dataset);	
	return T;	
}

vptree * getInner(vptree * T){
	unsigned max_size = (int)pow(2, ((int)log2(T->n) + 1)) - 1;
	if(max_size <= (2 * T->index + 1)){
		return NULL;
	}else{	
		T->index = (2 * T->index + 1);
		if((int)T->tree[T->index][T->d+1] ==  -1){
			return NULL;
		}else{
			return T;
		}
	}
}

vptree * getOuter(vptree * T){
	unsigned max_size = (int)pow(2, ((int)log2(T->n) + 1)) - 1;
	if(max_size <= (2 * T->index + 2)){
		return NULL;
	}else{	
		T->index = (2 * T->index + 2);
		if((int)T->tree[T->index][T->d+1] ==  -1){
			return NULL;
		}else{
			return T;
		}
	}
}

double getMD(vptree * T){
	return (T->tree[T->index][T->d]);
}

double *getVP(vptree * T){
	return (T->tree[T->index]);
}

int getIDX(vptree * T){
	return (T->index);
}

void createTree(vptree *T, double **arraySet, double *distances) {
	unsigned int index_id = 0;
	createRecursivly(T, index_id, arraySet, distances, 0, T->n);
}

void createRecursivly(vptree *T, int index, double **arrayset, double *distances, int start, int end) {
	unsigned int size = end - start;
	//printf("Node: %i ,with size: %i from index [%i, %i)\n",index,size,start,end);
	if(size <= 0){
		return;
	}else{
		T->tree[index][T->d + 1] = index;  //store index id
		memcpy(T->tree[index], arrayset[end - 1], T->d * sizeof(double)); //copy v_point
		int mean = (end - start - 1) / 2;
		
		//printf("Distances before quickselect: ");
		for(unsigned i = start; i < (end - 1); ++i){
			getDistanceSquared(arrayset[i], T->tree[index], &distances[i], T->d);
			//printf("%lf, ",distances[i]);
		}//printf("\n");
		
		int k = (size - 1)/2;
		if(k == 0)
			k = 1;
		T->tree[index][T->d] = kthSmallest(arrayset, distances, start, end - 2, k);
		//printf("Median: %lf\n",T->tree[index][T->d]);
		
		//printf("Vantage point: %lf\n",T->tree[index][0]);
		
		/*
		printf("Distances after quickselect: ");
		for (unsigned i = start; i < (end - 1); ++i) {
			printf("%lf, ",distances[i]);
		}printf("\n");
		
		printf("Dataset change to: ");
		for (unsigned i = start; i < (end - 1); ++i) {
			printf("%lf, ",arrayset[i][0]);
		}printf("\n");
		
		
		printf("=========================================\n");*/
		createRecursivly(T, (2 * index + 1), arrayset, distances, start, start + mean);
		createRecursivly(T, (2 * index + 2), arrayset, distances, start + mean , (end - 1));
	}
}

void swap(double **n1, double **n2) {
	double *temp = *n1;
	*n1 = *n2;
	*n2 = temp;
}

void swap2(double *array, int i, int j) {
	double temp = array[i];
	array[i] = array[j];
	array[j] = temp;
}

double partition(double **arraySet, double *distances, int low, int high) 
{ 
    int pivot = distances[high]; 
    int i = (low - 1); 
    for (int j = low; j <= high - 1; j++) { 
        if (distances[j] <= pivot) { 
            i++; 
            swap2(distances, i, j); 
			//swap(&arraySet[i], &arraySet[j]); 
        } 
    } 
    swap2(distances, (i+1), high); 
	//swap(&arraySet[i + 1], &arraySet[high]); 
    return (i + 1); 
} 
  
double kthSmallest(double **arraySet, double *distances, int left, int right, int k) 
{ 
	int temp = left;
    while (left <= right) { 
        int pivotIndex = partition(arraySet, distances, left, right); 
		//printf("inside k largest %i %i\n",pivotIndex,k-1);
        if ((pivotIndex - temp)  == k-1) 
            return (distances[pivotIndex]); 
        else if ((pivotIndex - temp) > k-1) 
            right = pivotIndex - 1; 
        else
            left = pivotIndex + 1; 
    } 
    return 0.0; 
} 

void getDistanceSquared(double *p1, double *p2, double *result,int dim)
{
	(*result) = 0;
	for (unsigned int i = 0;i < dim; ++i) {
		(*result) = (*result) + pow((p1[i] - p2[i]), 2);
	}
	(*result) = sqrt((*result));
}


