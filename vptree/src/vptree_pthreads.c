#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "vptree.h"

#define NUM_OF_THREAD 50

int curr_thread;
pthread_mutex_t mutexId;
pthread_attr_t attr;
pthread_t thread[NUM_OF_THREAD];


/* Usefull struct to pass the data from the recursive function to the thread*/
typedef struct{
	vptree *T;
	int index,
		start,
		end;
	double **arrayset,
			*distances;
}TreeData;


//====================================Vptree functions declaration===========================================
void createRecursivly(vptree *T, int index, double **arrayset, double *distances, int start, int end);
void createTree(vptree *T, double **arraySet, double *distances);
double partition(double **arraySet, double *distances, int low, int high);
double kthSmallest(double **arraySet, double *distances, int left, int right, int k);
void getDistanceSquared(double *p1, double *p2, double *result,int dim);
void swap2(double *array, int i, int j);
void swap(double **n1, double **n2);
void *recursiveJob(void *data);


//====================================Vptree functions implementation===========================================
/*Builds a tree from a dataset with n points of m dimentions*/
vptree *buildvp(double *X, int n, int d){
	unsigned i,j;
	double *ptr;
	
	/* We copy the data to another 2 dimentional array which is stored in one row.
	In the first n positions we set the pointers to the n points of d dimentions.
	However in this implementation a swap between 2 points always happens in O(1)
	cause we simply swap the pointers in the first n potision segment*/
	int length = n * sizeof(double *) + n * d * sizeof(double);
	double **dataset = (double **)malloc(length);
	if(dataset == NULL){
		printf("Not enought memory for creation of dataset table.Exiting....\n");
		exit(-1);
	}
	ptr = (double *)(dataset + n);
	for (i = 0; i < n; i++) {
		dataset[i] = (double *)(ptr + d * i);
		
		for (j = 0;j < d; j++) {
			dataset[i][j] = X[i * d + j];
		}
	}
	
	/* We allocate the necessary memory for the tree structure in an array representation
	For n element we will have log2(n) levels so the total Node will be 2^0 + 2^1 + ... + 2^log2(n) = 2^(log2(n) + 1) - 1
	We also keep another exta two for each point to store the median and its position in the tree*/
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
	
	/* This an auxilary table for the calculations of the median so as we dont need to
	allocate and deallocate memory constantly. Also it eases the need for offset in the partitioning*/
	double *distances = (double *)malloc(max_nodes * sizeof(double)); 
	if(distances == NULL){
		printf("Not enought memory for the creation of distances.Exiting....\n");
		exit(-1);
	}
	
	//==============================parallel segment===================================
	void *status;
	pthread_mutex_init(&mutexId, NULL);
	pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	
	curr_thread = 0;
	printf("Starting...\n");
	createTree(T, dataset ,distances);
	
	/* Release the threads with a join*/
	pthread_attr_destroy(&attr);
	for(i=0; i< curr_thread; i++) {
        int rc = pthread_join(thread[i], &status);
        if (rc) {
			printf("ERROR; return code from pthread_join() is %d\n", rc);
			exit(1);
        }
    }
	
	/* free the auxilary table */
	free(distances);
	free(dataset);	
	
	return T;	
}

/*return the inner subtree with all its contents: vpoint, median, index.
In an array representation of a tree the left(inside) child of the Node i is the (2 * i + 1) position.
We also check if the index gets out of boundary*/
vptree * getInner(vptree * T){
	unsigned max_size = (int)pow(2, ((int)log2(T->n) + 1)) - 1;
	if(max_size <= (2 * T->index + 1)){
		return NULL;
	}else{	
		T->index = (2 * T->index + 1);
		//We check if the node has been "created" which will given it a position different from -1
		if((int)T->tree[T->index][T->d+1] ==  -1){
			return NULL;
		}else{
			return T;
		}
	}
}


/*return the outer subtree with all its contents: vpoint, median, index.
In an array representation of a tree the left(inside) child of the Node i is the (2 * i + 2) position.
We also check if the index gets out of boundary*/
vptree * getOuter(vptree * T){
	unsigned max_size = (int)pow(2, ((int)log2(T->n) + 1)) - 1;
	if(max_size <= (2 * T->index + 2)){
		return NULL;
	}else{	
		T->index = (2 * T->index + 2);
		//We check if the node has been "created" which will given it a position different from -1
		if((int)T->tree[T->index][T->d+1] ==  -1){
			return NULL;
		}else{
			return T;
		}
	}
}

/*Returns the median of the subtree in the position index*/
double getMD(vptree * T){
	return (T->tree[T->index][T->d]);
}

/*Returns a pointer to the v point of the subtree in the position index*/
double *getVP(vptree * T){
	return (T->tree[T->index]);
}

int getIDX(vptree * T){
	return (T->index);
}

void createTree(vptree *T, double **arraySet, double *distances) {
	unsigned int index_id = 0;
	createRecursivly(T, index_id, arraySet, distances, 0, T->n );
}


void createRecursivly(vptree *T, int index, double **arrayset, double *distances, int start, int end) {
	unsigned int size = end - start;
	if(size <= 0){
		return;
	}else{
		/*Select the last point of the set as vantage point use memcpy for effiecient copying*/
		T->tree[index][T->d + 1] = index;  
		memcpy(T->tree[index], arrayset[end - 1], T->d * sizeof(double)); 
		int mean = (end - start - 1) / 2;
		
		//to add in parallel
		/*calculate distances between all points and  the vantage point*/
		for(unsigned i = start; i < (end - 1); ++i){
			getDistanceSquared(arrayset[i], T->tree[index], &distances[i], T->d);
		}
		
		/*find the median*/
		int k = (size - 1)/2;
		if(k == 0)
			k = 1;
		T->tree[index][T->d] = kthSmallest(arrayset, distances, start, end - 2, k);
		
		/*Check how many threads exist if there are enough pass the data to another thread
		to continue the recursive creation in parallel*/
		pthread_mutex_lock (&mutexId);
		if(curr_thread < NUM_OF_THREAD){
			/*create a structure to pass the necessary data*/
			TreeData *leftSubTree = (TreeData *)malloc(sizeof(TreeData));
			leftSubTree->T = T;
			leftSubTree->index = (2 * index + 1);
			leftSubTree->arrayset = arrayset;
			leftSubTree->distances = distances;
			leftSubTree->start = start;
			leftSubTree->end = start + mean;
			
			pthread_create(&thread[curr_thread], &attr, recursiveJob, (void *)leftSubTree);
			curr_thread = curr_thread + 1;
			pthread_mutex_unlock (&mutexId);
		}else{
			/* If there are not enough threads continue in serial*/
			pthread_mutex_unlock (&mutexId);
			createRecursivly(T, (2 * index + 1), arrayset, distances, start, start + mean);
		}
		
		/*caller thread takes over the other subtree (outside)*/
		createRecursivly(T, (2 * index + 2), arrayset, distances, start + mean , (end - 1));		
	}
}

/* Work as connector between the recursive function and the created thread 
by passing all the appropriate arguments*/
void *recursiveJob(void *data){
	TreeData *d = (TreeData *)data;
	createRecursivly(d->T, d->index, d->arrayset, d->distances, d->start, d->end);
	//free(d);
	pthread_exit(NULL);
}

/*Swaps two pointers*/
void swap(double **n1, double **n2) {
	double *temp = *n1;
	*n1 = *n2;
	*n2 = temp;
}

/* Swaps two element of a one dimentional array (array)*/
void swap2(double *array, int i, int j) {
	double temp = array[i];
	array[i] = array[j];
	array[j] = temp;
}

/* partition the distance array and also swaps the points of the dataset accordingly*/
double partition(double **arraySet, double *distances, int low, int high) 
{ 
    int pivot = distances[high]; 
    int i = (low - 1); 
    for (int j = low; j <= high - 1; j++) { 
        if (distances[j] <= pivot) { 
            i++; 
            swap2(distances, i, j); 
			swap(&arraySet[i], &arraySet[j]); 
        } 
    } 
    swap2(distances, (i+1), high); 
	swap(&arraySet[i + 1], &arraySet[high]); 
    return (i + 1); 
} 
 
/* Find the K largest element and its used as a median finder
In this implementation there is a one to one relation between distances 
and points of data, so it hanldes both both swaps simultaneously
Lastly special care must be given to the indexes because the array may start 
at any index due to the recursive calls*/ 
double kthSmallest(double **arraySet, double *distances, int left, int right, int k) 
{ 
	int temp = left;
    while (left <= right) { 
        int pivotIndex = partition(arraySet, distances, left, right); 
		//printf("inside k largest %i %i\n",pivotIndex,k-1);
        if ((pivotIndex - temp)  == k - 1) 
            return sqrt(distances[pivotIndex]); 
        else if ((pivotIndex - temp) > k - 1) 
            right = pivotIndex - 1; 
        else
            left = pivotIndex + 1; 
    } 
    return 0.0; 
} 

/* Given the pointer of the two points, where to store the result 
and their dimentions it calculates the Euclidean distance between them*/
void getDistanceSquared(double *p1, double *p2, double *result,int dim)
{
	(*result) = 0;
	for (unsigned int i = 0;i < dim; ++i) {
		(*result) = (*result) + pow((p1[i] - p2[i]), 2);
	}
}

