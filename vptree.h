#ifndef VPTREE_H
#define VPTREE_H


/* Definition of the tree structure in an array representaion form
n: number of nodes
d: dimentions of a point
index: identifier of a subtree
tree: 2 dimentinal array [n x (d + 2)] where tree[i][d] is the radius and tree[i][d+1] position of the node in the tree
*/
typedef struct{
	double **tree; 
	unsigned index;
	unsigned d;
	unsigned n;
}vptree;

// ========== LIST OF ACCESSORS ====================


//! Build vantage-point tree given input dataset X
/*!
\param X Input data points, stored as [n-by-d] array
\param n Number of data points (rows of X)
\param d Number of dimensions (columns of X)
\return The vantage-point tree
*/
vptree * buildvp(double *X, int n, int d);


//! Return vantage-point subtree with points inside radius
/*!
\param node A vantage-point tree
\return The vantage-point subtree
*/
vptree * getInner(vptree * T);


//! Return vantage-point subtree with points outside radius
/*!
\param node A vantage-point tree
\return The vantage-point subtree
*/
vptree * getOuter(vptree * T);


//! Return median of distances to vantage point
/*!
\param node A vantage-point tree
\return The median distance
*/
double getMD(vptree * T);


//! Return the coordinates of the vantage point
/*!
\param node A vantage-point tree
\return The coordinates [d-dimensional vector]
*/
double * getVP(vptree * T);



//! Return the index of the vantage point
/*!
\param node A vantage-point tree
\return The index to the input vector of data points
*/
int getIDX(vptree * T);


#endif