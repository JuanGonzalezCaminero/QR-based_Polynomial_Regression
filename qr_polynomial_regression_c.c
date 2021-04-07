#include "mex.h"
#include "lapack.h"
#include "blas.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <lapack.h>


#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))

void qr_polynomial_regression_c(double *x, double *y, size_t k, size_t n, double *out)
{
	
	double *A, *WORK, *R, *Q, *aux_vec, *aux_mat, WORK_AUX, ALPHA=1, BETA=1;
	double TAU[min(k, n+1)];
	size_t LWORK, INFO, n_plus_one=n+1, INCX=1, INCY=1;
	char UPLO, TRANS, DIAG;

	///////////////////////////////////////////////////////
	//				  Memory Allocation		  			 //
	///////////////////////////////////////////////////////

	A = (double*) malloc(k*(n+1)*sizeof(double));
	R = (double*)malloc((n+1)*(n+1)*sizeof(double));
	Q = (double*)malloc(k*k*sizeof(double));
	aux_vec = (double*) malloc(k*sizeof(double));

	///////////////////////////////////////////////////////
	//			    Memory Initialization		  		 //
	///////////////////////////////////////////////////////

	for(int i=0; i<(n+1)*(n+1); i++)
		R[i]=0;
	for(int i=0; i<k*k; i++)
		Q[i]=0;
	for(int i=0; i<k; i++)
		aux_vec[i]=0;

	///////////////////////////////////////////////////////
	//			 Fill A in column-major order		  	 //
	///////////////////////////////////////////////////////
	
	for(int j=0; j<n+1; j++)
	{
		for(int i=0; i<k; i++)
		{	
			//row i, column j: x[i]^(n)-j 
			//(The higest power is the leftmost element of the row)
			A[j*k + i]=pow(x[i], (n)-j);
		}
	}

	///////////////////////////////////////////////////////
	//				Compute QR using dgeqrf				 //
	///////////////////////////////////////////////////////

	//Workspace query to determine the optimal size of WORK
	LWORK=-1;
	dgeqrf(&k, &n_plus_one, A, &k, TAU, &WORK_AUX, &LWORK, &INFO);

	LWORK=(size_t)WORK_AUX;
	WORK=(double*)malloc((size_t)WORK_AUX*sizeof(double));

	//Returns R avove the diagonal of A, returns Q as a product of householder reflectors
	dgeqrf(&k, &n_plus_one, A, &k, TAU, WORK, &LWORK, &INFO);
	free(WORK);

	///////////////////////////////////////////////////////
	//				   Extract R from A			  	     //
	///////////////////////////////////////////////////////

	UPLO='U';
	dlacpy(&UPLO, &k, &n_plus_one, A, &k, R, &n_plus_one);

	///////////////////////////////////////////////////////
	//		Generate the matrix representation of Q		 //
	///////////////////////////////////////////////////////

	size_t aux_min = min(k, n+1);

	//We need the lower triangular matrix of A, each column 
	//is one of the vectors needed to compute Q

	//Copy the lower triangular matrix of A in Q
	UPLO='L';
	dlacpy(&UPLO, &k, &n_plus_one, A, &k, Q, &k);

	//Workspace query to determine the optimal size of WORK
	LWORK=-1;
	dorgqr(&k, &k, &aux_min, Q, &k, TAU, &WORK_AUX, &LWORK, &INFO);

	LWORK=(size_t)WORK_AUX;
	WORK=(double*)malloc((size_t)WORK_AUX*sizeof(double));

	//Finally, compute Q
	dorgqr(&k, &k, &aux_min, Q, &k, TAU, WORK, &LWORK, &INFO);

	///////////////////////////////////////////////////////
	//				    Compute Q'*b					 //
	///////////////////////////////////////////////////////

	//Q'*b
	TRANS='T';
	dgemv(&TRANS, &k, &k, &ALPHA, Q, &k, y, &INCX, &BETA, aux_vec, &INCY);

	///////////////////////////////////////////////////////
	//				   Solve Rx = Q'*b					 //
	///////////////////////////////////////////////////////

	UPLO='U';
	TRANS='N';
	DIAG='N';
	dtrsv(&UPLO, &TRANS, &DIAG, &n_plus_one, R, &n_plus_one, aux_vec, &INCX);

	///////////////////////////////////////////////////////
	//					Print results					 //
	///////////////////////////////////////////////////////

	printf("R =\n");
	for(int i=0; i<n+1; i++)
	{
		for(int j=0; j<n+1; j++)
		{
			printf("%10.4lf", R[j*(n+1) + i]);
		}
		printf("\n");
	}
	printf("\n");
	
	printf("Q =\n");
	for(int i=0; i<k; i++)
	{
		for(int j=0; j<k; j++)
		{
			printf("%10.4lf", Q[j*k + i]);
		}
		printf("\n");
	}
	printf("\n");

	printf("Rx = Q'b\n");
	for(int i=0; i<n+1; i++)
	{
		printf("%10.4lf", aux_vec[i]);
		printf("\n");
	}
	printf("\n");

	///////////////////////////////////////////////////////
	//           Copy result to output array			 //
	///////////////////////////////////////////////////////

	for(int i=0; i<n+1; i++)
	{
		out[i]=aux_vec[i];
	}

	///////////////////////////////////////////////////////
	//					Free memory						 //
	///////////////////////////////////////////////////////

	free(WORK);
	free(aux_vec);
	free(A);
	free(R);
	free(Q);
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

	//Coordinates of the points
	double *x, *y;
	//Output
	double *out;
	//Degree of the desired polynomial
	size_t n;
	//Number of points
	size_t k;

	//Check for proper number of arguments.
	if(nrhs!=3) 
	{
		mexErrMsgIdAndTxt( "MATLAB:qr_polynomial_regression_c:invalidNumInputs", "One input required.");
	}
	else if(nlhs>3) 
	{
		mexErrMsgIdAndTxt( "MATLAB:qr_polynomial_regression_c:maxlhs", "Too many output arguments.");
	}
	//Check the type of the arguments
	if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
	{
		mexErrMsgIdAndTxt("MATLAB:qr_polynomial_regression_c:inputNotRealScalarDouble", "x must be a noncomplex scalar double.");
	}
	if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]))
	{
		mexErrMsgIdAndTxt("MATLAB:qr_polynomial_regression_c:inputNotRealScalarDouble", "y must be a noncomplex scalar double.");
	}

	n = (size_t)((double*)mxGetPr(prhs[2]))[0];

	//Allocate the output array
	plhs[0] = mxCreateDoubleMatrix(n+1, 1, mxREAL);

	//Assign the inputs and otputs to C pointers
	x = mxGetPr(prhs[0]);
	y = mxGetPr(prhs[1]);
	out = mxGetPr(plhs[0]);

	//Get the number of points
	k = mxGetM(prhs[0]);

	qr_polynomial_regression_c(x, y, k, n, out);
}