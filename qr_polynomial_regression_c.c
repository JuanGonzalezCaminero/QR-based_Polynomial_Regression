#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "lapack.h"

#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))

void qr_polynomial_regression_c(double *x, double *y, size_t k, size_t n)
{
	//Matrix A
	double *A;

	A = (double*) malloc(k*n*sizeof(double));

	//Fill A in column-major order
	for(int j=0; j<n+1; j++)
	{
		for(int i=0; i<k; i++)
		{	
			//row i, column j: x[i]^(n)-j (The higest power is the leftmost element of the row)
			A[j*k + i]=pow(x[i], (n)-j);
		}
	}

	double TAU[min(k, n+1)];
	double WORK_AUX[1];
	size_t LWORK, INFO;
	size_t n_plus_one=n+1;

	//Workspace query to determine the optimal size of WORK
	LWORK=-1;
	dgeqrf(&k, &n_plus_one, A, &k, TAU, WORK_AUX, &LWORK, &INFO);

	LWORK=(size_t)WORK_AUX[0];
	double *WORK=(double*)malloc((size_t)WORK_AUX[0]*sizeof(double));

	//Returns R avove the diagonal of A, returns Q as a product of householder reflectors
	dgeqrf(&k, &n_plus_one, A, &k, TAU, WORK, &LWORK, &INFO);

	//Generate the matrix representation of Q:

	size_t aux_min = min(k, n+1);

	//We need the lower triangular matrix of A, each column is one of the vectors needed to compute Q
	//Q will be square, but its dimensions depend on the dimensions of A. (Specifically the number of rows)
	double *Q;
	Q = (double*)malloc(k*k*sizeof(double));

	//Copy the lower triangular matrix of A in Q
	char UPLO='L';
	dlacpy(&UPLO, &k, &n_plus_one, A, &k, Q, &k);

	//Workspace query to determine the optimal size of WORK
	LWORK=-1;
	dorgqr(&k, &k, &aux_min, Q, &k, TAU, WORK_AUX, &LWORK, &INFO);

	LWORK=(size_t)WORK_AUX[0];
	free(WORK);
	WORK=(double*)malloc((size_t)WORK_AUX[0]*sizeof(double));

	//Finally, compute Q
	dorgqr(&k, &k, &aux_min, Q, &k, TAU, WORK, &LWORK, &INFO);

	//Extract R from A
	double *R;
	R = (double*)malloc((n+1)*(n+1)*sizeof(double));
	UPLO='U';
	dlacpy(&UPLO, &k, &n_plus_one, A, &k, R, &n_plus_one);

	for(int i=0; i<n+1; i++)
	{
		for(int j=0; j<n+1; j++)
		{
			printf("%lf, ", R[j*(n+1) + i]);
		}
		printf("\n");
	}
	printf("\n");

	for(int i=0; i<k; i++)
	{
		for(int j=0; j<k; j++)
		{
			printf("%lf, ", Q[j*k + i]);
		}
		printf("\n");
	}
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

	/* Check for proper number of arguments. */
	if(nrhs!=3) 
	{
		mexErrMsgIdAndTxt( "MATLAB:qr_polynomial_regression_c:invalidNumInputs", "One input required.");
	}
	else if(nlhs>3) 
	{
		mexErrMsgIdAndTxt( "MATLAB:qr_polynomial_regression_c:maxlhs", "Too many output arguments.");
	}
	/* Check the type of the arguments */
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

	//Assign the inputs to C pointers
	x = mxGetPr(prhs[0]);
	y = mxGetPr(prhs[1]);

	//Get the number of points
	k = mxGetM(prhs[0]);

	qr_polynomial_regression_c(x, y, k, n);
}