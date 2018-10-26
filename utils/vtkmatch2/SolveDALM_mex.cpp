/*************************************************************************
multifield:   Surface-based Structural Analysis and Visualization 
              of Multifield Datasets

Author: Samer S. Barakat

Copyright (c) 2010-2012, Purdue University

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**************************************************************************/
/*
% MEX/C code for DALM l1-minimization

% Copyright ©2010. The Regents of the University of California (Regents).
% All Rights Reserved. Contact The Office of Technology Licensing,
% UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620,
% (510) 643-7201, for commercial licensing opportunities.

% Created by Victor Shia, Mark Murphy, Allen Y. Yang, Department of EECS, University of California,
% Berkeley.

% IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
% SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
% ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
% REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY,
% PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO
% PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <f2c.h>
#include <clapack.h>

#include "SolveDALM.h"

/*
int max(int a, int b)
{
	if (a > b) return a;
	return b;
}*/

double const eps = 1e-15;

//#if !defined (COMPILE_MEX)
//# include <stdio.h>
//#else
//# include <mex.h>
//# undef printf
//# define printf printf
//
//extern "C" void
//mexFunction (int nl, mxArray *pl[], int nr, mxArray const *pr[])
//{
//	if (nr < 4){
//	  mexErrMsgTxt ("[x nIter] = SolveDALM(A, b, nu, tol, stop, xG)");
//	}
//
//
//	double *b = mxGetPr (pr[1]);
//	double *A = mxGetPr (pr[0]);
//	int m = mxGetM (pr[0]);
//	int n = mxGetN (pr[0]);
//	
//
//	if (mxGetM (pr[1]) * mxGetN (pr[1])  != m){
//		mexErrMsgTxt ("SolveDALM: min |x|1 + |e|1 s.t. Ax + e = b\n");
//	}
//
//	double nu = mxGetScalar (pr[2]);
//	double tol = mxGetScalar (pr[3]);
//	double *xG;	int stop;
//
//	if (nr < 6)
//	  xG = NULL;
//	else
//	  xG = mxGetPr(pr[5]);
//	if(nr < 5)
//	  stop = 5;
//	else
//	  stop = (int)mxGetScalar(pr[4]);
//	
//	double *x = new double[n];
//	int nIter;
//	int maxIter = 5000;
//
//	SolveDALM(x, nIter, b, A, nu, tol, maxIter, m, n, stop, xG);
//	
//	if (nl > 0){
//		pl[0] = mxCreateNumericMatrix (n, 1, mxDOUBLE_CLASS, mxREAL);
//		memcpy (mxGetData (pl[0]), (void*)x, n*sizeof(double));
//	}
//	delete [] x;
//
//	if (nl > 1){
//		pl[1] = mxCreateDoubleScalar (nIter);
//	}
//}
//
//#endif

extern "C"
{
	int  MAIN__(int argc, char** argv) {  printf("%d %d\n", argc, __LINE__); return 0; };
	
/* Subroutine */ int dgemm_(char *transa, char *transb, integer *m, integer *
	n, integer *k, doublereal *alpha, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, 
	integer *ldc);
/* Subroutine */ int dgetrf_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info);
/* Subroutine */ int dgetrf_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info);
/* Subroutine */ int dgetri_(integer *n, doublereal *a, integer *lda, integer 
	*ipiv, doublereal *work, integer *lwork, integer *info);
/* Subroutine */ int dgemv_(char *trans, integer *m, integer *n, doublereal *
	alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, 
	doublereal *beta, doublereal *y, integer *incy);
}

enum stoppingCriteria{
  STOPPING_GROUND_TRUTH   = -1,
  STOPPING_DUALITY_GAP    = 1,
  STOPPING_SPARSE_SUPPORT = 2,
  STOPPING_OBJECTIVE_VALUE = 3,
  STOPPING_SUBGRADIENT    = 4,
  STOPPING_INCREMENTS     = 5,
  STOPPING_DEFAULT        = STOPPING_INCREMENTS
};

int dgemm(char transa, char transb, long int m, long int n, long int k, double alpha, double *a, long int lda, double *b, long int ldb, double beta, double *c__, long int ldc)
{
	return dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c__, &ldc);
}

int dgetrf(long int m, long int n, double *a, long int lda, long int *ipiv, long int& info)
{
	return dgetrf_(&m, &n, a, &lda, ipiv, &info);
}

int dgetri(long int m, double *a, long int lda, long int *ipiv, long int& info)
{
	long int lwork = 2 * m;
	double * work = (double*) malloc(2 * m * sizeof(double));
	int ret = dgetri_(&m, a, &lda, ipiv, work, &lwork, &info);
	free(work);
	return ret;
}

int dgemv(char trans, long int m, long int n, doublereal alpha, double *a, long int lda, double *x, long int incx, double beta, double *y, long int incy)
{
	return dgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

void SolveDALM (double *x, int&  nIter, double *b, double *A, double nu, double tol, int maxIter, int m, int n, int stoppingCriterion, double *xG)
{
	int ldA = m;

	enum stoppingCriteria stop;
	switch (stoppingCriterion){
		case -1:
		stop = STOPPING_GROUND_TRUTH;
		break;
		case 1:
		stop = STOPPING_DUALITY_GAP;
		break;
		case 2:
		stop = STOPPING_SPARSE_SUPPORT;
		break;
		case 3:
		stop = STOPPING_OBJECTIVE_VALUE;
		break;
		case 4:
		stop = STOPPING_SUBGRADIENT;
		break;
		case 5:
		stop = STOPPING_INCREMENTS;
		break;
	}
	
	bool verbose = false;

// beta = norm(b,1) / m;
	double beta = 0;
	for (long i = 0; i < m; i++){
		beta += fabs (b[i]);
	}
	beta = beta / m;

	double betaInv = 1 / beta;

	nIter = 0;

// G = A * A' + eye(m) * lambda (or nu) / beta;
	double *G = new double[m*m];
	int ldG = m;
	double tempInt = nu/beta;

	dgemm ('N', 'T', m, m, n, 1.0, A, ldA, A, ldA, 0.0, G, ldG); 

	for (long i = 0; i < m; i++)
		G[i + ldG*i] += tempInt;

// invG = inv (G)
	double *invG = new double[m*m];
	int ldinvG = m;
	long *ipiv = new long[m], info = 0;
	for (long i = 0; i < m*m; i++)
		invG[i] = G[i];

	dgetrf (m, m, invG, ldinvG, ipiv, info);
	dgetri (m, invG, ldinvG, ipiv, info);

	delete[] ipiv;

// A_invG_b = A' * invG * b
	double *A_invG_b = new double [n];
	double *tmp = new double [m];

	dgemv ('N', m, m, 1.0, invG, ldinvG, b, 1, 0.0, tmp, 1);
	dgemv ('T', m, n, 1.0, A, ldA, tmp, 1, 0.0, A_invG_b, 1);

	delete [] tmp;

// y = zeros(m,1)
	double *y = new double[m];
	for (long i = 0; i < m; i++)
		y[i] = 0;

// x = zeros(n,1)
//	x = new double[n];
	for (long i = 0; i < n; i++)
		x[i] = 0;

// z = zeros (m+n,1);
	double *z = new double[n]; 
	for (long i = 0; i < n; i++)
		z[i] = 0;

	bool converged_main = false;

// temp = A' * y;
	double *temp = new double [max(m,n)];
	dgemv ('T', m, n, 1.0, A, ldA, y, 1, 0.0, temp, 1);
 
	double *x_old  = new double[n];
	double *temp1 = new double[max(m,n)];
	tmp = new double[max(m,n)];

// f = norm(x,1);  x is 0 at this point
	double f = 0;
	double prev_f = 0;
	double total = 0;
	double nxo, nx, dx;
	
	do {
		nIter = nIter + 1;
		if(verbose) printf("==== [%d] ====\n", nIter);

	// x_old = x
		for (long i = 0; i < n; i++)
			x_old[i] = x[i];

	// % update z
	// temp1 = temp + x * betaInv
	// z = sign(temp1) .* min(1, abs(temp1));
		for (long i = 0; i < n; i++){
			temp1[i] = temp[i] + x[i] * betaInv;
		}
		for (long i = 0; i < n; i++){
			z[i] = (temp1[i] > 0 ? 1 : -1)
			       * ((fabs(temp1[i]) > 1) ? 1 : fabs(temp1[i]));
		}

	// temp = A' * (invG * (A * (z - xv * betaInv))) + A_invG_b * betaInv
		for (long i = 0; i < n; i++)
			temp1[i] = z[i] - x[i]/beta;

		dgemv ('N', m, n, 1.0, A, ldA, temp1, 1, 0.0, tmp, 1);
		dgemv ('N', m, m, 1.0, invG, ldinvG, tmp, 1, 0.0, temp1, 1);
		dgemv ('T', m, n, 1.0, A, ldA, temp1, 1, 0.0, tmp, 1);

		for (long i = 0; i < n; i++)
			temp[i] = tmp[i] + A_invG_b[i]/beta;

	// % update x
	// x = x - beta * (z - temp);
		for (long i = 0; i < n; i++)
			x[i] = x[i] - beta * (z[i] - temp[i]);

		switch (stop){
		case STOPPING_GROUND_TRUTH:
		  total = 0;
		  for(int i = 0 ; i < n; i++){
		    total += (xG[i] - x[i])*(xG[i] - x[i]);
		  }
		
		  if (total < tol * tol)
		    converged_main = true;
		  break;
		case STOPPING_SUBGRADIENT:
		  printf("Duality gap is not a valid stopping criterion for ALM.");
		  break;
		case STOPPING_SPARSE_SUPPORT:
		  printf("DALM does not have a support set.");
		  break;
		case STOPPING_OBJECTIVE_VALUE:
		  prev_f = f;
		  f = 0;
		  for(int i = 0 ; i < n; i++){
		    f += fabs(x[i]);
		  }
		  if (fabs(f-prev_f)/prev_f <= tol){
		    converged_main = true;
		  }
		  break;
		case STOPPING_DUALITY_GAP:
		  printf("Duality gap is not a valid stopping criterion for ALM.");
		  break;
		case STOPPING_INCREMENTS:
		  // if norm(x_old - x) < tol * norm(x_old)
		  //     converged_main = true;
		  
		  nxo = 0;
		  for (int i = 0; i < n; i++)
		    nxo = nxo + x_old[i]*x_old[i];
		  nxo = sqrt (nxo);
		  
		  nx = 0;
		  for (int i = 0; i < n; i++)
		    nx = nx + x[i]*x[i];
		  nx = sqrt (nx);
		  
		  dx = 0;
		  for (int i = 0; i < n; i++)
		    dx = dx + (x_old[i]-x[i])*(x_old[i]-x[i]);
		  dx = sqrt(dx);
		  
		  if (dx < tol*nxo)
		    converged_main = true;
		  
		  if (verbose){
		    printf("  ||x|| = %f\n", nx);
		  }
		  
		  if (verbose){
		    if (nIter > 1){
		      printf ("  ||dx|| = %f (= %f * ||x_old||)\n",
				 dx, dx/(nxo+eps));
		    } else {
		      printf ("  ||dx|| = %f\n", dx);
		    }
		  }
		  break;
		default:
		  printf("Undefined stopping criterion.");
		  break;
		}

 		if (nIter >= maxIter){
			if (verbose)
				printf ("Maximum Iterations Reached\n");
			converged_main = true;
		}
	}
	while (!converged_main);

	if (verbose) printf("==== CONVERGED ==== \n", nIter);

	delete [] G;
	delete [] invG;
	delete [] A_invG_b;
	delete [] tmp;
	delete [] y;
	delete [] z;
	delete [] x_old;
	delete [] temp;
	delete [] temp1;
}

void SolveDALM_fast (double *&x, int&  nIter, double *b, double *A, double lambda, double tol, int maxIter, int m, int n, int stoppingCriterion, double *xG)
{
	int ldA = m;

	enum stoppingCriteria stop;
	switch (stoppingCriterion){
		case -1:
		stop = STOPPING_GROUND_TRUTH;
		break;
		case 1:
		stop = STOPPING_DUALITY_GAP;
		break;
		case 2:
		stop = STOPPING_SPARSE_SUPPORT;
		break;
		case 3:
		stop = STOPPING_OBJECTIVE_VALUE;
		break;
		case 4:
		stop = STOPPING_SUBGRADIENT;
		break;
		case 5:
		stop = STOPPING_INCREMENTS;
		break;
	}

	bool verbose = false;

// beta = norm(b,1) / m;
	double beta = 0;
	for (long i = 0; i < m; i++){
		beta += fabs (b[i]);
	}
	beta = beta / m;

	double betaInv = 1 / beta;

	nIter = 0;

// y = zeros(m,1)
	double *y = new double[m];
	for (long i = 0; i < m; i++)
		y[i] = 0;

// x = zeros(n,1)
	//x = new double[n];
	for (long i = 0; i < n; i++)
		x[i] = 0;

// z = zeros (m+n,1);
	double *z = new double[n]; 
	for (long i = 0; i < n; i++)
		z[i] = 0;

	bool converged_main = false;

// temp = A' * y;
	double *temp = new double [max(m,n)];
	dgemv ('T', m, n, 1.0, A, ldA, y, 1, 0.0, temp, 1);
 
	double *x_old  = new double[n];
	double *temp1 = new double[max(m,n)];
	double *tmp = new double[max(m,n)];
	double *g = new double[m];

// f = norm(x,1);  x is 0 at this point
	double f = 0;
	double prev_f = 0;
	double total = 0;
	double nxo, nx, dx;
	
	double *Ag = new double[n];
	
	double dg, dAg, alpha;
	
	do {
		nIter = nIter + 1;
		if(verbose) printf("==== [%d] ====\n", nIter);

	// x_old = x
		for (long i = 0; i < n; i++)
			x_old[i] = x[i];

	// % update z
	// temp1 = temp + x * betaInv
	// z = sign(temp1) .* min(1, abs(temp1));
		for (long i = 0; i < n; i++){
			temp1[i] = temp[i] + x[i] * betaInv;
		}
		for (long i = 0; i < n; i++){
			z[i] = (temp1[i] > 0 ? 1 : -1)
			       * ((fabs(temp1[i]) > 1) ? 1 : fabs(temp1[i]));
		}

	//    %compute A' * y    
	//    g = lambda * y - b + A * (beta * (temp - z) + x);
		for(long i = 0 ; i < n; i++){
			tmp[i] = beta * (temp[i] - z[i]) + x[i];
		}
		for(long i = 0; i < m; i++){
			g[i] = lambda * y[i] - b[i];
		}
		dgemv('N', m, n, 1.0, A, ldA, tmp, 1, 1, g, 1);
	
	//    %alpha = g' * g / (g' * G * g);
	//    Ag = A' * g;
		dgemv('T', m, n, 1.0, A, ldA, g, 1, 0.0, Ag, 1);
	
	//    alpha = g' * g / (lambda * g' * g + beta * Ag' * Ag);
		dg = 0;
		dAg = 0;
		for(long i = 0 ; i < n; i++){
			dAg += Ag[i] * Ag[i];
		}
		for(long i = 0 ; i < m; i++){
			dg += g[i] * g[i];
		}
		alpha = dg / (lambda * dg + beta * dAg);
	
	//    y = y - alpha * g;
		for(long i = 0 ; i < m; i++){
			y[i] = y[i] - alpha * g[i];
		}
	
	//    temp = A' * y;
		dgemv ('T', m, n, 1.0, A, ldA, y, 1, 0.0, temp, 1);

	// % update x
	// x = x - beta * (z - temp);
		for (long i = 0; i < n; i++)
			x[i] = x[i] - beta * (z[i] - temp[i]);

		switch (stop){
		case STOPPING_GROUND_TRUTH:
		  total = 0;
		  for(int i = 0 ; i < n; i++){
		    total += (xG[i] - x[i])*(xG[i] - x[i]);
		  }
		  if (total < tol * tol)
		    converged_main = true;
		  break;
		case STOPPING_SUBGRADIENT:
		  printf("Duality gap is not a valid stopping criterion for ALM.");
		  break;
		case STOPPING_SPARSE_SUPPORT:
		  printf("DALM does not have a support set.");
		  break;
		case STOPPING_OBJECTIVE_VALUE:
		  prev_f = f;
		  f = 0;
		  for(int i = 0 ; i < n; i++){
		    f += fabs(x[i]);
		  }
		  if (fabs(f-prev_f)/prev_f <= tol){
		    converged_main = true;
		  }
		  break;
		case STOPPING_DUALITY_GAP:
		  printf("Duality gap is not a valid stopping criterion for ALM.");
		  break;
		case STOPPING_INCREMENTS:
		  // if norm(x_old - x) < tol * norm(x_old)
		  //     converged_main = true;
		  
		  nxo = 0;
		  for (int i = 0; i < n; i++)
		    nxo = nxo + x_old[i]*x_old[i];
		  nxo = sqrt (nxo);
		  
		  nx = 0;
		  for (int i = 0; i < n; i++)
		    nx = nx + x[i]*x[i];
		  nx = sqrt (nx);
		  
		  dx = 0;
		  for (int i = 0; i < n; i++)
		    dx = dx + (x_old[i]-x[i])*(x_old[i]-x[i]);
		  dx = sqrt(dx);
		  
		  if (dx < tol*nxo)
		    converged_main = true;
		  
		  if (verbose){
		    printf("  ||x|| = %f\n", nx);
		  }
		  
		  if (verbose){
		    if (nIter > 1){
		      printf ("  ||dx|| = %f (= %f * ||x_old||)\n",
				 dx, dx/(nxo+eps));
		    } else {
		      printf ("  ||dx|| = %f\n", dx);
		    }
		  }
		  break;
		default:
		  printf("Undefined stopping criterion.");
		  break;
		}

 		if (nIter >= maxIter){
			if (verbose)
				printf ("Maximum Iterations Reached\n");
			converged_main = true;
		}
	}
	while (!converged_main);

	if (verbose) printf("==== CONVERGED ==== \n", nIter);

	delete [] tmp;
	delete [] g;
	delete [] Ag;
	delete [] y;
	delete [] z;
	delete [] x_old;
	delete [] temp;
	delete [] temp1;
}

double SumOver(double*  arr, int n)
{
	double s = 0.0;
	for (int i = 0; i < n; i++)
		s += abs(arr[i]);
	return s;
}

void SetXZero(double*  arr, int n)
{
	for (int i = 0; i < n; i++)
		arr[i] = 0.0;
}

void SolveDALML1(double *A, double *x, double *b, int m, int n)
{
	int maxIter = 50000;
	double tol = 0.001;
	int iter;
	double lambda = 0.01;
	double *xG = new double[n];

	double* rA = new double[m * n];
	double* rb = new double[m * n];

	memcpy(rA, A, m * n * sizeof(double));
	for (int i = 0; i < m; i++)
		rb[i] = b[i * n + 0];
	if (SumOver(rb, m) == 0.0)
		SetXZero(x + 0 * n, n);
	else
		SolveDALM(x + 0 * n, iter, rb, rA, lambda, tol, maxIter, m, n, 5, xG);

	memcpy(rA, A, m * n * sizeof(double));
	for (int i = 0; i < m; i++)
		rb[i] = b[i * n + 1];
	if (SumOver(rb, m) == 0.0)
		SetXZero(x + 1 * n, n);
	else
		SolveDALM(x + 1 * n, iter, rb, rA, lambda, tol, maxIter, m, n, 5, xG);

	memcpy(rA, A, m * n * sizeof(double));
	for (int i = 0; i < m; i++)
		rb[i] = b[i * n + 2];
	if (SumOver(rb, m) == 0.0)
		SetXZero(x + 2 * n, n);
	else
		SolveDALM(x + 2 * n, iter, rb, rA, lambda, tol, maxIter, m, n, 5, xG);

	memcpy(rA, A, m * n * sizeof(double));
	for (int i = 0; i < m; i++)
		rb[i] = b[i * n + 3];
	if (SumOver(rb, m) == 0.0)
		SetXZero(x + 3 * n, n);
	else
		SolveDALM(x + 3 * n, iter, rb, rA, lambda, tol, maxIter, m, n, 5, xG);
	

	delete rA;
	delete rb;
	delete xG;
}

void SolveDALML1_fast(double *A, double *x, double *b, int m, int n)
{
	int maxIter = 50000;
	double tol = 0.0000001;
	int iter;
	double lambda = 0.0001;
	double *xG = new double[n];

	double* rA = new double[m * n];
	double* rb = new double[m * n];

	double* xp;

	xp = x + 0 * n;
	memcpy(rA, A, m * n * sizeof(double));
	for (int i = 0; i < m; i++)
		rb[i] = b[i * n + 0];
	if (SumOver(rb, m) == 0.0)
		SetXZero(x + 0 * n, n);
	else
		SolveDALM_fast(xp, iter, rb, rA, lambda, tol, maxIter, m, n, 5, xG);

	xp = x + 1 * n;
	memcpy(rA, A, m * n * sizeof(double));
	for (int i = 0; i < m; i++)
		rb[i] = b[i * n + 1];
	if (SumOver(rb, m) == 0.0)
		SetXZero(x + 1 * n, n);
	else
		SolveDALM_fast(xp, iter, rb, rA, lambda, tol, maxIter, m, n, 5, xG);

	xp = x + 2 * n;
	memcpy(rA, A, m * n * sizeof(double));
	for (int i = 0; i < m; i++)
		rb[i] = b[i * n + 2];
	if (SumOver(rb, m) == 0.0)
		SetXZero(x + 2 * n, n);
	else
		SolveDALM_fast(xp, iter, rb, rA, lambda, tol, maxIter, m, n, 5, xG);

	xp = x + 3 * n;
	memcpy(rA, A, m * n * sizeof(double));
	for (int i = 0; i < m; i++)
		rb[i] = b[i * n + 3];
	if (SumOver(rb, m) == 0.0)
		SetXZero(x + 3 * n, n);
	else
		SolveDALM_fast(xp, iter, rb, rA, lambda, tol, maxIter, m, n, 5, xG);
	

	delete rA;
	delete rb;
	delete xG;
}