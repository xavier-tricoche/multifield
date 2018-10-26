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
 *  Paul's GSL wrappers for the GNU Linear Programming Kit
 *
 *	Author: Paul O'Grady
 *	www.hamilton.ie/paul
 *
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
 *
 */

#include "GLPK_GSL_wrap.h"


/* L1-norm minimisation (minimization) */

/* Set up L1-norm minimisation using GLPK solver */
LPX  *setUpLOneNormMinimisation(gsl_matrix * basis)
{

	int basisRows = basis->size1;
	int basisCols = basis->size2;
	int row,col, basisSize;

	basisSize = (basisRows * basisCols)*2;/* [A -A] */

	double point;

	LPX *lp;
	/*int count, dataPointCount, rn[1+basisSize], cn[1+basisSize];
	double a[1+basisSize], Z, coeff[basisCols*2], matrixValue;*/
	int count, dataPointCount;
	int* rn = (int *) malloc((1+basisSize)*sizeof(int));
	int* cn = (int *) malloc((1+basisSize)*sizeof(int));
	double Z, matrixValue;
	double* a = (double *) malloc((1+basisSize)*sizeof(double));
	double* coeff = (double *) malloc(basisCols*sizeof(double));

	/* Check memory allocation */
	if ((rn == NULL)|(cn == NULL)|(a == NULL)|(coeff == NULL))
	{
		/* Memory could not be allocated, so print an error and exit. */
		fprintf(stderr, "Couldn't allocate memory\n");
		exit(EXIT_FAILURE);
	}

	lp = lpx_create_prob();

	/* Switch off display */
	lpx_set_int_parm(lp,LPX_K_MSGLEV,0);

	lpx_set_prob_name(lp, "L1min");

	lpx_add_rows(lp, basisRows);
	lpx_add_cols(lp, (basisCols*2));

	/* fix lower bound of coefficients to 0 */
	for(count = 0; count < (basisCols*2); count++)
	{
		lpx_set_col_bnds(lp,(count+1), LPX_LO, 0.0, 0.0);
	}

	/* [A] */
	count = 1; /* GLPK doesn't use 0 index */
	for(row = 0; row < basisRows; row++)
	{
		for(col = 0; col< basisCols; col++)
		{

			matrixValue = gsl_matrix_get(basis,row,col);

			if(matrixValue != 0) /* Only load non-zero elements */
			{
				rn[count] = (row+1);
				cn[count] = (col+1);
				a[count] = matrixValue;
				count++;
			}
			else
				printf("zero element\n");
		}
	}

	/* [A -A] */
	for(row = 0; row < basisRows; row++)
	{
		for(col = basisCols; col< (basisCols*2); col++)
		{

			matrixValue = -1*(gsl_matrix_get(basis,row,(col-basisCols)));

			if(matrixValue != 0)/* Only load non-zero elements */
			{
				rn[count] = (row+1);
				cn[count] = (col+1);
				a[count] = matrixValue;
				count++;
			}
			else
				printf("zero element\n");
		}
	}

	/* lpx_load_mat3(lp, basisSize, rn, cn, a); */
	lpx_load_matrix(lp, (count-1), rn, cn, a);

	lpx_set_obj_dir(lp, LPX_MIN);

	for (count = 0; count < (basisCols*2); count++)
	{
		/* lpx_set_col_coef(lp, (count+1), 1.0); */
		lpx_set_obj_coef(lp, (count+1), 1.0);
	}

	//free(rn);
	//free(cn);
	//free(a);
	//free(coeff);

	return (lp);

}

/* Perform L1-norm minimisation on a data point */
void LOneNormMinimisationVector(LPX *lp, gsl_vector * dataVector, gsl_vector *  minL1NormSolution)
{

	int lengthOfBasisVectors = dataVector->size;
	int noOfBasisVectors = minL1NormSolution->size;
	int count;
	double* coeff = (double*) malloc(noOfBasisVectors*2 * sizeof(double));


	/* Perform Simplex Algorithm */
	if(gsl_vector_isnull(dataVector))
	{
		gsl_vector_set_zero(minL1NormSolution);
	}
	else
	{
		for(count = 0; count < lengthOfBasisVectors; count++)
		{
			lpx_set_row_bnds(lp, (count+1), LPX_FX, dataVector->data[count], 0.0);
		}

		//lpx_simplex(lp);
		if (lpx_simplex(lp) != LPX_E_OK)
		{
			printf(">>>>LOneNormMinimisationVector: Problem has not been sucessfully solved\n");
		}

		for (count = 0; count < (noOfBasisVectors*2); count++)
		{
			lpx_get_col_info(lp, (count+1), NULL, &coeff[count], NULL);
		}

		for(count = 0; count < noOfBasisVectors; count++)
		{
			gsl_vector_set(minL1NormSolution,count,(coeff[count] - coeff[count+(noOfBasisVectors)]));
		}
	}

	//free(coeff);
}


/* Perform L1-norm minimisation on a matrix of data */
void LOneNormMinimisationMatrix(gsl_matrix* basis, gsl_matrix* data, gsl_matrix* minL1NormSolutions)
{

	int length = data->size2;
	int dataPointLength = data->size1;
	int solutionLength = minL1NormSolutions->size1;
	int count;

	gsl_vector* minL1NormSolution = gsl_vector_alloc(solutionLength);
	gsl_vector* dataPoint = gsl_vector_alloc(dataPointLength);

	/* Set up L1-norm minimisation */
	LPX *lp;
	lp = setUpLOneNormMinimisation(basis);

	/* Perform L1-norm minimisation */
	for(count = 0; count < length; count++)
	{
		gsl_matrix_get_col(dataPoint,data, count);
		LOneNormMinimisationVector(lp,dataPoint,minL1NormSolution);
		gsl_matrix_set_col(minL1NormSolutions,count,minL1NormSolution);
	}

	/* Free up allocated memory */
	gsl_vector_free(dataPoint);
	gsl_vector_free(minL1NormSolution);
	//lpx_delete_prob(lp);

}
