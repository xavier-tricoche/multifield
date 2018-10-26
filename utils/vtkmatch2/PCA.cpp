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
///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2003-2005 Neuroshare Project                                                         
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// A copy of the GNU Lesser General Public License can be obtained by writing to:
//  Free Software Foundation, Inc.,
//  59 Temple Place, Suite 330,
//  Boston, MA  02111-1307
//  USA
//
// Contact information:
//  Kirk Korver
//  CyberKinetics, Inc.,
//  391 G Chipeta Way
//  Salt Lake City,  UT  84108
//  USA
//  kkorver@cyberkineticsinc.com
//
// Website:
//  www.neuroshare.org
//
// All other copyrights on this material are replaced by this license agreeement.
//
///////////////////////////////////////////////////////////////////////////////////////////////////
//
// $Author: kirkkorver $
// $Date: 2005/07/11 22:42:39 $
// $Revision: 1.1 $
// $Source: /cvsroot/neuroshare/Suite/PowerNAP/PCA.cpp,v $
//
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "PCA.h"


void meancols(gsl_vector *meanc, gsl_matrix *data);
void sumcols(gsl_vector *sumc, gsl_matrix *data);

// Author & Date:   Kirk Korver     06 Jul 2005
// Purpose: Given the data, center it about 0
// Inputs:
//  pData - the data matrix to adjust the values
void SubtractMeans(gsl_matrix * data)
{
    int m = data->size1;        // # rows   
    int n = data->size2;        // # cols

    gsl_vector * pAvg = gsl_vector_alloc(n);
    meancols(pAvg, data);

    for (int col = 0; col < n; ++col) 
    {
        double dVal = gsl_vector_get(pAvg, col);
        for (int i = 0; i < m; ++i)
        {
            double dNew = gsl_matrix_get(data, i, col) - dVal;
            gsl_matrix_set(data, i, col, dNew);
        }
    }
    gsl_vector_free(pAvg);
}

// Author & Date:   Kirk Korver 07 Jul 2005
// Purpose: Transpose the matrix, even if not square
void PCA::non_square_transpose(gsl_matrix * data)
{
    assert(0);
    // This should work....but I got some weirdness when I used it, so there
    // may be a problem
    gsl_matrix * tmp = gsl_matrix_alloc(data->size2, data->size1);
    gsl_matrix_transpose_memcpy(tmp, data);
    gsl_matrix_free(data);
    *data = *tmp;
}

// Author & Date:   Kirk Korver     07 Jul 2005
//
// Algorithm from 
//   A TUTORIAL ON PRINCIPAL COMPONENT ANALYSIS
//   Derivation, Discussion and Singular Value Decomposition
// by
//   Jon Shlens | jonshlens@ucsd.edu
//
// Purpose: do the principle component analysis by acting directly on
//  the covariance matrix.
// Inputs:
//  data - the data matrix. It will be "clobbered" by this activity (n x m)
//         m = num of dimensions     and     n = num of trials
// Outputs:
//  pc - the eigen vectors used for the conversion (m x m)
//  score - the data remapped into PC space (n x m)
//          NOTICE, that score data data are transposed in dimension from data
void PCA::FindPCA(gsl_matrix *pc, gsl_matrix *score, gsl_matrix *data)
{
    int n = data->size1;        // # trials
    int m = data->size2;        // # dimensions

    //TRACE("N:%d M:%d\n", n, m);

    SubtractMeans(data);            // convert the data, so it is symmetric about 0

    // covariance = 1 / N-1 * data * data' <== MATLAB code
    gsl_matrix * covariance = gsl_matrix_alloc(m, m);
    double dScale = 1.0 / (n - 1.0);

    // Use BLAS to perform the multiplication 
    // C = \alpha * A * B + \beta * C...but  \alpha = 1 / N-1 and \beta = 0    so we get
    // C = 1/(N-1) AB.....but   A could be A or A_transform and  B could be B or B_transform   so....
    //
    // Remember...our "data" is already transformed
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, dScale, data, data, 0.0, covariance);

    gsl_vector * eval = gsl_vector_alloc(m);                        // the eigen values
    {
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(m);        // temp
    gsl_eigen_symmv(covariance, eval, pc, w);
    gsl_eigen_symmv_free(w);
    }

    // Sort the eigen vectors and eigen values
    gsl_eigen_symmv_sort(eval, pc, GSL_EIGEN_SORT_ABS_DESC);      // sort


    // signals = PC' * data;        // <=== MATLAB matlab code

    // Use BLAS to perform the multiplication 
    // C = \alpha * A * B + \beta * C...but  \alpha = 1 and \beta = 0    so we get
    // C = AB.....but   A = A_transform and  B = B_transform   so....
    // C = A_transform * B_transform
    //
    // score = pc' * data;  <== (MATLAB code)
    //
    // The reason why we do data' instead of data is that our "data" is already transformed
    gsl_blas_dgemm(CblasTrans, CblasTrans, 1.0, pc, data, 0.0, score);

    gsl_vector_free(eval);
    gsl_matrix_free(covariance);
}

// Purpose: calculate the average of each column
// Inputs:
//  data - the matrix of interesting data
// Outputs:
//  meanc - the average of each column
void meancols(gsl_vector * meanc, gsl_matrix * data)
{
    sumcols(meanc, data);

    // Now divide by the number of entries
    double dScale = 1.0 / double(data->size1);
    gsl_vector_scale(meanc, dScale);
}


// Purpose: calculate the sum of each column
// Inputs:
//  data - the matrix of interesting data
// Outputs:
//  sumc - the vector that contains the sum of each column
void sumcols(gsl_vector * sumc, gsl_matrix * data)
{
    int rows = data->size1;
    int cols = data->size2;

    for (int j = 0; j < cols; ++j)
    {
        double dVal;
        dVal = 0;
        for (int i = 0; i < rows; ++i)
        {
            dVal += gsl_matrix_get(data, i, j);
        }
        gsl_vector_set(sumc, j, dVal);
    }
}

void PCA::trace(gsl_matrix *mtrace)
{
    double temp;
    int maxx = (mtrace->size1 < 10) ? mtrace->size1 : 10;
    int maxy = (mtrace->size2 < 10) ? mtrace->size2 : 10;


    for (int x = 0; x < maxx; ++x)
    {
        for (int y = 0; y < maxy; ++y)
        {
            temp = gsl_matrix_get(mtrace,x,y);
            printf("%4.4f  \t", temp);
        }
        TRACE("\n");
    }
    TRACE("\n----------------------\n");
}

void PCA::trace(gsl_vector * v)
{
	int I = (v->size < 10) ? v->size : 10;
    for (int i = 0; i < I; ++i)
    {
        //TRACE("%4.4f  \t", gsl_vector_get(v, i));
    }
    TRACE("\n----------------------\n");
}