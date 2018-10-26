///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Neuroshare Project                                                         
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
// $Source: /cvsroot/neuroshare/Suite/PowerNAP/PCA.h,v $
//
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PCA_H_INCLUDED
#define PCA_H_INCLUDED

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_vector.h>	
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#define TRACE(x) {printf("%s\n", x); fflush(stdout);}

// Wrapper class to avoid polution of the global namespace
class PCA
{
public:

// Purpose: do the principle component analysis by acting directly on
//  the covariance matrix.
// Inputs:
//  data - the data matrix. It will be "clobbered" by this activity (n x m)
//         m = num of dimensions     and     n = num of trials
// Outputs:
//  pc - the eigen vectors used for the conversion (m x m)
//  score - the data remapped into PC space (n x m)
//          NOTICE, that score data data are transposed in dimension from data
static void FindPCA(gsl_matrix *pc, gsl_matrix *score, gsl_matrix *data);

// Purpose: Transpose the matrix, even if not square
static void non_square_transpose(gsl_matrix * data);

static void trace(gsl_matrix * m);     // Utility function to "output" a matrix
static void trace(gsl_vector * v);     // Utility function to "output" a vector

};


#endif // include guards