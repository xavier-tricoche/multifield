#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <string>
#include <iostream>
#include <vector_types.h>
#include <vector_functions.h>
#include <cutil_inline.h>
#include <cutil_math.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <vector>
#include <list>
#include <teem/nrrd.h>
#include <teem/ten.h>
#include <teem/seek.h>
#include "mymath.h"

using namespace std;

typedef unsigned int uint;
typedef unsigned char uchar;
typedef struct { float4 m[3]; } double3x4;
typedef struct { double3 m[3]; } double3x3;


int cint(double x);
void normalize(double* vec);


class RegularGrid
{
public:
	RegularGrid(Nrrd* _ni, int _width, int _height, int _depth, double spcx, double spcy, double spcz);
	~RegularGrid(void);


	////// Properties
	float* data; 
	int width;
	int height;
	int depth;
	double* xs;
	double* ys;
	double* zs;
	double msp;
	double xspc;
	double yspc;
	double zspc;

	////// Methods
	bool IsBoundary(const int& x, const int& y, const int& z);
	bool IsValid(const int& x, const int& y, const int& z);
	bool InDomain(const float3& spt);
	int Coord2Addr(const int& x, const int& y, const int& z);
	int Coord2Addr(const int3& coord);
	int3 Addr2Coord(const int& i);
	float3 Space2Grid(const float3& spt);
	float3 Grid2Space(const float3& gpt);
	int3 Space2Regular(float3 spt, int w, int h, int d);
	int3 Grid2Regular(float3 gpt, int w, int h, int d);


	bool GetValueAt(float3 v, double* scalar);
	double GetValueAt(double x, double y, double z);
};

