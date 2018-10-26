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
#include "FlowGrid.h"

int cint(double x){

	double fractpart, intpart;
	fractpart = modf (x , &intpart);

	if (fractpart>=.5)
		return x>=0?ceil(x):floor(x);
	else
		return x<0?ceil(x):floor(x);
}

void normalize(double* vec)
{
	double a = sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
	if (a == 0.0) return;
	vec[0] /= a;
	vec[1] /= a;
}

inline nvis::vec2 array2vec(const double* a, int i) {
	return nvis::vec2(a[2*i], a[2*i+1]);
}

inline nvis::mat2 tensor2mat(const double* a) {
	nvis::mat2 A;
	A(0,0) = a[0];
	A(0,1) = A(1,0) = a[1];
	A(1,1) = a[3];
	return A;
}

inline nvis::mat2 tensorgrad2mat(const double *a, int i) {
	nvis::mat2 A;
	A(0,0) = a[i*4];
	A(0,1) = A(1,0) = a[i*4+1];
	A(1,1) = a[i*4+3];
	return A;
}

nvis::mat2 PseudoInverse(nvis::mat2 m) {

	int n_fil = 2;
	int n_col = 2;

	unsigned i = 0;
	unsigned j = 0;
	gsl_matrix * gA = gsl_matrix_alloc (n_fil, n_col);
	for (i = 0; i < n_fil; i++)
		for (j = 0; j < n_col; j++)
	   		gsl_matrix_set (gA, i, j, m(i,j));


	gsl_matrix * gA_t = gsl_matrix_alloc (n_col, n_fil);
	gsl_matrix_transpose_memcpy (gA_t, gA);					// Computing the transpose of gA
		
	gsl_matrix * U = gsl_matrix_alloc (n_col, n_fil);
	gsl_matrix * V= gsl_matrix_alloc (n_fil, n_fil);
	gsl_vector * S = gsl_vector_alloc (n_fil);


	// Computing the SVD of the transpose of A
	// The matrix 'gA_t' will contain 'U' after the function is called
	gsl_vector * work = gsl_vector_alloc (n_fil);
	gsl_linalg_SV_decomp (gA_t, V, S, work);
	gsl_vector_free(work);
		
	gsl_matrix_memcpy (U, gA_t);
	
				
	//Inverting S//
	//----------------------------------------------------------
	// Matrix 'S' is diagonal, so it is contained in a vector.
	// We operate to convert the vector 'S' into the matrix 'Sp'.
	//Then we invert 'Sp' to 'Spu'
	//----------------------------------------------------------
	gsl_matrix * Sp = gsl_matrix_alloc (n_fil, n_fil);
	gsl_matrix_set_zero (Sp);
	for (i = 0; i < n_fil; i++)
		gsl_matrix_set (Sp, i, i, gsl_vector_get(S, i));	// Vector 'S' to matrix 'Sp'
	
	gsl_permutation * p = gsl_permutation_alloc (n_fil);
	int signum;
	gsl_linalg_LU_decomp (Sp, p, &signum);				// Computing the LU decomposition

	// Compute the inverse like in the MATLAB script

	gsl_matrix * SI = gsl_matrix_calloc (n_fil, n_fil);

	for (i = 0; i < n_fil; i++) {
	  //std::cout << "S [" << i << "] = " << gsl_vector_get (S, i) << std::endl;

	  if (gsl_vector_get (S, i) > 0.0000000001)
	    gsl_matrix_set (SI, i, i, 1.0 / gsl_vector_get (S, i));
	}		
	
	gsl_matrix * VT = gsl_matrix_alloc (n_fil, n_fil);
	gsl_matrix_transpose_memcpy (VT, V);					// Tranpose of V
		
		
	//THE PSEUDOINVERSE//
	//----------------------------------------------------------
	//Computation of the pseudoinverse of trans(A) as pinv(A) = U·inv(S).trans(V)   with trans(A) = U.S.trans(V)
	//----------------------------------------------------------
	gsl_matrix * SIpVT = gsl_matrix_alloc (n_fil, n_fil);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,				// Calculating  inv(S).trans(V)
                	1.0, SI, VT,
                	0.0, SIpVT);

			
	gsl_matrix * pinv = gsl_matrix_alloc (n_col, n_fil);	// Calculating  U·inv(S).trans(V)
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                	1.0, U, SIpVT,
                	0.0, pinv);

	//end THE PSEUDOINVERSE//
	
	nvis::mat2 ret;
   	for (i = 0; i < n_col; i++)  
    	for (j = 0; j < n_fil; j++)
    		ret(i,j) = gsl_matrix_get (pinv, i, j);

	gsl_matrix_free(VT);
	gsl_matrix_free(SI);
	gsl_matrix_free(SIpVT);
	gsl_matrix_free(gA_t);
	gsl_matrix_free(U);
	gsl_matrix_free(gA);
	gsl_matrix_free(V);
	gsl_vector_free(S);

	return ret;
}

void writeNrrd2(void* data, const string& filename, int data_type, const vector< size_t >& dims, const vector< double >& spacing)
{
	Nrrd *nout = nrrdNew();

	if (nrrdWrap_nva(nout, data, data_type, dims.size(), &dims[0])) {
		cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
		exit(-1);
	}
	if (spacing.size() == dims.size()) {
		nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, (const void*)&spacing[0]);
	}
	if (nrrdSave(filename.c_str(), nout, NULL)) {
		cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
		exit(-1);
	}
}

FlowGrid::FlowGrid(Nrrd* _ni, double _time)
{
	this->ni = _ni;
	this->time = _time;

	this->width = ni->axis[1].size;
	this->height = ni->axis[2].size;
	this->depth = ni->axis[3].size;

	xs = (double*) malloc(width * sizeof(double));
	ys = (double*) malloc(height * sizeof(double));
	zs = (double*) malloc(depth * sizeof(double));
	
	if (ISNAN(_ni->axis[1].min))
		_ni->axis[1].min = 0;
	if (ISNAN(_ni->axis[2].min))
		_ni->axis[2].min = 0;
	if (ISNAN(_ni->axis[3].min))
		_ni->axis[3].min = 0;
		
	double spcx = ni->axis[1].spacing;
	double spcy = ni->axis[2].spacing;
	double spcz = ni->axis[3].spacing;

	for (int i = 0; i < width; i++)
		xs[i] = spcx * i + _ni->axis[1].min;

	for (int j = 0; j < height; j++)
		ys[j] = spcy * j + _ni->axis[2].min;

	for (int k = 0; k < depth; k++)
		zs[k] = spcz * k + _ni->axis[3].min;

	msp = (spcx < spcy)? spcx : spcy;
	msp = (msp < spcz)? msp : spcz;
	msp *= 1;
}

FlowGrid::~FlowGrid(void)
{
	free(xs);
	free(ys);
	free(zs);

	nrrdNuke(ni);
}

bool FlowGrid::IsBoundary(const int& x, const int& y, const int& z)
{
	if ((x < 1) || (y < 1) || (z < 1) || (x > width - 2) || (y > height - 2) || (z > depth - 2))
		return true;
	else return false;
}

bool FlowGrid::IsValid(const int& x, const int& y, const int& z)
{
	if ((x < 0) || (y < 0) || (z < 0) || (x > width - 1) || (y > height - 1) || (z > depth - 1))
		return false;
	else return true;
}

bool FlowGrid::InDomain(float3 spt)
{
	if ((spt.x < xs[0]) || (spt.x > xs[8 - 1]) 
	 || (spt.y < ys[0]) || (spt.y > ys[8 - 1])
	 || (spt.z < zs[0]) || (spt.z > zs[8 - 1]))
		return false;
	return true;
}

int FlowGrid::Coord2Addr(const int& x, const int& y, const int& z)
{
	return x + width * y + width * height * z;
}

int FlowGrid::Coord2Addr(const int3& coord)
{
	return Coord2Addr(coord.x, coord.y, coord.z);
}

int3 FlowGrid::Addr2Coord(const int& i)
{
	int x = i % width;
	int y = ((i - x) / width) % height;
	int z = (((i - x) / width) - y) / height;

	return make_int3(x, y, z);
}

float3 FlowGrid::Space2Grid(float3 spt)
{
	float3 gpt;

	for (int i = 0; i < width - 1; i++)
	{
		if ((spt.x >= xs[i]) && (spt.x < xs[i+1]))
		{
			gpt.x = i + (spt.x - xs[i]) / (xs[i+1] - xs[i]);
			break;
		}
	}

	for (int i = 0; i < height - 1; i++)
	{
		if ((spt.y >= ys[i]) && (spt.y < ys[i+1]))
		{
			gpt.y = i + (spt.y - ys[i]) / (ys[i+1] - ys[i]);
			break;
		}
	}

	for (int i = 0; i < depth - 1; i++)
	{
		if ((spt.z >= zs[i]) && (spt.z < zs[i+1]))
		{
			gpt.z = i + (spt.z - zs[i]) / (zs[i+1] - zs[i]);
			break;
		}
	}

	return gpt;
}

float3 FlowGrid::Grid2Space(float3 gpt)
{
	float3 spt;

	for (int i = 0; i < width - 1; i++)
	{
		if ((gpt.x >= i) && (gpt.x <= i+1))
		{
			spt.x = xs[i] + (gpt.x - i) * (xs[i+1] - xs[i]);
			break;
		}
	}

	for (int i = 0; i < height - 1; i++)
	{
		if ((gpt.y >= i) && (gpt.y <= i+1))
		{
			spt.y = ys[i] + (gpt.y - i) * (ys[i+1] - ys[i]);
			break;
		}
	}

	for (int i = 0; i < depth - 1; i++)
	{
		if ((gpt.z >= i) && (gpt.z <= i+1))
		{
			spt.z = zs[i] + (gpt.z - i) * (zs[i+1] - zs[i]);
			break;
		}
	}

	return spt;
}

int3 FlowGrid::Space2Regular(float3 spt, int w, int h, int d)
{
	int3 index;
	index.x = cint(w * (spt.x - xs[0]) / (xs[width - 1] - xs[0]));
	index.y = cint(h * (spt.y - ys[0]) / (ys[height - 1] - ys[0]));
	index.z = cint(d * (spt.z - zs[0]) / (zs[depth - 1] - zs[0]));
	return index;
}

int3 FlowGrid::Grid2Regular(float3 gpt, int w, int h, int d)
{
	float3 spt = Grid2Space(gpt);
	return Space2Regular(spt, w, h, d);
}

float3 FlowGrid::GetFlowAt(float3 spt, bool& failed)
{	
	float3 gpt = Space2Grid(spt);
	int3 l = make_int3(floor(gpt.x), floor(gpt.y), floor(gpt.z));
	int3 u = make_int3(ceil(gpt.x), ceil(gpt.y), ceil(gpt.z));
	
	if ((!IsValid(l.x, l.y, l.z)) || (!IsValid(u.x, u.y, u.z)))
	{
		failed = true;
		return make_float3(0.0);
	}
	
	double x = gpt.x - l.x;
	double y = gpt.y - l.y;
	double z = gpt.z - l.z;
	
	float* data = (float*) ni->data;
	float* ptr;

	ptr = &(data[3 * Coord2Addr(l.x, l.y, l.z)]);
	float3 V000 = make_float3(ptr[0], ptr[1], ptr[2]);
	ptr = &(data[3 * Coord2Addr(l.x, l.y, u.z)]);
	float3 V001 = make_float3(ptr[0], ptr[1], ptr[2]);
	ptr = &(data[3 * Coord2Addr(l.x, u.y, l.z)]);
	float3 V010 = make_float3(ptr[0], ptr[1], ptr[2]);
	ptr = &(data[3 * Coord2Addr(l.x, u.y, u.z)]);
	float3 V011 = make_float3(ptr[0], ptr[1], ptr[2]);
	ptr = &(data[3 * Coord2Addr(u.x, l.y, l.z)]);
	float3 V100 = make_float3(ptr[0], ptr[1], ptr[2]);
	ptr = &(data[3 * Coord2Addr(u.x, l.y, u.z)]);
	float3 V101 = make_float3(ptr[0], ptr[1], ptr[2]);
	ptr = &(data[3 * Coord2Addr(u.x, u.y, l.z)]);
	float3 V110 = make_float3(ptr[0], ptr[1], ptr[2]);
	ptr = &(data[3 * Coord2Addr(u.x, u.y, u.z)]);
	float3 V111 = make_float3(ptr[0], ptr[1], ptr[2]);
	
	float3 V =  V000 * (1 - x) * (1 - y) * (1 - z) +
				V100 * x * (1 - y) * (1 - z) +
				V010 * (1 - x) * y * (1 - z) +
				V001 * (1 - x) * (1 - y) * z +
				V101 * x * (1 - y) * z +
				V011 * (1 - x) * y * z +
				V110 * x * y * (1 - z) +
				V111 * x * y * z; 
	
	return V;
}

float3 TimeInterpolate(double t0, double t1, float3 V0, float3 V1, double t)
{
	double frac = (t1 - t) / (t1 - t0);
	return V0 * frac + V1 * (1.0 - frac);
}
