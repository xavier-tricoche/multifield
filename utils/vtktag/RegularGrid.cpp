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
#include "RegularGrid.h"

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

RegularGrid::RegularGrid(Nrrd* _ni, int _width, int _height, int _depth, double spcx, double spcy, double spcz)
{
	this->ni = _ni;

	this->width = _width;
	this->height = _height;
	this->depth = _depth;

	xs = (double*) malloc(width * sizeof(double));
	ys = (double*) malloc(height * sizeof(double));
	zs = (double*) malloc(depth * sizeof(double));
	
	if (ISNAN(_ni->axis[0].min))
		_ni->axis[0].min = 0;
	if (ISNAN(_ni->axis[1].min))
		_ni->axis[1].min = 0;
	if (ISNAN(_ni->axis[2].min))
		_ni->axis[2].min = 0;

	for (int i = 0; i < width; i++)
		xs[i] = spcx * i + _ni->axis[0].min;

	for (int j = 0; j < height; j++)
		ys[j] = spcy * j + _ni->axis[1].min;

	for (int k = 0; k < depth; k++)
		zs[k] = spcz * k + _ni->axis[2].min;

	msp = (spcx < spcy)? spcx : spcy;
	msp = (msp < spcz)? msp : spcz;
	msp *= 1;

	gctx = gageContextNew();
	double kparm[3]={1, 1.0, 0.0}; /* parameters of the cubic kernel */
	gagePerVolume *pvl;
	if (!(pvl = gagePerVolumeNew(gctx, ni, gageKindScl))
		|| gagePerVolumeAttach(gctx, pvl)
		|| gageKernelSet(gctx, gageKernel00, nrrdKernelBCCubic, kparm)
		|| gageKernelSet(gctx, gageKernel11, nrrdKernelBCCubicD, kparm)
		|| gageKernelSet(gctx, gageKernel22, nrrdKernelBCCubicDD, kparm)
		|| gageQueryItemOn(gctx, pvl, gageSclValue)
		|| gageQueryItemOn(gctx, pvl, gageSclGradVec)
		|| gageQueryItemOn(gctx, pvl, gageSclNormal)
		|| gageQueryItemOn(gctx, pvl, gageSclHessian)
		|| gageUpdate(gctx)) 
	{
		char *err;
		fprintf(stderr, "ERROR while setting up Gage:\n%s\n", err);
	}

	gctx_sclr = gageAnswerPointer(gctx, pvl, gageSclValue);
	gctx_grad = gageAnswerPointer(gctx, pvl, gageSclGradVec);
	
	errors = 0;
}

RegularGrid::~RegularGrid(void)
{
	free(xs);
	free(ys);
	free(zs);

	gageContextNix(gctx); 
	nrrdNuke(ni);
}

bool RegularGrid::IsBoundary(const int& x, const int& y, const int& z)
{
	if ((x < 1) || (y < 1) || (z < 1) || (x > width - 2) || (y > height - 2) || (z > depth - 2))
		return true;
	else return false;
}

bool RegularGrid::IsValid(const int& x, const int& y, const int& z)
{
	if ((x < 0) || (y < 0) || (z < 0) || (x > width - 1) || (y > height - 1) || (z > depth - 1))
		return false;
	else return true;
}

bool RegularGrid::InDomain(float3 spt)
{
	if ((spt.x < xs[0]) || (spt.x > xs[8 - 1]) 
	 || (spt.y < ys[0]) || (spt.y > ys[8 - 1])
	 || (spt.z < zs[0]) || (spt.z > zs[8 - 1]))
		return false;
	return true;
}

int RegularGrid::Coord2Addr(const int& x, const int& y, const int& z)
{
	return x + width * y + width * height * z;
}

int RegularGrid::Coord2Addr(const int3& coord)
{
	return Coord2Addr(coord.x, coord.y, coord.z);
}

int3 RegularGrid::Addr2Coord(const int& i)
{
	int x = i % width;
	int y = ((i - x) / width) % height;
	int z = (((i - x) / width) - y) / height;

	return make_int3(x, y, z);
}

float3 RegularGrid::Space2Grid(float3 spt)
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

float3 RegularGrid::Grid2Space(float3 gpt)
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

int3 RegularGrid::Space2Regular(float3 spt, int w, int h, int d)
{
	int3 index;
	index.x = cint(w * (spt.x - xs[0]) / (xs[width - 1] - xs[0]));
	index.y = cint(h * (spt.y - ys[0]) / (ys[height - 1] - ys[0]));
	index.z = cint(d * (spt.z - zs[0]) / (zs[depth - 1] - zs[0]));
	return index;
}

int3 RegularGrid::Grid2Regular(float3 gpt, int w, int h, int d)
{
	float3 spt = Grid2Space(gpt);
	return Space2Regular(spt, w, h, d);
}

bool RegularGrid::GetValueAt(float3 pt, double* scalar)
{	
	/*if (!InDomain(pt))
		return false;

	// convert from space to grid coordinates
	float3 gpt = Space2Grid(pt);

	int i = floor(gpt.x);
	int j = floor(gpt.y);
	int k = floor(gpt.z);

	// if out of boundaries stop
	if ((!IsValid(i, j, k)) || (!IsValid(i+1, j, k)) || (!IsValid(i, j+1, k)) || (!IsValid(i+1, j+1, k))
		|| (!IsValid(i, j, k+1)) || (!IsValid(i+1, j, k+1)) || (!IsValid(i, j+1, k+1)) || (!IsValid(i+1, j+1, k+1)))
		return false;*/

	if (gageProbe(gctx, pt.x, pt.y, pt.z)) {
		fprintf(stderr, "probe: trouble:\n(%d) %s\n", gctx->errNum, gctx->errStr);
		return false;
	}

	*scalar = *gctx_sclr;

	return true;
}

double RegularGrid::GetValueAt(double px, double py, double pz)
{
	float3 spt = make_float3((float) px, (float) py, (float)pz);
	float3 gpt = Space2Grid(spt);
	int3 l = make_int3(floor(gpt.x), floor(gpt.y), floor(gpt.z));
	int3 u = make_int3(ceil(gpt.x), ceil(gpt.y), ceil(gpt.z));
	
	if ((!IsValid(l.x, l.y, l.z)) || (!IsValid(u.x, u.y, u.z)))
	{
		//printf("Out of boundary.\n");
		errors++;
		return 0.0;
	}
	
	/*double x = gpt.x - l.x;
	double y = gpt.y - l.y;
	double z = gpt.z - l.z;
	
	float* data = (float*) ni->data;
	float* ptr;

	ptr = &(data[Coord2Addr(l.x, l.y, l.z)]);
	float V000 = ptr[0];
	ptr = &(data[Coord2Addr(l.x, l.y, u.z)]);
	float V001 = ptr[0];
	ptr = &(data[Coord2Addr(l.x, u.y, l.z)]);
	float V010 = ptr[0];
	ptr = &(data[Coord2Addr(l.x, u.y, u.z)]);
	float V011 = ptr[0];
	ptr = &(data[Coord2Addr(u.x, l.y, l.z)]);
	float V100 = ptr[0];
	ptr = &(data[Coord2Addr(u.x, l.y, u.z)]);
	float V101 = ptr[0];
	ptr = &(data[Coord2Addr(u.x, u.y, l.z)]);
	float V110 = ptr[0];
	ptr = &(data[Coord2Addr(u.x, u.y, u.z)]);
	float V111 = ptr[0];
	
	float V =   V000 * (1 - x) * (1 - y) * (1 - z) +
				V100 * x * (1 - y) * (1 - z) +
				V010 * (1 - x) * y * (1 - z) +
				V001 * (1 - x) * (1 - y) * z +
				V101 * x * (1 - y) * z +
				V011 * (1 - x) * y * z +
				V110 * x * y * (1 - z) +
				V111 * x * y * z; 
	
	return V;*/
	
	if (gageProbe(gctx, gpt.x, gpt.y, gpt.z)) {
		fprintf(stderr, "probe: trouble:\n(%d) %s\n", gctx->errNum, gctx->errStr);
		return 0.0;
	}

	return *gctx_sclr;
}
