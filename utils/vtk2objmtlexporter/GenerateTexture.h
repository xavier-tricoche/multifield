#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector_types.h>
#include <vector_functions.h>
#include <cutil_inline.h>
#include <cutil_math.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPNGWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkActor.h>
#include <vtkOBJExporter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>

#include "EasyBMP.h"

using namespace std;

typedef unsigned char byte;
typedef struct
{
	byte red, green, blue, alpha;
} RGB_t;

int cint(double x);
double dist_Point_to_Line(float3 P, float3 LP0, float3 LP1);
int write_truecolor_tga(string filename, RGB_t* data, unsigned width, unsigned height );
void WriteBMPImage(string filename, RGB_t* data, int w, int h);
void GenerateColorTexture(vtkPolyData* vtkMesh, string texprefix, FILE* objfile);


class Triangle
{
	vtkPolyData* mesh;
	vtkIdType pts[3];
	float3 ptsc[3];

	double height;
	double baselength;
	
	int texture;
	float2 uv[3];
	
	
public:
	
	Triangle(vtkPolyData* _mesh, vtkIdType* _pts)
	{
		mesh = _mesh;
		pts[0] = _pts[0];
		pts[1] = _pts[1];
		pts[2] = _pts[2];
		
		// get points coordinates
		double c[3];
		_mesh->GetPoint(pts[0], c);
		ptsc[0] = make_float3(c[0], c[1], c[2]);
		_mesh->GetPoint(pts[1], c);
		ptsc[1] = make_float3(c[0], c[1], c[2]);
		_mesh->GetPoint(pts[2], c);
		ptsc[2] = make_float3(c[0], c[1], c[2]);
		
		// compute the length of the triangle
		double l1 = length(ptsc[0] - ptsc[1]);
		double l2 = length(ptsc[1] - ptsc[2]);
		double l3 = length(ptsc[0] - ptsc[2]);
		double ml = max(l1, max(l2, l3));
		if (l2 == ml)
		{
			float3 tmp = ptsc[0];
			ptsc[0] = ptsc[2];
			ptsc[2] = tmp;
			
			vtkIdType tmpi = pts[0];
			pts[0] = pts[2];
			pts[2] = tmpi;
		}
		if (l3 == ml)
		{
			float3 tmp = ptsc[1];
			ptsc[1] = ptsc[2];
			ptsc[2] = tmp;
			
			vtkIdType tmpi = pts[1];
			pts[1] = pts[2];
			pts[2] = tmpi;
		}
		baselength = ml;
		
		// compute the height of the triangle
		height = dist_Point_to_Line(ptsc[2], ptsc[0], ptsc[1]);
	}
	
	double GetHeight()
	{
		return height;
	}
	
	double GetBaseLength()
	{
		return baselength;
	}
	
	int GetPixelHeight(double scale)
	{
		return ceil(height * scale) + 3;
	}
	
	int GetPixelBaseLength(double scale)
	{
		return ceil(baselength * scale) + 3;
	}
	
	void RasterizeColors(int id, RGB_t* data, int xsp, int ysp, double scale, 
		FILE* objfile, int locx, int locy, int& vtcount,
		int w, int h)
	{
		texture = id;
		
		// get points colors
		vtkUnsignedCharArray* colors = (vtkUnsignedCharArray*) mesh->GetPointData()->GetScalars();
		float4 c[3];
		double color[4];
		colors->GetTuple(pts[0], color);
		c[0] = make_float4(color[0], color[1], color[2], color[3]);
		colors->GetTuple(pts[1], color);
		c[1] = make_float4(color[0], color[1], color[2], color[3]);
		colors->GetTuple(pts[2], color);
		c[2] = make_float4(color[0], color[1], color[2], color[3]);
		
		// compute the uv coordinates of the three points
		uv[0] = make_float2(1.0, 1.0);
		uv[1] = make_float2(1.0 + scale * baselength, 1.0);
		uv[2] = make_float2(1.0 + scale * dot(ptsc[2] - ptsc[0], normalize(ptsc[1] - ptsc[0])), scale * height + 1.0);
		
		// write to file
		fprintf(objfile, "vt %f %f\n", float(float(locx)+uv[0].x) / float(w - 1), float(float(locy)+uv[0].y) / float(h - 1));
		fprintf(objfile, "vt %f %f\n", float(float(locx)+uv[1].x) / float(w - 1), float(float(locy)+uv[1].y) / float(h - 1));
		fprintf(objfile, "vt %f %f\n", float(float(locx)+uv[2].x) / float(w - 1), float(float(locy)+uv[2].y) / float(h - 1));
		fprintf(objfile, "f %d/%d %d/%d %d/%d\n", pts[0] + 1, vtcount + 1, pts[1] + 1, vtcount + 2, pts[2] + 1, vtcount + 3);
		vtcount += 3;
		
		// no loop adding the colors
		for (int y = 0; y < GetPixelHeight(scale); y++)
		{
			for (int x = 0; x < GetPixelBaseLength(scale); x++)
			{
				float2 uvc = make_float2(x, y);

				// get barycentric coordinates
				double l1 = ((uv[1].y - uv[2].y) * (uvc.x - uv[2].x) + (uv[2].x - uv[1].x) * (uvc.y - uv[2].y)) / 
							((uv[1].y - uv[2].y) * (uv[0].x - uv[2].x) + (uv[2].x - uv[1].x) * (uv[0].y - uv[2].y));
				double l2 = ((uv[2].y - uv[0].y) * (uvc.x - uv[2].x) + (uv[0].x - uv[2].x) * (uvc.y - uv[2].y)) / 
							((uv[1].y - uv[2].y) * (uv[0].x - uv[2].x) + (uv[2].x - uv[1].x) * (uv[0].y - uv[2].y));
				double l3 = 1 - l1 - l2;
				if (l1 < 0.0)
					l1 = 0.0;
				if (l2 < 0.0)
					l2 = 0.0;
				if (l3 < 0.0)
					l3 = 0.0;
				double sum = l1 + l2 + l3;
				l1 /= sum;
				l2 /= sum;
				l3 /= sum;
				
				float4 cc = make_float4(255.0);
				cc = l1 * c[0] + l2 * c[1] + l3 * c[2];
				data[x * xsp + y * ysp].red = cint(cc.x);
				data[x * xsp + y * ysp].green = cint(cc.y);
				data[x * xsp + y * ysp].blue = cint(cc.z);
				data[x * xsp + y * ysp].alpha = cint(cc.w);
			}
		}
	}
};

class QuadraticTexture
{
public:

	int id;
	string filename;
	RGB_t* data;
	int width;
	int height;

	int locx;
	int locy;
	
	int hs;
	
	QuadraticTexture(int _id, string prefix, int _width, int _height)
	{
		id = _id;
		
		char buffer [33];
		sprintf(buffer,"%d",_id);
  
		filename = string(prefix) + string("_") + string(buffer);
		width = _width;
		height = _height;
		data = (RGB_t*) malloc(width * height * sizeof(RGB_t));
		memset(data, 255, width * height * sizeof(RGB_t));
		
		locx = 0;
		locy = 0;
		
		hs = 0;
	}
	
	void InsertTriangleInNewRow(Triangle& t, double scale, FILE* objfile, int& vtcount)
	{
		//printf("new: %d %d\n", locx, locy);//getchar();
		
		locx = 0;
		locy += hs;
		
		// rasterize to texture
		t.RasterizeColors(id, &(data[locx + width * locy]), 1, width, scale, 
							objfile, locx, locy, vtcount,
							width, height);
		locx += t.GetPixelBaseLength(scale);
		
		// save height information for next row
		hs = t.GetPixelHeight(scale);
	}
	
	void InsertTriangleInSameRow(Triangle& t, double scale, FILE* objfile, int& vtcount)
	{
		//printf("cont: %d %d\n", locx, locy);//getchar();
		
		// rasterize to texture
		t.RasterizeColors(id, &(data[locx + width * locy]), 1, width, scale,
							objfile, locx, locy, vtcount,
							width, height);
		locx += t.GetPixelBaseLength(scale);
	}
	
	void Destruct()
	{
		free(data);
	}
};

class CompareHeight
{
public:
    // Compare two Foo structs.
    bool operator()(Triangle& x, Triangle& y) const
    {
        return x.GetHeight() < y.GetHeight();
    }
};
