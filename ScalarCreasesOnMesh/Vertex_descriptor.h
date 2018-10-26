#ifndef VERTEXDESCRIPTION
#define VERTEXDESCRIPTION

#define CGAL_NO_AUTOLINK_GMP 1
#define CGAL_NO_AUTOLINK_MPFR 1

#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <queue>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <omp.h>
#include <vector_types.h>
#include <vector_functions.h>
#include <cutil_inline.h>
#include <cutil_math.h>
#include <math/fixed_vector.hpp>
#include <math/bezier_spline.hpp>
#include <math/dopri5.hpp>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <teem/nrrd.h>
#include <teem/ten.h>
#include <teem/seek.h>
#include <teem/air.h>
#include <teem/limn.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Delaunay_Triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/nearest_neighbor_delaunay_2.h>

class Vertex_descriptor;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_descriptor*, K> Vb;
typedef CGAL::Triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds> CGAL_Triangulation;
typedef CGAL_Triangulation::Vertex_handle CGAL_Vertex_handle;
typedef CGAL_Triangulation::Face_handle CGAL_Face_handle;
typedef CGAL_Triangulation::Point  CGAL_Point;
typedef CGAL_Triangulation::Edge  CGAL_Edge;
typedef CGAL_Triangulation::Finite_faces_iterator CGAL_Finite_faces_iterator;
typedef CGAL_Triangulation::Finite_vertices_iterator CGAL_Finite_vertices_iterator;
typedef CGAL_Triangulation::All_edges_iterator CGAL_All_edges_iterator;
typedef CGAL_Triangulation::All_faces_iterator CGAL_All_faces_iterator;
typedef CGAL_Triangulation::Edge_iterator CGAL_Edge_iterator;
typedef CGAL_Triangulation::Segment CGAL_Segment;
typedef CGAL_Triangulation::Vertex_iterator CGAL_Vertex_iterator;

enum Vertex_type {CENTER, LINK, TIP};


class Vertex_descriptor
{
public:
	static int contour_counter;

	int contourid;
	Vertex_type type;
	CGAL_Vertex_handle vh;
	Vertex_descriptor* frwrd;
	Vertex_descriptor* bkwrd;
	float2 coordinates;
	double value;

	Vertex_descriptor(Vertex_type _type, Vertex_descriptor* _frwrd, Vertex_descriptor* _bkwrd, float2 _coordinates, double _value)
	{
		type = _type;
		frwrd = _frwrd;
		bkwrd = _bkwrd;
		coordinates = _coordinates;
		value = _value;
		vh = NULL;
	}
};

int Vertex_descriptor::contour_counter = 0;

struct TargetPoint
{
	float2 origin;
	float2 target;
	double oldprob;
	double newprob;
};

inline bool operator==(const float2& a1, const float2& a2)
{
    return ((a1.x == a2.x) && (a1.y == a2.y));
}


#endif