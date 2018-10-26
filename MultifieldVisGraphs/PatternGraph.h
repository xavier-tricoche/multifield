#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <string.h>
#include <math.h>
#include <time.h>
#include <map>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

using namespace std;

class PatternNode
{ 
public:
	
	int timestep;
	int pattern;
	int size;

	PatternNode()
	{
	}

	PatternNode(int _timestep, int _pattern) 
	{
		timestep = _timestep;
		pattern = _pattern;
		size = 0;
	}
	
	bool operator<(const PatternNode& node) const
	{
		if (timestep < node.timestep)
			return true;

		if ((timestep == node.timestep) && (pattern < node.pattern))
			return true;

		return false;
	};

	bool operator==(const PatternNode& node) const
	{
		if ((timestep == node.timestep) && (pattern == node.pattern))
			return true;
		
		return false;
	}
};




////// (1) undirected graph 
typedef boost::adjacency_list
  <
    //Store all edges as a std::vector
    boost::setS,
    //Store all vertices in a std::vector
    boost::setS,
    //Relations are both ways (in this example)
    //(note: but you can freely change it to boost::directedS)
    boost::undirectedS,
    //All vertices are person names of type std::string
    PatternNode,
    //All edges are weights equal to the encounter frequencies
    boost::property<boost::edge_weight_t,double>,
    //Graph itself has a std::string name
    boost::property<boost::graph_name_t,std::string>
  > PatternGraph_t;
typedef boost::graph_traits<PatternGraph_t>::vertex_descriptor patternvertex_t;
typedef boost::graph_traits<PatternGraph_t>::edge_descriptor patternedge_t;

////// (2) graph to be written in file
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, 
	  PatternNode,
	  boost::property<boost::edge_weight_t, double>,
	  boost::property<boost::graph_name_t,std::string>
> PatternGraph_o;
typedef boost::graph_traits<PatternGraph_o>::vertex_descriptor patternvertex_o;
typedef boost::graph_traits<PatternGraph_o>::edge_descriptor patternedge_o;

struct PatternGraph_writer {

	// my_node_writer() {}
	PatternGraph_writer(PatternGraph_o* g_) : g (g_) {};

	void operator()(std::ostream& out, const patternvertex_o& v) {
		string str;
		int pattern = (*g)[v].pattern;

		for (int i = 0; i < 10; i++)
		{
			if (pattern & 1 == 1)
				str = "1" + str;
			else
				str = "0" + str;

			pattern >>= 1;
		}

		out << "[label=\"" << (*g)[v].timestep << "/" << str << "\",nvert=" << (*g)[v].size << "]";
	};

	void operator()(std::ostream& out, const patternedge_o& e) {
		boost::property_map<PatternGraph_o, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, *g);
		out << "[weight=" << (get(weightmap, e) / 10000.0) << "]";
	};

	// bleah.  I can't get references right...
	// according to http://www.knowledgesearch.org/doc/examples.html
	// it should be a reference here, not the type itself.
	// but g++ either barfs, or the program segfaults.
	PatternGraph_o* g;
};




// main class definition
class PatternGraph
{
public:

	PatternGraph_t graph;
	map<PatternNode, patternvertex_t> nodes;

	PatternGraph(void);
	~PatternGraph(void);

	vector<int>* FindAllPatterns(vtkPolyData* vtkMesh);
	void UpdateGraphWithInput(int ts1, vtkPolyData* mesh1, int ts2, vtkPolyData* mesh2);
	void WriteFile(string filename);
};

