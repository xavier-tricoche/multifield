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

double mylog2( double n );

class ClusterNode {
public:
	string cluster;
	int size;
	int nvar;
	double entropy;

	string mpclusterstr;
	int mpclusterint;
	
	ClusterNode()
	{
		cluster = string("");
		size = 0;
		nvar = 0;
		entropy = 0.0;
	}

	ClusterNode(int _nvar)
	{
		size = 0;
		nvar = _nvar;
		entropy = 0;

		for (int i = 0; i < 2 * nvar; i++)
			cluster += "x";
	}
	
	ClusterNode(ClusterNode c, int bit, int val)
	{
		cluster = string(c.cluster);
		nvar = c.nvar;
		if (val == 0)
			cluster[bit] = '0';
		else
			cluster[bit] = '1';
		entropy = 0;

		mpclusterint = c.mpclusterint;
		mpclusterstr = c.mpclusterstr;
	}
	
	bool TagMatch(int tag)
	{
		bool match = true;
		for (int i = 0; i < 2 * nvar; i++)
		{
			int bit = ((tag >> i) & 1);
			
			if (((bit == 1) && (cluster[i] == '0')) || 
				((bit == 0) && (cluster[i] == '1')))
			{
				match = false;
				return match;
			}
		}
		return match;
	}
	
	void ComputeSize(int* tagcnt)
	{
		size = 0;
		
		int ntags = std::pow(2.0, 2 * nvar);
		for (int i = 0; i < ntags; i++)
		{
			if (TagMatch(i))
				size += tagcnt[i];
		}
	}

	void ComputeEntropy(int* tagcnt)
	{
		int ntags = std::pow(2.0, 2 * nvar);
		for (int i = 0; i < ntags; i++)
		{
			if (!TagMatch(i))
				continue;

			double p = ((double) tagcnt[i]) / size;

			if (p != 0.0)
				entropy += -p * mylog2(p);
		}
	}

	void ComputeMostProbableString(int* tagcnt)
	{
		string mps;

		int ntags = std::pow(2.0, 2 * nvar);

		// loop on bits
		for (int bit = 0; bit < 2 * nvar; bit++)
		{
			if (cluster[bit] != 'x')
			{
				mps += cluster[bit];
				continue;
			}
				
			int sum0 = 0;
			int sum1 = 0;
			
			// loop on all matching tags
			for (int tag = 0; tag < ntags; tag++)
			{
				if (!TagMatch(tag)) continue;
				
				if (((tag >> bit) & 1) == 0)
					sum0 += tagcnt[tag];
				else
					sum1 += tagcnt[tag];
			}
			
			if (sum0 > sum1)
				mps += "'0";
			else
				mps += "'1";

			cout << bit << ": " << sum0 << " " << sum1 << "\n";
		}

		mpclusterstr = mps;

		cout << " -> " << string(mpclusterstr.rbegin(), mpclusterstr.rend()) << "\n";
	}

	void ComputeMostProbableInt(int* tagcnt)
	{
		int mps = 0;

		int ntags = std::pow(2.0, 2 * nvar);

		// loop on bits
		for (int bit = 0; bit < 2 * nvar; bit++)
		{
			if (cluster[bit] != 'x')
			{
				if (cluster[bit] == '1')
					mps |= (1 << bit);
				continue;
			}
				
			int sum0 = 0;
			int sum1 = 0;
			
			// loop on all matching tags
			for (int tag = 0; tag < ntags; tag++)
			{
				if (!TagMatch(tag)) continue;
				
				if (((tag >> bit) & 1) == 0)
					sum0 += tagcnt[tag];
				else
					sum1 += tagcnt[tag];
			}
			
			if (sum0 < sum1)
				mps |= (1 << bit);
		}

		mpclusterint = mps;
	}
	
	int FindBestBitForSplit(int* tagcnt)
	{
		double minentropy = std::numeric_limits<double>::max();
		int bestbit = -1;
		int ntags = std::pow(2.0, 2 * nvar);
		
		// loop on bits
		for (int bit = 0; bit < 2 * nvar; bit++)
		{
			if (cluster[bit] != 'x')
				continue;
				
			// S1
			ClusterNode c1(*this, bit, 0);
			c1.ComputeSize(tagcnt);
			c1.ComputeEntropy(tagcnt);
		
			// S2
			ClusterNode c2(*this, bit, 1);
			c2.ComputeSize(tagcnt);
			c2.ComputeEntropy(tagcnt);

			// H(bit;S)
			double has = ((double) c1.size / size) * c1.entropy + ((double) c2.size / size) * c2.entropy;

			if (has < minentropy)
			{
				bestbit = bit;
				minentropy = has;
			}
		}
		
		return bestbit;
	}
	
	void Print()
	{
		cout << string(cluster.rbegin(), cluster.rend()) << " " << string(mpclusterstr.rbegin(), mpclusterstr.rend()) << " " << size << " " << entropy << "\n";
	}

	bool operator<(const ClusterNode& e2) const
	{
		if (entropy == e2.entropy)
		{
			return (cluster.compare(e2.cluster) < 0);
		}

		if (entropy < e2.entropy)
			return true;
		else
			return false;
	}

	bool operator==(const ClusterNode& node) const
	{
		return (cluster.compare(node.cluster) == 0);
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
    boost::directedS,
    //All vertices are person names of type std::string
    ClusterNode,
    //All edges are weights equal to the encounter frequencies
    boost::property<boost::edge_weight_t,double>,
    //Graph itself has a std::string name
    boost::property<boost::graph_name_t,std::string>
  > ClusterTransitionGraph_t;
typedef boost::graph_traits<ClusterTransitionGraph_t>::vertex_descriptor clustertransitionvertex_t;
typedef boost::graph_traits<ClusterTransitionGraph_t>::edge_descriptor clustertransitionedge_t;

////// (2) graph to be written in file
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, 
	  ClusterNode,
	  boost::property<boost::edge_weight_t, double>,
	  boost::property<boost::graph_name_t,std::string>
> ClusterTransitionGraph_o;
typedef boost::graph_traits<ClusterTransitionGraph_o>::vertex_descriptor clustertransitionvertex_o;
typedef boost::graph_traits<ClusterTransitionGraph_o>::edge_descriptor clustertransitionedge_o;

struct ClusterTransitionGraph_writer {

	// my_node_writer() {}
	ClusterTransitionGraph_writer(ClusterTransitionGraph_o* g_) : g (g_) {};

	void operator()(std::ostream& out, const clustertransitionvertex_o& v) {
		out << "[label=\"" << string((*g)[v].mpclusterstr.rbegin(), (*g)[v].mpclusterstr.rend()) << "\",nvert=" << (*g)[v].size << ",entropy=" << (*g)[v].entropy << "]";
	};

	void operator()(std::ostream& out, const clustertransitionedge_o& e) {
		boost::property_map<ClusterTransitionGraph_o, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, *g);
		out << "[weight=" << (get(weightmap, e) / 10000.0) << "]";
	};

	// bleah.  I can't get references right...
	// according to http://www.knowledgesearch.org/doc/examples.html
	// it should be a reference here, not the type itself.
	// but g++ either barfs, or the program segfaults.
	ClusterTransitionGraph_o* g;
};


// main class definition
class ClusterTransitionGraph
{
public:

	ClusterTransitionGraph_t graph;
	map<ClusterNode, clustertransitionvertex_t> nodes;

	ClusterTransitionGraph(void);
	~ClusterTransitionGraph(void);

	void CreateGraphNodes(vector<ClusterNode>* clusters);
	void UpdateGraphWithInput(int ts1, vtkPolyData* mesh1, int ts2, vtkPolyData* mesh2, vector<ClusterNode>* clusters);
	void WriteFile(string filename);
};

