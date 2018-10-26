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
#include "ClusterTransitionGraph.h"

double mylog2( double n )
{
    // log(n)/log(2) is log2.
    return log( n ) / log( 2.0 );
}

ClusterTransitionGraph::ClusterTransitionGraph(void)
{
}


ClusterTransitionGraph::~ClusterTransitionGraph(void)
{
}

void ClusterTransitionGraph::CreateGraphNodes(vector<ClusterNode>* clusters)
{
	for (int i = 0; i < clusters->size(); i++)
	{
		clustertransitionvertex_t v = boost::add_vertex(graph);
		graph[v].cluster = clusters->at(i).cluster;
		graph[v].size = clusters->at(i).size;
		graph[v].nvar = clusters->at(i).nvar;
		graph[v].entropy = clusters->at(i).entropy;
		graph[v].mpclusterint = clusters->at(i).mpclusterint;
		graph[v].mpclusterstr = clusters->at(i).mpclusterstr;
		nodes[graph[v]] = v;
	}
}

void ClusterTransitionGraph::UpdateGraphWithInput(int ts1, vtkPolyData* mesh1, int ts2, vtkPolyData* mesh2, vector<ClusterNode>* clusters)
{
	// get the arrays
	vtkIntArray* mesh1clustertransitions = (vtkIntArray*) mesh1->GetPointData()->GetArray("CLUSTBL");
	vtkIntArray* mesh2clustertransitions = (vtkIntArray*) mesh2->GetPointData()->GetArray("CLUSTBL");

	// get links weights map
	boost::property_map<ClusterTransitionGraph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, graph);

	// add links
	vtkIntArray* m12 = (vtkIntArray*) mesh1->GetPointData()->GetArray("FMTBL");
	for (int id1 = 0; id1 < mesh1->GetNumberOfPoints(); id1++)
	{
		int id2 = m12->GetValue(id1);
		if (id2 == -1) continue;

		int comp1 = mesh1clustertransitions->GetValue(id1);
		int comp2 = mesh2clustertransitions->GetValue(id2);

		clustertransitionvertex_t v1 = nodes[clusters->at(comp1)];
		clustertransitionvertex_t v2 = nodes[clusters->at(comp2)];

		bool bRet;
		clustertransitionedge_t e;
		tie(e,bRet) = boost::edge(v1, v2, graph);

		if (bRet == true)
		{
			put(weightmap, e, get(weightmap, e) + 1.0);
		}
		else
		{
			boost::add_edge(v1, v2, 1.0, graph);
		}
	}
}

void ClusterTransitionGraph::WriteFile(string filename)
{
	// initialize Graph
	ClusterTransitionGraph_o ng;

	// copy the vertices
	map<ClusterNode, clustertransitionvertex_o> nodes_o;
	boost::graph_traits<ClusterTransitionGraph_t>::vertex_iterator vi, vend;
	for (boost::tie(vi, vend) = vertices(graph); vi != vend; ++vi) 
	{
		clustertransitionvertex_o v = boost::add_vertex(ng);
		ng[v].cluster = graph[*vi].cluster;
		ng[v].size = graph[*vi].size;
		ng[v].nvar = graph[*vi].nvar;
		ng[v].entropy = graph[*vi].entropy;
		ng[v].mpclusterint = graph[*vi].mpclusterint;
		ng[v].mpclusterstr = graph[*vi].mpclusterstr;
		nodes_o[ng[v]] = v;
	}
	
	// copy the edges
	boost::property_map<ClusterTransitionGraph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, graph);
	boost::graph_traits<ClusterTransitionGraph_t>::edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei) 
	{
		clustertransitionvertex_t gv1 = source(*ei, graph);
		clustertransitionvertex_t gv2 = target(*ei, graph);
		clustertransitionvertex_o v1 = nodes_o[graph[gv1]];
		clustertransitionvertex_o v2 = nodes_o[graph[gv2]];
		boost::add_edge(v1, v2, get(weightmap, *ei), ng);
	}
		
	// write the graph to a file
	std::ofstream output(filename.c_str(), std::ios::out);
	boost::write_graphviz(output, ng, ClusterTransitionGraph_writer(&ng), ClusterTransitionGraph_writer(&ng));
	output.close();
}