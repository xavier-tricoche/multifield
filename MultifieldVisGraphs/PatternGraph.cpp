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
#include "PatternGraph.h"


PatternGraph::PatternGraph(void)
{
}


PatternGraph::~PatternGraph(void)
{
}

vector<int>* PatternGraph::FindAllPatterns(vtkPolyData* vtkMesh)
{
	vtkIntArray* tagarr = (vtkIntArray*) vtkMesh->GetPointData()->GetArray("TAGTBL");
	
	// find all tags
	set<int> tags;
	for (vtkIdType i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
		tags.insert(tagarr->GetValue(i));
	printf("Number of distinct tags is %d\n", tags.size());	
		
	// initialize counts to zero
	vector<int> counts;
	for (set<int>::iterator it = tags.begin(); it != tags.end(); it++)
		counts.push_back(0);
		
	// now assign colors
	for (vtkIdType i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
	{
		int tag = tagarr->GetValue(i);
		set<int>::iterator it = tags.find(tag);
		int cc = std::distance(tags.begin(), it);
		
		counts[cc] ++;
	}
	
	// sort counts
	vector<int> countscp = counts;
	sort(countscp.begin(), countscp.end());

	// find the sizelimit for the tags to include
	int sizelimit = 0;
	int sum = 0;
	int i = countscp.size() - 1;
	while ((i >= 0) && (sum < 0.25 * vtkMesh->GetNumberOfPoints()))
	{
		sum += countscp[i];
		sizelimit = countscp[i];
		i--;
	}

	// find tags with a certain size limit
	vector<int>* allpatterns = new vector<int>();
	for (set<int>::iterator it = tags.begin(); it != tags.end(); it++)
	{
		int cc = std::distance(tags.begin(), it);
		if (counts[cc] >= sizelimit)
		{
			allpatterns->push_back(*it);
		}
	}

	// assign patterns to node
	vtkSmartPointer<vtkIntArray> newScalars = vtkSmartPointer<vtkIntArray>::New();
	newScalars->SetName("PATTBL");
	int max_dist = 0;
	for (vtkIdType i = 0; i < vtkMesh->GetNumberOfPoints(); i++)
	{
		int tag = tagarr->GetValue(i);
		int closest_pattern = 0;
		int closest_dist = 100000000;
		for (vector<int>::iterator it = allpatterns->begin(); it != allpatterns->end(); it++)
		{
			int dist = 0;
			int xor = (*it) ^ tag;
			for (int k = 0; k < 8 * sizeof(int); k++)
			{
				dist += (xor & 1);
				xor >>= 1;
			}
			if (dist < closest_dist)
			{
				closest_dist = dist;
				closest_pattern = *it;
			}
		}

		max_dist = max(max_dist, closest_dist);
		newScalars->InsertValue(i, closest_pattern);
	}
	vtkMesh->GetPointData()->AddArray(newScalars);
	printf("Maximum distance seen is %d\n", max_dist);

	return allpatterns;
}

void PatternGraph::UpdateGraphWithInput(int ts1, vtkPolyData* mesh1, int ts2, vtkPolyData* mesh2)
{
	// find all patterns for the two meshes
	vector<int>* allpatterns1 = FindAllPatterns(mesh1);
	vector<int>* allpatterns2 = FindAllPatterns(mesh2);

	// add nodes for the first mesh
	vtkIntArray* mesh1patterns = (vtkIntArray*) mesh1->GetPointData()->GetArray("PATTBL");
	for (int i = 0; i < mesh1->GetNumberOfPoints(); i++)
	{
		PatternNode cn(ts1, mesh1patterns->GetValue(i));

		if (nodes.find(cn) == nodes.end())
		{
			patternvertex_t v = boost::add_vertex(graph);
			graph[v].timestep = ts1;
			graph[v].pattern = mesh1patterns->GetValue(i);
			graph[v].size = 0;
			nodes[cn] = v;
		}
		graph[nodes[cn]].size++;
	}

	// add nodes for the second mesh
	vtkIntArray* mesh2patterns = (vtkIntArray*) mesh2->GetPointData()->GetArray("PATTBL");
	for (int i = 0; i < mesh2->GetNumberOfPoints(); i++)
	{
		PatternNode cn(ts2, mesh2patterns->GetValue(i));

		if (nodes.find(cn) == nodes.end())
		{
			patternvertex_t v = boost::add_vertex(graph);
			graph[v].timestep = ts2;
			graph[v].pattern = mesh2patterns->GetValue(i);
			graph[v].size = 0;
			nodes[cn] = v;
		}
		graph[nodes[cn]].size++;
	}

	// get links weights map
	boost::property_map<PatternGraph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, graph);

	// add links
	vtkIntArray* m12 = (vtkIntArray*) mesh1->GetPointData()->GetArray("FMTBL");
	for (int id1 = 0; id1 < mesh1->GetNumberOfPoints(); id1++)
	{
		int id2 = m12->GetValue(id1);
		int comp1 = mesh1patterns->GetValue(id1);
		int comp2 = mesh2patterns->GetValue(id2);

		patternvertex_t v1 = nodes[PatternNode(ts1, comp1)];
		patternvertex_t v2 = nodes[PatternNode(ts2, comp2)];

		bool bRet;
		patternedge_t e;
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

void PatternGraph::WriteFile(string filename)
{
	// initialize Graph
	PatternGraph_o ng;

	// copy the vertices
	map<PatternNode, patternvertex_o> nodes_o;
	boost::graph_traits<PatternGraph_t>::vertex_iterator vi, vend;
	for (boost::tie(vi, vend) = vertices(graph); vi != vend; ++vi) 
	{
		patternvertex_o v = boost::add_vertex(ng);
		ng[v].timestep = graph[*vi].timestep;
		ng[v].pattern = graph[*vi].pattern;
		ng[v].size = graph[*vi].size;
		nodes_o[ng[v]] = v;
	}
	
	// copy the edges
	boost::property_map<PatternGraph_t, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, graph);
	boost::graph_traits<PatternGraph_t>::edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei) 
	{
		patternvertex_t gv1 = source(*ei, graph);
		patternvertex_t gv2 = target(*ei, graph);
		patternvertex_o v1 = nodes_o[graph[gv1]];
		patternvertex_o v2 = nodes_o[graph[gv2]];
		boost::add_edge(v1, v2, get(weightmap, *ei), ng);
	}
		
	// write the graph to a file
	std::ofstream output(filename.c_str(), std::ios::out);
	boost::write_graphviz(output, ng, PatternGraph_writer(&ng), PatternGraph_writer(&ng));
	output.close();
}