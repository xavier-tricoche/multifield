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
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <bitset>
#include <iostream>
#include <fstream>
#include <sstream>
#include <teem/nrrd.h>
#include <vector_types.h>
#include <vector_functions.h>
#include <cutil_inline.h>
#include <cutil_math.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include "PatternGraph.h"
#include "FieldTypeGraph.h"
#include "ClusterTransitionGraph.h"
#include "ClusterTimeTransitionGraph.h"

using namespace std;

int ts_count = 0;
int ts_str = 0;
int ts_end = 0;
int fs_count = 0;
char** meshes_files;

string convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

string convertBinaryInt(int number)
{
	char buf[512];
	itoa(number, buf, 2);
   stringstream ss;//create a stringstream
   ss << buf;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

void BuildPatternGraph()
{
	// build the component graph
	PatternGraph cg;
	vtkPolyData* m1;
	vtkPolyData* m2;
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader1;
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader2;
	for (int i = ts_str; i < ts_end; i++)
	{
		cout << "time step " << i << "\n";
		// read mesh 1
		vtkMeshReader1 = vtkSmartPointer<vtkPolyDataReader>::New();
		vtkMeshReader1->SetFileName(meshes_files[i]);
		vtkMeshReader1->Update();
		m1 = vtkMeshReader1->GetOutput();
		
		// read mesh 2
		vtkMeshReader2 = vtkSmartPointer<vtkPolyDataReader>::New();
		vtkMeshReader2->SetFileName(meshes_files[i + 1]);
		vtkMeshReader2->Update();
		m2 = vtkMeshReader2->GetOutput();

		// update the graph information
		cg.UpdateGraphWithInput(i, m1, i + 1, m2);
	}
	cg.WriteFile("patterngraph.dot");
}

void BuildFieldTypeGraph(int nvar)
{
	// build the component graph
	FieldTypeGraph cg(nvar);
	vtkPolyData* m1;
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader1;
	for (int i = ts_str; i <= ts_end; i++)
	{
		cout << "time step " << i << "\n";

		vtkMeshReader1 = vtkSmartPointer<vtkPolyDataReader>::New();
		vtkMeshReader1->SetFileName(meshes_files[i]);
		vtkMeshReader1->Update();
		m1 = vtkMeshReader1->GetOutput();

		// update the graph information
		cg.UpdateGraphWithInput2(m1);
	}
	cg.BuildGraph();
	cg.WriteFile("fieldtypegraph.dot");
}

void BuildLDADocuments(int nvar)
{
	vtkPolyData* m1;
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader1;
	for (int i = ts_str; i <= ts_end; i++)
	{
		cout << "time step " << i << "\n";

		vtkMeshReader1 = vtkSmartPointer<vtkPolyDataReader>::New();
		vtkMeshReader1->SetFileName(meshes_files[i]);
		vtkMeshReader1->Update();
		m1 = vtkMeshReader1->GetOutput();

		// file name
		string filestr = convertInt(i);
		while(filestr.length() < 4)
		{
			filestr = string("0"+filestr);
		}
		filestr = "Z:\\nabla\\Samer\\Datasets\\Multifield\\LDA_documents\\LDA_document_" + filestr + ".txt";

		// create file
		ofstream myfile;
		myfile.open(filestr.c_str());

		// write file contents
		vtkIntArray* arr = (vtkIntArray*) m1->GetPointData()->GetArray("TAGTBL");
		string buffer;
		char replac[128];
		for (vtkIdType i = 0; i < m1->GetNumberOfPoints(); i++)
		{
			int tag = arr->GetValue(i);
			buffer = convertBinaryInt(tag);
			while(buffer.length() < 2*nvar)
			{
				buffer = string("0"+buffer);
			}

			int j = 0;
			while (j < buffer.length())
			{
				switch (buffer[j])
				{
				case '0':
					replac[j] = 'a';
					break;
				case '1':
					replac[j] = 'b';
					break;
				};

				j++;
			}
			replac[j] = '\0';

			myfile << replac << " ";
		}
		

		// close file
		myfile.close();
	}

	
}

void BuildClusterTransitionGraph(int nvar, int ncluster)
{
	cout << "Number of variables is " << nvar << "\n";
	
	// create an array for all patterns counts
	int ntags = std::pow(2.0, 2 * nvar);
	printf("Number of possible tags is %d\n", ntags);
	int* tagcnt = (int*) malloc(ntags * sizeof(int));
	memset(tagcnt, 0, ntags * sizeof(int));
	
	// loop loading meshes and counting
	vtkPolyData* mesh;
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader;
	for (int i = ts_str; i <= ts_end; i++)
	{
		cout << "time step " << i << "\n";
		
		// read mesh
		cout << "Loading " << meshes_files[i] << "\n";
		vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
		vtkMeshReader->SetFileName(meshes_files[i]);
		vtkMeshReader->Update();
		mesh = vtkMeshReader->GetOutput();
		
		// get table for the tags
		vtkIntArray* tagarr = (vtkIntArray*) mesh->GetPointData()->GetArray("TAGTBL");
		
		// set count
		for (vtkIdType j = 0; j < mesh->GetNumberOfPoints(); j++)
		{
			int tag = tagarr->GetValue(j);
			tagcnt[tag]++;
		}
	}
	cout << "Done counting\n";
	
	// now do the clustering 
	printf("Start the clustering\n");
	vector<ClusterNode> clusters;
	ClusterNode fc(nvar);
	fc.ComputeSize(tagcnt);
	clusters.push_back(fc);
	make_heap (clusters.begin(),clusters.end());
	while (clusters.size() < ncluster)
	{
		// find biggest cluster
		pop_heap(clusters.begin(),clusters.end()); 
		ClusterNode c = clusters.back();
		clusters.pop_back();
		
		cout << "split "; c.Print();
		
		// find the bit with highest entropy 
		int bit = c.FindBestBitForSplit(tagcnt);
		
		// do split
		ClusterNode c1(c, bit, 0);
		c1.ComputeSize(tagcnt);
		c1.ComputeEntropy(tagcnt);
		c1.ComputeMostProbableString(tagcnt);
		c1.ComputeMostProbableInt(tagcnt);
		clusters.push_back(c1); 
		push_heap(clusters.begin(), clusters.end());
		//cout << "new 0: "; c1.Print();
		
		ClusterNode c2(c, bit, 1);
		c2.ComputeSize(tagcnt);
		c2.ComputeEntropy(tagcnt);
		c2.ComputeMostProbableString(tagcnt);
		c2.ComputeMostProbableInt(tagcnt);
		clusters.push_back(c2); 
		push_heap(clusters.begin(), clusters.end());
		//cout << "new 1: "; c2.Print();
	}
	
	// print the clusters
	for (vector<ClusterNode>:: iterator it = clusters.begin(); it != clusters.end(); it++)
	{
		(*it).Print();
	}

	// CSV file for the time chart
	ofstream myfile;
	myfile.open ("timechart.csv");
	for (int  j = 0; j < ncluster; j++)
	{
		if (j < ncluster - 1)
			myfile << string(clusters[j].mpclusterstr.rbegin(), clusters[j].mpclusterstr.rend()) << ",";
		else
			myfile << string(clusters[j].mpclusterstr.rbegin(), clusters[j].mpclusterstr.rend()) << "\n";
	}

	// write the result to the files
	for (int i = ts_str; i <= ts_end; i++)
	{
		cout << "time step " << i << "\n";
		
		// read mesh
		cout << "Loading " << meshes_files[i] << "\n";
		vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
		vtkMeshReader->SetFileName(meshes_files[i]);
		vtkMeshReader->Update();
		mesh = vtkMeshReader->GetOutput();
		
		// get table for the tags
		vtkIntArray* tagarr = (vtkIntArray*) mesh->GetPointData()->GetArray("TAGTBL");
		
		// create cluster size array
		int* clssize = (int*) malloc(ncluster * sizeof(int));
		memset(clssize, 0, ncluster * sizeof(int));

		// set count
		vtkPointData* pointsData = mesh->GetPointData();
		vtkSmartPointer<vtkIntArray> newScalars = vtkSmartPointer<vtkIntArray>::New();
		newScalars->SetName("CLUSTBL");
		for (vtkIdType j = 0; j < mesh->GetNumberOfPoints(); j++)
		{
			int tag = tagarr->GetValue(j);
			
			for (int k = 0; k < clusters.size(); k++)
			{
				if (clusters[k].TagMatch(tag))
				{
					newScalars->InsertValue(j, k);
					clssize[k]++;
				}
			}
		}
		pointsData->AddArray(newScalars);

		// write line to file
		for (int  j = 0; j < ncluster; j++)
		{
			if (j < ncluster - 1)
				myfile << clssize[j] << ",";
			else
				myfile << clssize[j] << "\n";
		}

		// free array
		free(clssize);
		
		// write file
		vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
		vtkMeshWriter->SetFileName(meshes_files[i]);
		vtkMeshWriter->SetInput(mesh);
		vtkMeshWriter->Write();
	}
	
	// close file
	myfile.close();

	// free memory
	free(tagcnt);
	

	// build the transition graph
	ClusterTransitionGraph cg;
	ClusterTimeTransitionGraph ctg;
	cg.CreateGraphNodes(&clusters);
	ctg.CreateGraphNodes(&clusters, ts_str, ts_end);
	vtkPolyData* m1;
	vtkPolyData* m2;
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader1;
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader2;
	for (int i = ts_str; i < ts_end; i++)
	{
		cout << "time step " << i << "\n";
		// read mesh 1
		vtkMeshReader1 = vtkSmartPointer<vtkPolyDataReader>::New();
		vtkMeshReader1->SetFileName(meshes_files[i]);
		vtkMeshReader1->Update();
		m1 = vtkMeshReader1->GetOutput();
		
		// read mesh 2
		vtkMeshReader2 = vtkSmartPointer<vtkPolyDataReader>::New();
		vtkMeshReader2->SetFileName(meshes_files[i + 1]);
		vtkMeshReader2->Update();
		m2 = vtkMeshReader2->GetOutput();

		// update the graph information
		cg.UpdateGraphWithInput(i, m1, i + 1, m2, &clusters);
		ctg.UpdateGraphWithInput(i, m1, i + 1, m2, &clusters);
	}
	cg.WriteFile("clustertransitiongraph.dot");
	ctg.WriteFile("clustertimetransitiongraph.dot");

	
}


void MultifieldDatasetCorrect3(int nvar)
{
	int total= 0;
	int* count = (int*) malloc(2 * nvar * sizeof(int));
	for (int bit = 0; bit < 2 * nvar; bit++)
		count[bit] = 0;

	// write the result to the files
	vtkPolyData* mesh;
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader;
	for (int i = ts_str; i <= ts_end; i++)
	{
		cout << "time step " << i << "\n";
		
		// read mesh
		cout << "Loading " << meshes_files[i] << "\n";
		vtkMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
		vtkMeshReader->SetFileName(meshes_files[i]);
		vtkMeshReader->Update();
		mesh = vtkMeshReader->GetOutput();
		
		// get table for the tags
		vtkIntArray* tagarr = (vtkIntArray*) mesh->GetPointData()->GetArray("TAGTBL");
		vtkIntArray* clsarr = (vtkIntArray*) mesh->GetPointData()->GetArray("CLUSTBL");

		// set count
		for (vtkIdType j = 0; j < mesh->GetNumberOfPoints(); j++)
		{
			int tag = tagarr->GetValue(j);
			int cls = clsarr->GetValue(j);

			if (cls != 3) continue;
			if (((tag >> 12) & 3) > 0) continue;

			total++;

			for (int bit = 0; bit < 2 * nvar; bit++)
			{
				int b = ((tag >> bit) & 1);
				count[bit] += b;
			}
			
		}


		for (int bit = 0; bit < 2 * nvar; bit++)
		{
			cout << (count[bit] / float(total)) << "\n";
		}
	}

	cout << "Total is " << total << "\n";
	for (int bit = 0; bit < 2 * nvar; bit++)
	{
		cout << (count[bit] / float(total)) << "\n";
	}
}

int main (int argc, char *argv[])
{
	/*if (argc != 5)
	{
		printf("Example: vtkflowgraph meshes.txt 50 87 network.txt\n");
		return 0;
	}*/

	cout << "============\n";
	cout << "vtkflowgraph\n";
	cout << "============\n";
	
	string directory("Z:\\nabla\\Samer\\Datasets\\Multifield\\multi_mf2\\");
	ts_str = 0;//atoi(argv[2]);
	ts_end = 0;//atoi(argv[3]);
	int nvar = 4;
	int ncluster = 8;

	// read the file for information
	fstream input(directory+"meshes.txt", std::ios::in);
	input >> ts_count;
	input >> fs_count;
	meshes_files = (char**) malloc(ts_count * sizeof(char*));
	for (int i = 0; i < ts_count; i++)
	{
		string buffer;
		input >> buffer;
		buffer = directory + string(buffer);
		meshes_files[i] = new char[buffer.size()+1];
		strcpy(meshes_files[i], buffer.c_str());
	}
	input.close();

	//MultifieldDatasetCorrect3(nvar);
	//BuildClusterTransitionGraph(nvar, ncluster);
	BuildFieldTypeGraph(nvar);

	return 0;
}