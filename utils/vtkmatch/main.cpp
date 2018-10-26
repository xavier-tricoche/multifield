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
#include <vtkCurvatures.h>
#include <vtkDoubleArray.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <gmmreg_api.h>
#include <gmmreg_utils.h>

using namespace std;

int ts_count = 0;
int fs_count = 0;
char** meshes_files;

int cint(double x){

	double fractpart, intpart;
	fractpart = modf (x , &intpart);

	if (fractpart>=.5)
		return x>=0?ceil(x):floor(x);
	else
		return x<0?ceil(x):floor(x);
}

Nrrd *readNrrd(const char* filename)
{
    Nrrd *nin = nrrdNew();
    if (nrrdLoad(nin, filename, NULL)) {
        printf("Error reading %s\n",filename);
                return NULL;
    }
    return nin;
}

void FillTable(vtkPolyData* m, string tbl, int value)
{
	vtkPointData* pointsData = m->GetPointData();
	vtkSmartPointer<vtkIntArray> newScalars = vtkSmartPointer<vtkIntArray>::New();
	newScalars->SetName(tbl.c_str());
	for (vtkIdType i = 0; i < m->GetNumberOfPoints(); i++)
	{
		newScalars->InsertValue(i, value);
	}
	pointsData->AddArray(newScalars);
}

void replace(string &str, const string &find_what, const string &replace_with)
{
	string::size_type pos=0;
	while((pos=str.find(find_what, pos))!=string::npos)
	{
		str.erase(pos, find_what.length());
		str.insert(pos, replace_with);
		pos+=replace_with.length();
	}
}

void WriteINI(string model, string scene, string cmodel, string omodel, string config, string oconfig)
{
	std::fstream input(config.c_str(), std::ios::in);
	std::fstream output(oconfig.c_str(), std::ios::out);
	
	string buffer;
	while (std::getline(input, buffer))
	{
		replace(buffer, "MODEL.txt", model);
		replace(buffer, "SCENE.txt", scene);
		replace(buffer, "CTRLPTS.txt", cmodel);
		replace(buffer, "TRANSFORMED.txt", omodel);		
		
		output << buffer << '\n';
	}
	
	input.close();
	output.close();
}

void f(const vnl_matrix<double>& model,
    const vnl_matrix<double>& scene, double threshold,
    vnl_matrix<double>& extracted_model,
    vnl_matrix<double>& extracted_scene,
	vnl_matrix<int>& pairs) 
{
	vnl_matrix<double> dist;
	ComputeSquaredDistanceMatrix(model, scene, dist);
	pick_indices(dist, pairs, threshold*threshold);
	std::cout << "distance threshold : " << threshold << std::endl;
	int j, n = pairs.cols();
	int d = model.cols();
	extracted_model.set_size(n,d);
	extracted_scene.set_size(n,d);
	std::cout << "# of matched point pairs : " << n << std::endl;
	for (j=0; j<n; ++j) {
		extracted_model.set_row(j,model.get_row(pairs(0,j)));
	}
	for (j=0; j<n; ++j) {
		extracted_scene.set_row(j,scene.get_row(pairs(1,j)));
	}
}

void GetPairs(string model_file, string scene_file, double msp, vnl_matrix<int>& pairs)
{
	std::ifstream infile1(model_file.c_str());
	vnl_matrix<double> model;
	model.read_ascii(infile1);

	std::ifstream infile2(scene_file.c_str());
	vnl_matrix<double> scene;
	scene.read_ascii(infile2);

	vnl_matrix<double> extracted_model, extracted_scene;
	f(model, scene, msp, extracted_model, extracted_scene, pairs);
}

void OutMesh(vtkPolyData* m, set<vtkIdType> l)
{
	vtkSmartPointer<vtkIntArray> newScalars1 = vtkSmartPointer<vtkIntArray>::New();
	newScalars1->SetName("TAGTBL");
	for (vtkIdType i = 0; i < m->GetNumberOfPoints(); i++)
	{
		if (l.find(i) != l.end())
			newScalars1->InsertValue(i, 15);
		else
			newScalars1->InsertValue(i, 0);
	}
	m->GetPointData()->AddArray(newScalars1);
	
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter2 = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter2->SetFileName("bla.vtk");
	vtkMeshWriter2->SetInput(m);
	vtkMeshWriter2->Write();
}

void WriteControlPoints(vtkPolyData* m, string fname)
{
	vtkSmartPointer<vtkCurvatures> curve1 = vtkSmartPointer<vtkCurvatures>::New();
	curve1->SetInput(m);
	curve1->SetCurvatureTypeToMean();
	curve1->Update();
	vtkPolyData* om = curve1->GetOutput();
	vtkPointData* pointsData = om->GetPointData();
	vtkDoubleArray* arr = (vtkDoubleArray*) pointsData->GetArray("Mean_Curvature");
	cout << arr->GetNumberOfTuples() << "\n";
	
	// find range
	double range[2];
	arr->GetRange(range);
	double minv = 0.0;
	double maxv = 0.001 * max(range[0] * -1, range[1]);
	
	// create histogram
	int bins[10000];
	memset(bins, 0, sizeof(int) * 10000);
	for (vtkIdType i = 0; i < m->GetNumberOfPoints(); i++)
	{
		int c = cint(9999 * (abs(arr->GetValue(i)) - minv) / (maxv - minv));
		c = min(9999, max(0, c));
		bins[c] ++;
	}
	//for (int i = 0; i < 10000; i++)
		//cout << bins[i] <<",";
	double* prob = (double*) malloc(m->GetNumberOfPoints() * sizeof(double));
	double total = 0;
	for (vtkIdType i = 0; i < m->GetNumberOfPoints(); i++)
	{
		int c = cint(9999 * (abs(arr->GetValue(i)) - minv) / (maxv - minv));
		c = min(9999, max(0, c));
		
		prob[i] = 0.0001 / bins[c];
		total += prob[i];
	}
	cout << total << "\n";
	
	// find distributions
	fstream filestr1;
	filestr1.open (fname.c_str(), fstream::in | fstream::out | fstream::app);
	if (!filestr1.is_open()) cout << "Could not open file\n";
	int number_needed = 500;
	int count = 0;
	double sum = 0.0;
	set<vtkIdType> l;
	for (vtkIdType i = 0; i < m->GetNumberOfPoints(); i++)
	{
		//sum += prob[i];
		//if (sum > 1.0/number_needed)
		if (i%(m->GetNumberOfPoints()/number_needed) == 0)
		{
			sum = 0.0;
			l.insert(i);
			double* point = m->GetPoints()->GetPoint(i);
			filestr1 << "   " << point[0];
			filestr1 << "   " << point[1];
			filestr1 << "   " << point[2] << "\n";
			count ++;
		}
	}
	cout << count << "\n";
	filestr1.close();
	
	OutMesh(m,l);
}

void MatchMeshes(vtkPolyData* m1, vtkPolyData* m2, char* config, char* method, double msp)
{
	// write first mesh to file
	fstream filestr1;
	filestr1.open ("tmp1.txt", fstream::in | fstream::out | fstream::app);
	if (!filestr1.is_open()) cout << "Could not open file\n";
	for (vtkIdType i = 0; i < m1->GetNumberOfPoints(); i++)
	{
		double* point = m1->GetPoints()->GetPoint(i);
		filestr1 << "   " << point[0];
		filestr1 << "   " << point[1];
		filestr1 << "   " << point[2] << "\n";
	}	
	filestr1.close();
	
	// write second mesh to file
	fstream filestr2;
	filestr2.open ("tmp2.txt", fstream::in | fstream::out | fstream::app);
	if (!filestr2.is_open()) cout << "Could not open file\n";
	for (vtkIdType i = 0; i < m2->GetNumberOfPoints(); i++)
	{
		double* point = m2->GetPoints()->GetPoint(i);
		filestr2 << "   " << point[0];
		filestr2 << "   " << point[1];
		filestr2 << "   " << point[2] << "\n";
	}	
	filestr2.close();
	
	// write control points
	WriteControlPoints(m1, "tmp4.txt");
	
	// write the ini
	vnl_matrix<int> pairs1;
	WriteINI("tmp2.txt", "tmp1.txt", "tmp4.txt", "tmp3.txt", string(config), "tmp.ini");
	gmmreg_api("tmp.ini", method);
	//GetPairs("tmp3.txt", "tmp1.txt", msp, pairs1); // 0 for m2, 1 for m1 => m1 MATCHTBL1
	
	// write table for m1
	/*int n1 = pairs1.cols();
	vtkPointData* pointsData1 = m1->GetPointData();
	vtkSmartPointer<vtkIntArray> newScalars1 = vtkSmartPointer<vtkIntArray>::New();
	newScalars1->SetName("MATCHTBL1");
	for (vtkIdType i = 0; i < m1->GetNumberOfPoints(); i++)
	{
		newScalars1->InsertValue(i, -1);
	}
	for (int j = 0; j < n1; ++j)
	{
		newScalars1->SetValue(pairs1(1,j), pairs1(0,j));
	}
	pointsData1->AddArray(newScalars1);*/
	
	// write the ini
	//vnl_matrix<int> pairs2;
	//WriteINI("tmp1.txt", "tmp2.txt", "tmp3.txt", string(config), "tmp.ini");
	/*gmmreg_api("tmp.ini", method);
	GetPairs("tmp3.txt", "tmp2.txt", msp, pairs2); // 0 for m1, 1 for m2 => m2 MATCHTBL0
	
	// write table for m2
	int n2 = pairs2.cols();
	vtkPointData* pointsData2 = m2->GetPointData();
	vtkSmartPointer<vtkIntArray> newScalars2 = vtkSmartPointer<vtkIntArray>::New();
	newScalars2->SetName("MATCHTBL0");
	for (vtkIdType i = 0; i < m2->GetNumberOfPoints(); i++)
	{
		newScalars2->InsertValue(i, -1);
	}
	for (int j = 0; j < n2; ++j)
	{
		newScalars2->SetValue(pairs2(1,j), pairs2(0,j));
	}
	pointsData2->AddArray(newScalars2);*/
	
}

int main (int argc, char *argv[])
{
	if (argc != 6)
	{
		printf("Example: vtkmatch mesh1.vtk mesh2.vtk config.ini TPS_L2 0.002 \n");
		return 0;
	}
	
	cout << "========\n";
	cout << "vtkmatch\n";
	cout << "========\n";
	
	
	// match the meshes
	vtkPolyData* m1;
	vtkPolyData* m2;
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader1;
	vtkSmartPointer<vtkPolyDataReader> vtkMeshReader2;
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter1;
	vtkSmartPointer<vtkPolyDataWriter> vtkMeshWriter2;
	
	// read mesh 1
	vtkMeshReader1 = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader1->SetFileName(argv[1]);
	vtkMeshReader1->Update();
	m1 = vtkMeshReader1->GetOutput();
	
	// read mesh 2
	vtkMeshReader2 = vtkSmartPointer<vtkPolyDataReader>::New();
	vtkMeshReader2->SetFileName(argv[2]);
	vtkMeshReader2->Update();
	m2 = vtkMeshReader2->GetOutput();
		
	// add empty tables 
	FillTable(m1, "MATCHTBL1", -1);
	FillTable(m2, "MATCHTBL0", -1);
	
	// do the matching process
	MatchMeshes(m1, m2, argv[3], argv[4], atof(argv[5]));
			
	// write the file m1
	vtkMeshWriter1 = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter1->SetFileName(argv[1]);
	vtkMeshWriter1->SetInput(m1);
	//vtkMeshWriter1->Write();
	
	// write the file m2
	vtkMeshWriter2 = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkMeshWriter2->SetFileName(argv[2]);
	vtkMeshWriter2->SetInput(m2);
	//vtkMeshWriter2->Write();
	
	return 0;
}