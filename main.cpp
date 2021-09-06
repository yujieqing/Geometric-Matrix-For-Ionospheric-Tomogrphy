#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "common.h"
#include "SIVT.h"
#include "AP.h"
#include <iomanip>
#include <stdlib.h>
#include <iostream>
#include<list>
#include <chrono>
#include<map>
using namespace std;

bool cmpfunction2(INDEX_LENGTH& a, INDEX_LENGTH & b) {
	return (a.id < b.id);
}
void read_rays_coords_from_stecfile(const char* filename, double*& rays, int& numOfRay)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function: read coordinates for rays from a file, and output the rays in an two-dimensional array, each row of the file records the two endpoints of the ray in longitude, latitude, and radial coordinates
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	ifstream ifile(filename);
	vector<double> vals;
	string linebuf;
	//int iLine = 0;
	while (getline(ifile, linebuf))
	{
		if (linebuf[0]=='#') continue;
		stringstream ss(linebuf);
		string tmp;
		int iCol = 0;
		while (getline(ss, tmp, ','))
		{
			if((iCol>3 && iCol<7) || (iCol>8 && iCol<12))
				vals.push_back(atof((tmp.c_str())));
			iCol++;
		}
		//iLine++;
	}

	numOfRay = vals.size() / 6;

	if (rays) delete[] rays;
	rays = new double[numOfRay * 6];
	for (int i = 0; i < numOfRay; i++)
	{
		int j = i;// +begin;
		rays[i * 6] = vals[j * 6];
		rays[i * 6 + 1] = vals[j * 6 + 1];
		rays[i * 6 + 2] = vals[j * 6 + 2];
		rays[i * 6 + 3] = vals[j * 6 + 3];
		rays[i * 6 + 4] = vals[j * 6 + 4];
		rays[i * 6 + 5] = vals[j * 6 + 5];
	}
	ifile.close();

}

long segments_compuation_by_SIVT(double* rays, int numOfRay, double* split_surf_coords[3], int num_of_split_surfs[3], std::list<RAY_VOX*>& segments)
//compute segements for a number of rays using SIVT method
{
	int count = 0;
	for (int i = 0; i < numOfRay; i++)
	{
		INDEX_LENGTH* newitc = NULL;
		int num = segments_of_a_ray_by_AP(rays + i * 6, rays + i * 6 + 3, split_surf_coords, num_of_split_surfs, newitc);
		if (num > 0) {
			RAY_VOX* new_ray_vox = new RAY_VOX();
			new_ray_vox->rayid = i;
			new_ray_vox->numOfVoxel = num;
			new_ray_vox->vox_itcs = newitc;
			segments.push_back(new_ray_vox);
			count++;
		}
	}
	return count;
}

void segments_computation_by_tradition(double* rays, int numOfRay, double* split_surf_coords[3], int num_of_split_surfs[3], std::list<RAY_VOX*>& segments)
{
	for (int i = 0; i < numOfRay; i++)
	{
		vector<INDEX_LENGTH*> segments_of_a_ray;
		int num = segments_of_a_ray_by_tradition(rays + i * 6, rays + i * 6 + 3, split_surf_coords, num_of_split_surfs, segments_of_a_ray);
		if (num > 0) {
			RAY_VOX* new_ray_vox = new RAY_VOX();
			new_ray_vox->rayid = i;
			new_ray_vox->numOfVoxel = num;
			INDEX_LENGTH* segments_of_a_ray_Ptr = new INDEX_LENGTH[segments_of_a_ray.size()];
			new_ray_vox->vox_itcs = segments_of_a_ray_Ptr;
			// the following codes can be removed when comparing with ray tracing method. They just copy data to a new structure
			for (int j = 0; j < segments_of_a_ray.size(); j++)
			{
				segments_of_a_ray_Ptr[j].intercept = segments_of_a_ray[j]->intercept;
				segments_of_a_ray_Ptr[j].t_vec_first_pt = segments_of_a_ray[j]->t_vec_first_pt;
				segments_of_a_ray_Ptr[j].index[0] = segments_of_a_ray[j]->index[0];
				segments_of_a_ray_Ptr[j].index[1] = segments_of_a_ray[j]->index[1];
				segments_of_a_ray_Ptr[j].index[2] = segments_of_a_ray[j]->index[2];
				segments_of_a_ray_Ptr[j].id = segments_of_a_ray[j]->id;
				delete[] segments_of_a_ray[j];
				segments_of_a_ray[j] = NULL;
			}			
			segments.push_back(new_ray_vox);
		}

	}
}
void segments_computation_by_AP(double* rays, int numOfRay, double* split_surf_coords[3], int num_of_split_surfs[3], std::list<RAY_VOX*>& segments)
//compute segements for a number of rays using AP method, see, Hong J, Kim Y H, Chung J K, et al. (2017) Tomography reconstruction of ionospheric electron density with empirical orthonormal functions using Korea GNSS network. J Astron Space Sci 34(1): 7-17

{
	for (int i = 0; i < numOfRay; i++)
	{
		vector<INDEX_LENGTH*> segments_of_a_ray;
		int num = segments_of_a_ray_by_AP(rays + i * 6, rays + i * 6 + 3, split_surf_coords, num_of_split_surfs, segments_of_a_ray);
		if (num > 0) {
			RAY_VOX* new_ray_vox = new RAY_VOX();
			new_ray_vox->rayid = i;
			new_ray_vox->numOfVoxel = num;
			INDEX_LENGTH* segments_of_a_ray_Ptr = new INDEX_LENGTH[segments_of_a_ray.size()];
			new_ray_vox->vox_itcs = segments_of_a_ray_Ptr;
			// the following codes can be removed when comparing with ray tracing method. They just copy data to a new structure
			for (int j = 0; j < segments_of_a_ray.size(); j++)
			{
				segments_of_a_ray_Ptr[j].intercept = segments_of_a_ray[j]->intercept;
				segments_of_a_ray_Ptr[j].t_vec_first_pt = segments_of_a_ray[j]->t_vec_first_pt;
				segments_of_a_ray_Ptr[j].index[0] = segments_of_a_ray[j]->index[0];
				segments_of_a_ray_Ptr[j].index[1] = segments_of_a_ray[j]->index[1];
				segments_of_a_ray_Ptr[j].index[2] = segments_of_a_ray[j]->index[2];
				segments_of_a_ray_Ptr[j].id = segments_of_a_ray[j]->id;
				delete[] segments_of_a_ray[j];
				segments_of_a_ray[j] = NULL;
			}
			segments.push_back(new_ray_vox);
		}

	}
}
void save_segments_by_matrix(const char* filename, std::list<RAY_VOX*>& segments,long voxDims)
// save the resutant segments in form of matix 
{
	ofstream ofile(filename);
	for (list<RAY_VOX*>::const_iterator it = segments.begin(); it != segments.end(); it++)
	{
		RAY_VOX* a_ray = *it;
		sort(a_ray->vox_itcs, a_ray->vox_itcs + a_ray->numOfVoxel, cmpfunction2);
		//ofile << index;
		long beginID = 0, endID = 0;
		for (long j = 0; j<a_ray->numOfVoxel; j++) {

			endID = a_ray->vox_itcs[j].id;
			for (int i = beginID; i < endID; i++)
			{
				if (i == voxDims - 1)
					ofile << 0;
				else
					ofile << 0 << ",";
			}
			if (endID == voxDims - 1)
			ofile << a_ray->vox_itcs[j].intercept;
			else
				ofile << a_ray->vox_itcs[j].intercept << ",";

			beginID = endID+1;
		}
		for (int i = beginID; i < voxDims; i++)
		{
			if (i == voxDims - 1)
				ofile << 0;
			else
				ofile << 0 << ",";
		}
		ofile << endl;
	}
	ofile.close();

}

void save_segments_by_cent_size_val(const char* filename, std::list<RAY_VOX*>& segments, double** split_surf_coords, int num_of_split_surfs[3])
// save the segments in center, size and length format.Each row of the file is a segment in a voxel. There are seven columns in a row
// The first six columns indicate the center and size of a voxel by spherical coordinates(in radians)
// the last column indicates the length of the segment associated with the voxel
{
	int numOfRay = segments.size();
	int dimVox10 = (num_of_split_surfs[LATDIM] - 1)*(num_of_split_surfs[LONDIM] - 1);
	ofstream ofile(filename);
	for (list<RAY_VOX*>::const_iterator it = segments.begin(); it != segments.end(); it++)
	{

		RAY_VOX* a_ray = *it;
		sort(a_ray->vox_itcs, a_ray->vox_itcs + a_ray->numOfVoxel, cmpfunction2);

		for (int j = 0; j < a_ray->numOfVoxel; j++)
		{
			long id = a_ray->vox_itcs[j].id;
			long rIndex = id / dimVox10;
			long restIndex = id % dimVox10;
			long latIndex = restIndex / (num_of_split_surfs[LONDIM] - 1);
			long lonIndex = restIndex % (num_of_split_surfs[LONDIM] - 1);
			ofile << (split_surf_coords[LONDIM][lonIndex + 1] + split_surf_coords[LONDIM][lonIndex]) / 2 << ", "
				<< (split_surf_coords[LATDIM][latIndex + 1] + split_surf_coords[LATDIM][latIndex]) / 2 << ", "
				<< (split_surf_coords[RDIM][rIndex + 1] + split_surf_coords[RDIM][rIndex]) / 2 << ", "
				<< (split_surf_coords[LONDIM][lonIndex + 1] - split_surf_coords[LONDIM][lonIndex]) << ", "
				<< (split_surf_coords[LATDIM][latIndex + 1] - split_surf_coords[LATDIM][latIndex]) << ", "
				<< (split_surf_coords[RDIM][rIndex + 1] - split_surf_coords[RDIM][rIndex]) << ", "
				<< a_ray->vox_itcs[j].intercept	<< endl;
		}
	}
	ofile.close();
}
void save_segments_by_sparse_matrix(const char* filename, std::list<RAY_VOX*>& segments)
// Each row of the file is a segment in a voxel. There are three colomns in a row, i.e., ray ID, voxel ID, segment's length
{
	int numOfRay = segments.size();
	ofstream ofile(filename);
	for (list<RAY_VOX*>::const_iterator it = segments.begin(); it != segments.end(); it++)
	{
		RAY_VOX* a_ray = *it;
		sort(a_ray->vox_itcs, a_ray->vox_itcs + a_ray->numOfVoxel, cmpfunction2);
		for (int j = 0; j < a_ray->numOfVoxel; j++)
		{
			long id = a_ray->vox_itcs[j].id;
			ofile << a_ray->rayid<<","<< id << ", "<< a_ray->vox_itcs[j].intercept << endl;
		}
	}
	ofile.close();
}
void save_segments_by_indice(const char* filename, std::list<RAY_VOX*>& segments, double** split_surf_coords, int num_of_split_surfs[3])
// ray ID, voxel ID, longitudinal index, latitudinal index, radial index,  segment's length
{
	int numOfRay = segments.size();
	int dimVox10 = (num_of_split_surfs[LATDIM] - 1)*(num_of_split_surfs[LONDIM] - 1);
	ofstream ofile(filename);
	for (list<RAY_VOX*>::const_iterator it = segments.begin(); it != segments.end(); it++)
	{

		RAY_VOX* a_ray = *it;
		sort(a_ray->vox_itcs, a_ray->vox_itcs + a_ray->numOfVoxel, cmpfunction2);
		for (int j = 0; j < a_ray->numOfVoxel; j++)
		{
			long id = a_ray->vox_itcs[j].id;
			
			long rIndex = id / dimVox10;
			long restIndex = id % dimVox10;
			long latIndex = restIndex / (num_of_split_surfs[LONDIM] - 1);
			long lonIndex = restIndex % (num_of_split_surfs[LONDIM] - 1);
			ofile << a_ray->rayid << ","
				<< id << ", "
				<< lonIndex << ", "
				<< latIndex << ", "
				<< rIndex << ", "
				<< a_ray->vox_itcs[j].intercept << endl;
		}
	}
	ofile.close();
}
void readConf(const string& configFile,std::map<int,double*>& voxModels, string& rayFile,string& outputPath)
{
	ifstream conf(configFile.c_str());
	string linebuf;
	//string ray_file_filter;
	int iVox = -1;
	while (getline(conf, linebuf))
	{
		if (linebuf[0] == '#') continue;
		std::size_t found = linebuf.find('=', 0);
		if (found != std::string::npos)
		{
			string key = linebuf.substr(0, found);
			if (key.compare("rayFile") == 0)
			{
				rayFile = linebuf.substr(found + 1, std::string::npos);
			}
			//if (key.compare("stecFilter") == 0)
			//{
			//	ray_file_filter = linebuf.substr(found + 1, std::string::npos);
			//}
			else if (key.compare("outPath") == 0)
			{
				outputPath = linebuf.substr(found + 1, std::string::npos);
			}
			else if (key.compare("VoxModel") == 0)
			{
				string suffix = linebuf.substr(found + 1, std::string::npos);
				stringstream ss(suffix);
				string tmp;
				double* voxModelInfo = new double[9];
				ss >> iVox >> voxModelInfo[0]
					>> voxModelInfo[1]
					>> voxModelInfo[2]
					>> voxModelInfo[3]
					>> voxModelInfo[4]
					>> voxModelInfo[5]
					>> voxModelInfo[6]
					>> voxModelInfo[7]
					>> voxModelInfo[8];
				voxModelInfo[2] = voxModelInfo[2] * 1000 + 6371000; // from km to m, and from height to radius
				voxModelInfo[5] = voxModelInfo[5] * 1000 + 6371000;// from km to m, and from height to radius
				voxModelInfo[8] *= 1000;// from height to radius
				if (voxModels.find(iVox) != voxModels.end())
					delete[] voxModelInfo;
				else
					voxModels.insert(std::pair<int, double*>(iVox, voxModelInfo));

			}
		}
	}
}
int main()
{

	string configFile = "input\\Conf.txt";
	string rayFile, outputPath;
	std::map<int, double*> voxModels;
	readConf(configFile, voxModels, rayFile, outputPath);	
	int numOfRay = 0;
	double* rays = NULL;
	read_rays_coords_from_stecfile(rayFile.c_str(), rays, numOfRay);
	char runtimeFile[2048]; // = "output/TimeCost.txt";
	sprintf_s(runtimeFile, 2048, "%s\\TimeCost.txt", outputPath.c_str());
	ofstream ofile(runtimeFile);
	ofile << "numOfRays,tnumOfValidRays,iVoxModel, t_VIST, t_AP, nSplitLon, nSplitLat, nSplitR, nVox" << endl;
	ofile.setf(ios::fixed, ios::floatfield);
	for (auto it = voxModels.begin(); it != voxModels.end(); it++)
	{
		//get the spliting surfaces information
		double* voxModelInfo = it->second;
		int iVoxModel = it->first;
		int num_of_split_surfs[3];
		num_of_split_surfs[LONDIM] = int((voxModelInfo[3] - voxModelInfo[0]) / voxModelInfo[6] + 0.5) + 1;
		num_of_split_surfs[LATDIM] = int((voxModelInfo[4] - voxModelInfo[1]) / voxModelInfo[7] + 0.5) + 1;
		num_of_split_surfs[RDIM] = int((voxModelInfo[5] - voxModelInfo[2]) / voxModelInfo[8] + 0.5) + 1;
		double* split_surf_coords[3];
		for (int iDim = 0; iDim < 2; iDim++)
		{
			split_surf_coords[iDim] = new double[num_of_split_surfs[iDim]];

			for (int i = 0; i < num_of_split_surfs[iDim]; i++)
			{
				split_surf_coords[iDim][i] = (voxModelInfo[iDim] + i * voxModelInfo[6 + iDim]) / 180 * PI;
			}
		}
		split_surf_coords[RDIM] = new double[num_of_split_surfs[RDIM]];
		for (int i = 0; i < num_of_split_surfs[RDIM]; i++)
		{
			split_surf_coords[RDIM][i] = (voxModelInfo[RDIM] + i * voxModelInfo[6 + RDIM]);
		}


		//SIVT approach
		std::cout << "Starting SIVT approach"<<endl;
		std::list<RAY_VOX*> segments;
		chrono::high_resolution_clock::time_point tt1 = chrono::high_resolution_clock::now();
		long numValid = segments_compuation_by_SIVT(rays, numOfRay, split_surf_coords, num_of_split_surfs, segments);
		chrono::high_resolution_clock::time_point tt2 = chrono::high_resolution_clock::now();
		auto time_span1 = chrono::duration_cast<chrono::microseconds>(tt2 - tt1);
		double msTrace = time_span1.count() / 1000000.0;

		char filename[1024];
		sprintf_s(filename, 1024, "%s\\VIST_centSize%d.txt", outputPath.c_str(), iVoxModel);
		save_segments_by_cent_size_val(filename, segments, split_surf_coords, num_of_split_surfs);
		sprintf_s(filename, 1024, "%s\\VIST_sparse%d.txt", outputPath.c_str(), iVoxModel);
		save_segments_by_sparse_matrix(filename, segments);
		for (auto it = segments.begin(); it != segments.end(); it++)
			delete[] * it;
		segments.clear();

		//AP approach
		std::cout << "Starting AP approach" << endl;
		tt1 = chrono::high_resolution_clock::now();
		segments_computation_by_AP(rays, numOfRay, split_surf_coords, num_of_split_surfs, segments);
		tt2 = chrono::high_resolution_clock::now();
		auto time_span3 = chrono::duration_cast<chrono::microseconds>(tt2 - tt1);
		double msPlane = time_span3.count() / 1000000.0;

		sprintf_s(filename, 1024, "%s\\AP_centSize%d.txt", outputPath.c_str(), iVoxModel);
		save_segments_by_cent_size_val(filename, segments, split_surf_coords, num_of_split_surfs);
		sprintf_s(filename, 1024, "%s\\AP_sparse%d.txt", outputPath.c_str(), iVoxModel);
		save_segments_by_sparse_matrix(filename, segments);
		for (auto it = segments.begin(); it != segments.end(); it++)
			delete[] * it;
		segments.clear();


		//output
		ofile << setprecision(15) << numOfRay << ", \t" << numValid << ",  \t"
			<< iVoxModel << ", \t"
			<< msTrace << ", \t" << msPlane << ", \t"
			<< num_of_split_surfs[LONDIM] << ", \t" << num_of_split_surfs[LATDIM] << ", \t" << num_of_split_surfs[RDIM]
			<< ", \t" << (num_of_split_surfs[LONDIM] - 1) * (num_of_split_surfs[LATDIM] - 1) * (num_of_split_surfs[RDIM] - 1)
			<< endl;
		cout << "numOfRays,\t\tnumOfValidRays,\t\tiVoxModel, \t\tt_VIST, \t\tt_AP, \t\tnSplitLon, \t\tnSplitLat, \t\tnSplitR, \t\tnVox" << endl;
		cout.setf(ios::fixed, ios::floatfield);
		cout << setprecision(15) << numOfRay << ", \t" << numValid << ",  \t"
			<< iVoxModel << ", \t"
			<< msTrace << ", \t" << msPlane << ", \t"
			<< num_of_split_surfs[LONDIM] << ", \t" << num_of_split_surfs[LATDIM] << ", \t" << num_of_split_surfs[RDIM]
			<< ", \t" << (num_of_split_surfs[LONDIM] - 1) * (num_of_split_surfs[LATDIM] - 1) * (num_of_split_surfs[RDIM] - 1)
			<< endl;
		for (int iDim = 0; iDim < 3; iDim++)
		{
			delete[] split_surf_coords[iDim];
		}
		delete[] it->second;
		it->second = NULL;
	}
	ofile.close();

	delete[] rays;

	return 0;
}