#pragma once
#define NULL 0
#define PI 3.1415926535897932384626433832795
#define THRESHOLD_T_OF_TWO_POINTS 1e-10 //threhold of t-value to distinguish two different intersections
#define LONDIM 0
#define LATDIM 1
#define RDIM 2
#define XDIM 0
#define YDIM 1
#define ZDIM 2
#include <string>
#include <vector>
using namespace std;
struct INDEX_LENGTH {
	int index[3]; //the three indice of a voxel, longitudinal, latitudinal, and radial index, respectively
	long id;
	double intercept;//the length of a segment
	double t_vec_first_pt;//the first intersection with a voxel along a ray
};
struct RAY_VOX {
	int rayid;
	int numOfVoxel;
	INDEX_LENGTH* vox_itcs;
};
struct RAYPOINT {
	double t;
	double xyz[3];
	double llr[3];
};
struct NEXTVISIT {
	double t_vector;
	int next_visit[3];//the next visitor
};
int lower_bound_index(double* vals, int num, double val);
int upper_bound_index(double* vals, int num, double val);
int t_vec_to_sphere(const double pt1[3], const double pt2[3], const double rCoord, double t_vector[2]);
int t_vec_to_latitude_surface(const double pt1[3], const double pt2[3], const double latCoord, double t_vector[2]);
int t_vec_to_longitude_surface(const double pt1[3], const double pt2[3], const double lonCoord, double& t_vector);
void XYZ_from_t_vec(const double pt1[3], const double pt2[3], const double t, double result[3]);
void  LLR_to_XYZ(const double oldCoord[3], double newCoord[3]);
void XYZ_to_LLR(const double oldCoord[3], double newCoord[3]);
bool is_between_lon_range(double lower_boundary, double upper_boundary, double lon);
bool is_a_valid_ray(double pt1[3], double pt2[3], double minCoords[3], double maxCoords[3], double* t_firstLast);
void get_list_from_filter(const string& filterStr, int defaultNum, vector<string>& list);