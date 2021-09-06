#pragma once
#include"common.h"
#include <math.h>
#include <regex>
#include <sstream> 

int lower_bound_index(double* vals, int num, double val)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function: obtain the lower bound index of a given value in an array
//Input:
//		vals: the array
//		num: the size of the array
//      val: the given value
//Output:
//      return the lower bound index of the given value in the array
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	int first = 0;
	int last = num - 1;
	int index = -1;
	if (val < vals[0])	return -1;
	else if (val >= vals[last])	return num - 1;
	while (first <= last) {
		int mid = (first + last) / 2; // 2
		if (val >= vals[mid] && val < vals[mid + 1])	return mid;
		else if (val < vals[mid])
			last = mid - 1;
		else//  val < vals[index]
			first = mid + 1;
	}
	return -1;
}
int upper_bound_index(double* vals, int num, double val)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function: obtain the upper bound index of a given value in an array
//Input:
//		vals: the array
//		num: the size of the array
//      val: the given value
//Output:
//      return the upper bound index of the given value in the array
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	int first = 0;
	int last = num - 1;
	int index = -1;
	if (val <= vals[0])	return 0;
	else if (val > vals[last])	return -1;
	while (first <= last) {
		int mid = (first + last) / 2; // 2
		if (val > vals[mid] && val <= vals[mid + 1])	return mid + 1;
		else if (val < vals[mid])
			last = mid - 1;
		else//  val > vals[index]
			first = mid + 1;
	}
	return -1;
}
int t_vec_to_sphere(const double pt1[3], const double pt2[3], const double rCoord, double t_vector[2])
{
	double a = pt2[XDIM] - pt1[XDIM];
	double b = pt2[YDIM] - pt1[YDIM];
	double c = pt2[ZDIM] - pt1[ZDIM];
	double A = a * a + b * b + c * c;
	double B = 2 * (a*pt1[XDIM] + b * pt1[YDIM] + c * pt1[ZDIM]);
	double C = pt1[XDIM] * pt1[XDIM] + pt1[YDIM] * pt1[YDIM] + pt1[ZDIM] * pt1[ZDIM] - rCoord * rCoord;
	double bac = B * B - 4 * A*C;
	if (bac < 0) return 0;
	double deta = sqrt(B*B - 4 * A*C);
	double t[2];
	t[0] = (-B - deta) / 2 / A;
	t[1] = (-B + deta) / 2 / A;
	int index = 0;
	if (t[0] >= -THRESHOLD_T_OF_TWO_POINTS && t[0] <= 1 + THRESHOLD_T_OF_TWO_POINTS)
		t_vector[index++] = t[0];
	if (t[1] >= -THRESHOLD_T_OF_TWO_POINTS && t[1] <= 1 + THRESHOLD_T_OF_TWO_POINTS)
		t_vector[index++] = t[1];
	return index;
}
int t_vec_to_latitude_surface(const double pt1[3], const double pt2[3], const double latCoord, double t_vector[2])
//latCoord: in radian
{
	double a = pt2[XDIM] - pt1[XDIM];
	double b = pt2[YDIM] - pt1[YDIM];
	double c = pt2[ZDIM] - pt1[ZDIM];
	double k = tan(latCoord);
	k = k * k;
	double A = k * (a*a + b * b) - c * c;
	double B = 2 * (k*(a*pt1[XDIM] + b * pt1[YDIM]) - c * pt1[ZDIM]);
	double C = k * pt1[XDIM] * pt1[XDIM] + k * pt1[YDIM] * pt1[YDIM] - pt1[ZDIM] * pt1[ZDIM];
	double bac = B * B - 4 * A*C;
	if (bac < 0) return 0;
	double deta = sqrt(bac);
	double t[2];
	t[0] = (-B - deta) / 2 / A;
	t[1] = (-B + deta) / 2 / A;
	int index = 0;
	if (t[0] >= -THRESHOLD_T_OF_TWO_POINTS && t[0] <= 1 + THRESHOLD_T_OF_TWO_POINTS) t_vector[index++] = t[0];
	if (t[1] >= -THRESHOLD_T_OF_TWO_POINTS && t[1] <= 1 + THRESHOLD_T_OF_TWO_POINTS) t_vector[index++] = t[1];
	return index;
}
int t_vec_to_longitude_surface(const double pt1[3], const double pt2[3], const double lonCoord, double& t_vector)
//lonCoord:in radian
{
	double a = pt2[XDIM] - pt1[XDIM];
	double b = pt2[YDIM] - pt1[YDIM];
	double c = pt2[ZDIM] - pt1[ZDIM];
	double cosa = cos(lonCoord);
	double sina = sin(lonCoord);
	double t = (pt1[YDIM] * cosa - pt1[XDIM] * sina) / (a*sina - b * cosa);
	if (t >= -THRESHOLD_T_OF_TWO_POINTS && t <= 1 + THRESHOLD_T_OF_TWO_POINTS) { t_vector = t; return 1; }
	return 0;
}
void XYZ_from_t_vec(const double pt1[3], const double pt2[3], const double t, double result[3])
{
	result[XDIM] = pt1[XDIM] + t * (pt2[XDIM] - pt1[XDIM]);
	result[YDIM] = pt1[YDIM] + t * (pt2[YDIM] - pt1[YDIM]);
	result[ZDIM] = pt1[ZDIM] + t * (pt2[ZDIM] - pt1[ZDIM]);
}
void  LLR_to_XYZ(const double oldCoord[3], double newCoord[3])
{
	newCoord[ZDIM] = oldCoord[RDIM] * sin(oldCoord[LATDIM]);
	newCoord[YDIM] = oldCoord[RDIM] * cos(oldCoord[LATDIM])*sin(oldCoord[LONDIM]);
	newCoord[XDIM] = oldCoord[RDIM] * cos(oldCoord[LATDIM])*cos(oldCoord[LONDIM]);
}

void XYZ_to_LLR(const double oldCoord[3], double newCoord[3])
{
	newCoord[RDIM] = sqrt(oldCoord[XDIM] * oldCoord[XDIM] + oldCoord[YDIM] * oldCoord[YDIM] + oldCoord[ZDIM] * oldCoord[ZDIM]);
	newCoord[LATDIM] = atan2(oldCoord[ZDIM], sqrt(oldCoord[XDIM] * oldCoord[XDIM] + oldCoord[YDIM] * oldCoord[YDIM]));
	newCoord[LONDIM] = atan2(oldCoord[YDIM], oldCoord[XDIM]);
}
bool is_between_lon_range(double lower_boundary, double upper_boundary, double lon)
//is a longitude value between two given longitude values
{
	double lon_10 = upper_boundary - lower_boundary;
	if (lon_10 < 0) lon_10 += 2 * PI;
	if (lon_10 >= PI + PI) lon_10 = lon_10 - PI - PI;
	double dlon = lon - lower_boundary;
	if (dlon < 0) dlon += 2 * PI;
	if (dlon >= PI + PI) dlon = dlon - PI - PI;

	return (dlon <= lon_10);
}

bool is_a_valid_ray(double pt1[3], double pt2[3], double minCoords[3], double maxCoords[3], double* t_firstLast)
{
	/////////////////////////////////////////////////////////////////////////////
	//1、计算射线与顶底球面的交点及对应的t vector
	//2、快速剔除无效射线
	/////////////////////////////////////////////////////////////////////////////
	double tt_low[2], tt_high[2], lt[2], low_pt_t, high_pt_t, low_pt_xyz[3], high_pt_xyz[3], high_pt_llr[3], low_pt_llr[3];
	int num_t_low_sphere = t_vec_to_sphere(pt1, pt2, minCoords[RDIM], tt_low);
	int num_t_high_sphere = t_vec_to_sphere(pt1, pt2, maxCoords[RDIM], tt_high);


	//贯穿顶底球面
	if (num_t_low_sphere == 1 && num_t_high_sphere == 1) { low_pt_t = tt_low[0];	high_pt_t = tt_high[0]; }
	//由内部射向顶部球面，或由顶部射向内部
	else if (num_t_low_sphere == 0 && num_t_high_sphere == 1) {
		double r1 = sqrt(pt1[0] * pt1[0] + pt1[1] * pt1[1] + pt1[2] * pt1[2]);
		if (r1 > minCoords[RDIM] && r1 < maxCoords[RDIM])
			//p1位于大气层区域
		{
			low_pt_t = 0;
			high_pt_t = tt_high[0];
		}
		else//p2位于大气层区域
		{
			low_pt_t = 1;
			high_pt_t = tt_low[0];
		}
	}
	//由内部射向底部球面，或由底部射向内部 
	else if (num_t_low_sphere == 1 && num_t_high_sphere == 0) {
		double r1 = sqrt(pt1[0] * pt1[0] + pt1[1] * pt1[1] + pt1[2] * pt1[2]);
		if (r1 > minCoords[RDIM] && r1 < maxCoords[RDIM])
			//当p1位于大气层区域
		{
			low_pt_t = tt_low[0];
			high_pt_t = 0;
		}
		else//p2位于大气层区域
		{
			low_pt_t = tt_low[0];
			high_pt_t = 1;
		}
	}
	//水平穿过顶部球面
	else if (num_t_low_sphere == 0 && num_t_high_sphere == 2) { low_pt_t = tt_high[0];	high_pt_t = tt_high[1]; }
	//其余均为无效射线
	else	return false;

	XYZ_from_t_vec(pt1, pt2, low_pt_t, low_pt_xyz);
	XYZ_from_t_vec(pt1, pt2, high_pt_t, high_pt_xyz);
	XYZ_to_LLR(high_pt_xyz, high_pt_llr);
	XYZ_to_LLR(low_pt_xyz, low_pt_llr);

	//满足上述筛选后，还需如下进一步筛选
	//1)上述两交点落在研究区域在球面投影范围之外，也为无效射线
	if (!(is_between_lon_range(minCoords[LONDIM], maxCoords[LONDIM], high_pt_llr[LONDIM])
		&& is_between_lon_range(minCoords[LONDIM], maxCoords[LONDIM], low_pt_llr[LONDIM])
		&& high_pt_llr[LATDIM] >= minCoords[LATDIM] && high_pt_llr[LATDIM] <= maxCoords[LATDIM]
		&& low_pt_llr[LATDIM] >= minCoords[LATDIM] && low_pt_llr[LATDIM] <= maxCoords[LATDIM]
		)) return false;

	//2)只要射线与南北纬度面存在任何2个交点即无效射线，主要为了排除当射线水平（起始点纬度相同）穿过某一纬度面，因为纬度面是曲面，射线是直线
	if (t_vec_to_latitude_surface(pt1, pt2, minCoords[LATDIM], lt) == 2) return false;
	if (t_vec_to_latitude_surface(pt1, pt2, maxCoords[LATDIM], lt) == 2) return false;

	if (t_firstLast)
	{
		t_firstLast[0] = low_pt_t;
		t_firstLast[1] = high_pt_t;
	}
	return true;
	/////////////////////////////////////////////////////////////////////////////
}
string replaceStringAll(string str, const char* replace, const char* with) {
	//替换字符串
	if (replace) {
		size_t pos = 0;
		while ((pos = str.find(replace, pos)) != string::npos) {
			str.replace(pos, strlen(replace), with);
			pos += strlen(with);
		}
	}
	return str;
}
string trim(const char* str)
//去两边空格
{
	string s = str;
	if (s.empty())
	{
		return s;
	}
	s.erase(0, s.find_first_not_of(" "));
	s.erase(s.find_last_not_of(" ") + 1);
	return s;
}
void get_list_from_filter(const string& filterStr, int defaultNum, vector<string>& list)
{
	string filter = trim(filterStr.c_str());
	regex sep("[{}]*[{}]");
	sregex_token_iterator p(filter.cbegin(), filter.cend(), // sequence
		sep, // separator
		{ -1 }); // -1: values between separators
	sregex_token_iterator e;
	vector<string> subNums;
	string outfilename;
	for (; p != e; ++p) {
		string str = *p;
		subNums.push_back(str);
	}
	int num = subNums.size() / 2;
	vector<string>* filenames = new vector<string>[2];
	int oldIndex = 0, newIndex = 1;
	filenames[0].push_back(subNums[0]);
	for (int i = 0; i < num; i++)
	{
		int min, max, step, nDigital;
		subNums[2 * i + 1] = replaceStringAll(subNums[2 * i + 1], ",", " ");
		stringstream ss(subNums[2 * i + 1]);
		ss >> min >> step >> max >> nDigital;
		if (max == -1) max = (defaultNum + min - 1) * step;
		for (int iOld = 0; iOld < filenames[oldIndex].size(); iOld++)
			for (int j = min; j <= max; j += step)
			{
				string fmt = filenames[oldIndex][iOld] + "%0*d";
				char ch[1024];
				sprintf_s(ch,1024, fmt.c_str(), nDigital, j);
				string tmp = 2 * (i + 1) >= subNums.size() ? string(ch) : string(ch) + subNums[2 * (i + 1)];
				filenames[newIndex].push_back(tmp);
				//filenames[newIndex].push_back(string(ch) + subNums[2 * (i + 1)]);
			}
		filenames[oldIndex].clear();
		int tmp = oldIndex;
		oldIndex = newIndex;
		newIndex = tmp;
	}
	if (num != 0)
	{
		for (int i = 0; i < filenames[oldIndex].size(); i++) {
			outfilename = filenames[oldIndex][i].c_str();
			list.push_back(outfilename);
		}
	}
	else
	{
		for (int i = 0; i < defaultNum; i++)
		{
			list.push_back(filter);
		}
	}
	delete[] filenames;
}