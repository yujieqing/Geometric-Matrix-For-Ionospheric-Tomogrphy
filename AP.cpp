#include "AP.h"
#include <assert.h>
#include <algorithm>
#include <vector>
bool cmpfunction1(double& a, double & b) {
	return (a < b);
}
bool LLR_box_of_ray_plane(const double startPt[3], const double endPt[3], double* split_surf_coords[3], const int num_of_split_surfs[3],double minLLR[3],double maxLLR[3])
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function: get the bounding box of the ray plane which is formed by startPt, endPt and the geocenter in longitude and latitude coordinates
//Input:
//		startPt: first point
//      endPt: second point
//      split_surf_coords: surfaces that split the voxel model
//		num_of_split_surfs: number of the surfaces in each dimension
// output:
//		minLLR: the minmal LLR coordinate of the bounding box
//		maxLLR: the maxmal LLR coordinate of the bounding box
//return:
//		true: success
//		false: failed
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

{
	double t_vector[2], cur_xyz[3],llr0[3],llr1[3],cur_t;

	int num = t_vec_to_sphere(startPt, endPt, split_surf_coords[RDIM][0], t_vector);
	if (num < 1) return false;
	cur_t = t_vector[0];
	if (num == 2 && t_vector[1]>=0)		cur_t = t_vector[1];
	if (cur_t < 0 || cur_t>1) return false;
	XYZ_from_t_vec(startPt, endPt, cur_t, cur_xyz);
	XYZ_to_LLR(cur_xyz, llr0);
	num = t_vec_to_sphere(startPt, endPt, split_surf_coords[RDIM][num_of_split_surfs[RDIM] - 1], t_vector);
	if (num < 1) return false;
	cur_t = t_vector[0];
	if (num == 2 && t_vector[1] >= 0)		cur_t = t_vector[1];
	if (cur_t < 0 || cur_t>1) return false;
	XYZ_from_t_vec(startPt, endPt, cur_t, cur_xyz);
	XYZ_to_LLR(cur_xyz, llr1);
	minLLR[LONDIM] = llr0[LONDIM] < llr1[LONDIM] ? llr0[LONDIM] : llr1[LONDIM];
	minLLR[LATDIM] = llr0[LATDIM] < llr1[LATDIM] ? llr0[LATDIM] : llr1[LATDIM];
	minLLR[RDIM] = split_surf_coords[RDIM][0];
	maxLLR[LONDIM] = llr0[LONDIM] > llr1[LONDIM] ? llr0[LONDIM] : llr1[LONDIM];
	maxLLR[LATDIM] = llr0[LATDIM] > llr1[LATDIM] ? llr0[LATDIM] : llr1[LATDIM];
	maxLLR[RDIM] = split_surf_coords[RDIM][num_of_split_surfs[RDIM] - 1];
	return true;

}
bool isIntersectedWithPlane(const double coeff[4], const double voxelMinLLR[3], const double voxelMaxLLR[3])
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function: test whether a voxel is interseced with the plane given by ax+by+cz+d=0
//Input:
//		coeff: the coeffcients of the plane equation, a,b,c,d
//      voxelMinLLR: the minmal LLR coordinates of the voxel to be tested
//      voxelMaxLLR: the maxmal LLR coordinates of the voxel to be tested
// output:
//		minLLR: the minimal LLR coordinate of the bounding box
//		maxLLR: the maximal LLR coordinate of the bounding box
//return:
//		true: intersected
//		false: not intersected
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	//the eight corners of the voxel in LLR coordinates
	double voxCorners[8][3] = {
		voxelMinLLR[0],voxelMinLLR[1],voxelMinLLR[2],
		voxelMinLLR[0],voxelMinLLR[1],voxelMaxLLR[2],
		voxelMinLLR[0],voxelMaxLLR[1],voxelMinLLR[2],
		voxelMinLLR[0],voxelMaxLLR[1],voxelMaxLLR[2],
		voxelMaxLLR[0],voxelMinLLR[1],voxelMinLLR[2],
		voxelMaxLLR[0],voxelMinLLR[1],voxelMaxLLR[2],
		voxelMaxLLR[0],voxelMaxLLR[1],voxelMinLLR[2],
		voxelMaxLLR[0],voxelMaxLLR[1],voxelMaxLLR[2]
	};
	int negtive = 0, positive = 0;
	double xyz[3];
	double xyzMin[3] = { 1e20,1e20,1e20 };
	double xyzMax[3]= { -1e20,-1e20,-1e20 };
	for (int i = 0; i < 8; i++)
	{
		LLR_to_XYZ(voxCorners[i],xyz);
		for (int j = 0; j < 3; j++)
		{
			xyzMin[j] = xyzMin[j] < xyz[j] ? xyzMin[j] : xyz[j];
			xyzMax[j] = xyzMax[j] > xyz[j] ? xyzMax[j] : xyz[j];
		}
	}
	double dx, dy,dz;
	dx = (xyzMax[XDIM] - xyzMin[XDIM]) / 10;
	dy= (xyzMax[YDIM] - xyzMin[YDIM]) / 10;
	dz = (xyzMax[ZDIM] - xyzMin[ZDIM]) / 10;
	//expand the voxel for 1/10 in case of wrong judge caused by floating operation 
	xyzMin[XDIM] -= dx; xyzMin[YDIM] -= dy; xyzMin[ZDIM] -= dz;
	xyzMax[XDIM] += dx; xyzMax[YDIM] += dy; xyzMax[ZDIM] += dz;

	//the eight corners of the voxel in xyz coordinates
	double voxCornersXYZ[8][3] = {
		xyzMin[0],xyzMin[1],xyzMin[2],
		xyzMin[0],xyzMin[1],xyzMax[2],
		xyzMin[0],xyzMax[1],xyzMin[2],
		xyzMin[0],xyzMax[1],xyzMax[2],
		xyzMax[0],xyzMin[1],xyzMin[2],
		xyzMax[0],xyzMin[1],xyzMax[2],
		xyzMax[0],xyzMax[1],xyzMin[2],
		xyzMax[0],xyzMax[1],xyzMax[2]
	};
	// test whether each corner is above (>=0) or below (<0) the plane
	for (int i = 0; i < 8; i++)
	{
		if (coeff[0]*voxCornersXYZ[i][0] + coeff[1] * voxCornersXYZ[i][1]+ coeff[2] * voxCornersXYZ[i][2]+ coeff[3] >= 0)
			positive++;
		else
			negtive++;
	}
	
	if (negtive == 8 || positive == 8) return false;// all eight corners are on the same side
	else return true;// one of corner are on the opposite side of the other corner
}
bool isCollinear(const double p1[3], const double p2[3], const double p3[3])
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function: test if three points are collinear
//Input:
//		p1: first point
//      p3: second point
//      p3: third point
//return:
//		true: collinear
//		false: non-collinear
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	//if A, B, C are collinear, then there is a ¦Ë that make AB = ¦Ë*AC
	//(x2 - x1, y2 - y1, z2 - z1) = ¦Ë(x3 - x1, y3 - y1, z3 - z1)
	//x2 - x1 = ¦Ë(x3 - x1), y2 - y1 = ¦Ë(y3 - y1), z2 - z1 = ¦Ë(z3 - z1).
	double err0 = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]);
	double err1 = (p2[2] - p1[2]) * (p3[1] - p1[1]) - (p3[2] - p1[2]) * (p2[1] - p1[1]);
	if (err0 > 1e-11 || err0 < -1e-11 || err1>1e-11 || err1 < -1e-11)
		return false;
	else
		return true;
}
bool planeCoefficent(const double B[3], const double C[3],double coeff[4])
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function: calucate the general equation(ax+by+cz+d=0) of the plane formed by the center of the earth and the two end-points of the ray
//https://au.mathworks.com/matlabcentral/fileexchange/97692-equation-of-a-plane-passing-through-three-points
//Input:
//		B: start point of a ray
//      C: end point of a ray
//Output:
//      coeff: a£¬b,c£¬d 
//return:
//		true: the 3 Points are Non-Collinear
//		false: the 3 points are Collinear
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	double A[3] = {0,0,0};
	if (isCollinear(B, C, A)) return false;
	coeff[0] = (B[2] - A[2]) * (C[3] - A[3]) - (C[2] - A[2]) * (B[3] - A[3]);
	coeff[1] = (B[3] - A[3]) * (C[1] - A[1]) - (C[3] - A[3]) * (B[1] - A[1]);
	coeff[2] = (B[1] - A[1]) * (C[2] - A[2]) - (C[1] - A[1]) * (B[2] - A[2]);
	coeff[3] = -(coeff[0] * A[1] + coeff[1] * A[2] + coeff[2] * A[3]);
	return true;
}
double segment_length_in_a_voxel(double minCoords[3], double maxCoords[3], double startPt[3], double endPt[3], double vectNormal)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function: calucate the length of a ray inside a voxel
//Input:
//		minCoords: the minimum value of each dimension for a voxel
//      maxCoords: the maximum value of each dimension for a voxel
//		startPt: the start point of a ray
//      endPt: the end point of a ray
//		vectNormal: the length of the ray's vector
//Output:
//      orientVector: the length of the ray inside the voxel

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	double t_vals[10];
	int nCount = 0;

	double tmp_t[2], tmp_XYZ[3], tmp_LLR[3];

	int tmp_num;
	// intersetion with the first longitudinal surface 
	if (t_vec_to_longitude_surface(startPt, endPt, minCoords[LONDIM], tmp_t[0]) > 0)
	{
		XYZ_from_t_vec(startPt, endPt, tmp_t[0], tmp_XYZ);
		XYZ_to_LLR(tmp_XYZ, tmp_LLR);
		if (tmp_LLR[LATDIM] >= minCoords[LATDIM] && tmp_LLR[LATDIM] <= maxCoords[LATDIM]
			&& tmp_LLR[RDIM] >= minCoords[RDIM] && tmp_LLR[RDIM] <= maxCoords[RDIM])
		{
			t_vals[nCount] = tmp_t[0];
			nCount++;
		}
	}
	// intersetion with the second longitudinal surface 
	if (t_vec_to_longitude_surface(startPt, endPt, maxCoords[LONDIM], tmp_t[0]) > 0)
	{
		XYZ_from_t_vec(startPt, endPt, tmp_t[0], tmp_XYZ);
		XYZ_to_LLR(tmp_XYZ, tmp_LLR);
		if (tmp_LLR[LATDIM] >= minCoords[LATDIM] && tmp_LLR[LATDIM] <= maxCoords[LATDIM]
			&& tmp_LLR[RDIM] >= minCoords[RDIM] && tmp_LLR[RDIM] <= maxCoords[RDIM])
		{
			t_vals[nCount] = tmp_t[0];
			nCount++;
		}
	}

	// intersetion with the first latitudinal surface 
	tmp_num = t_vec_to_latitude_surface(startPt, endPt, minCoords[LATDIM], tmp_t);
	for (int i = 0; i < tmp_num; i++)
	{
		XYZ_from_t_vec(startPt, endPt, tmp_t[i], tmp_XYZ);
		XYZ_to_LLR(tmp_XYZ, tmp_LLR);
		if (tmp_LLR[RDIM] >= minCoords[RDIM] && tmp_LLR[RDIM] <= maxCoords[RDIM]
			&& is_between_lon_range(minCoords[LONDIM], maxCoords[LONDIM], tmp_LLR[LONDIM]))
		{
			t_vals[nCount] = tmp_t[i];
			nCount++;
		}
	}

	// intersetion with the second latitudinal surface 
	tmp_num = t_vec_to_latitude_surface(startPt, endPt, maxCoords[LATDIM], tmp_t);
	for (int i = 0; i < tmp_num; i++)
	{
		XYZ_from_t_vec(startPt, endPt, tmp_t[i], tmp_XYZ);
		XYZ_to_LLR(tmp_XYZ, tmp_LLR);
		if (tmp_LLR[RDIM] >= minCoords[RDIM] && tmp_LLR[RDIM] <= maxCoords[RDIM]
			&& is_between_lon_range(minCoords[LONDIM], maxCoords[LONDIM], tmp_LLR[LONDIM]))
		{
			t_vals[nCount] = tmp_t[i];
			nCount++;
		}
	}

	// intersetion with the first spherical surface 
	tmp_num = t_vec_to_sphere(startPt, endPt, minCoords[RDIM], tmp_t);
	for (int i = 0; i < tmp_num; i++)
	{
		XYZ_from_t_vec(startPt, endPt, tmp_t[i], tmp_XYZ);
		XYZ_to_LLR(tmp_XYZ, tmp_LLR);
		if (tmp_LLR[LATDIM] >= minCoords[LATDIM] && tmp_LLR[LATDIM] <= maxCoords[LATDIM]
			&& is_between_lon_range(minCoords[LONDIM], maxCoords[LONDIM], tmp_LLR[LONDIM]))
		{
			t_vals[nCount] = tmp_t[i];
			nCount++;
		}
	}

	// intersetion with the second spherical surface 
	tmp_num = t_vec_to_sphere(startPt, endPt, maxCoords[RDIM], tmp_t);
	for (int i = 0; i < tmp_num; i++)
	{
		XYZ_from_t_vec(startPt, endPt, tmp_t[i], tmp_XYZ);
		XYZ_to_LLR(tmp_XYZ, tmp_LLR);
		if (tmp_LLR[LATDIM] >= minCoords[LATDIM] && tmp_LLR[LATDIM] <= maxCoords[LATDIM]
			&& is_between_lon_range(minCoords[LONDIM], maxCoords[LONDIM], tmp_LLR[LONDIM]))
		{
			t_vals[nCount] = tmp_t[i];
			nCount++;
		}
	}
	if (nCount < 2) return 0;
	if (nCount > 6) return 0;//there are six possible intersections maximally. Under this case, the two intersections are two vertexs of the intersecting voxel themself
	if (nCount == 2) 
		return fabs(t_vals[1] - t_vals[0])*vectNormal;
	else //merge those intersections that are very close to each other
	{
		std::sort(t_vals, t_vals + nCount, cmpfunction1);
		double t_finals[6];
		int iCurr_final = 0;
		t_finals[iCurr_final] = t_vals[0];
		for (int i = 1; i < nCount; i++)
		{
			if (fabs(t_finals[iCurr_final] - t_vals[i]) > THRESHOLD_T_OF_TWO_POINTS)
			{
				iCurr_final++;
				t_finals[iCurr_final] = t_vals[i];
			}
			i++;
		}
		assert(iCurr_final > 1);// If this is happened, we should increase the THRESHOLD_T_OF_TWO_POINTS
		return fabs(t_finals[1] - t_finals[0])*vectNormal;
	}
}
int segments_of_a_ray_by_tradition(double startPt[3], double endPt[3], double* split_surf_coords[3], int num_of_split_surfs[3], std::vector<INDEX_LENGTH*>& segments)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function: calucate the segments' length using tradition method, that is, test whether the ray is interseting with every voxel. If intersected, compute the length of the ray inside the voxel.
//Input:
//		startPt: start point (x,y,z)of a ray
//      endPt: end point (x,y,z)of a ray
//      split_surf_coords: an two dimension array for the splitting surfaces' coordinates 
//      num_of_split_surfs: the number of splitting surfaces along each dimension
//Output:
//      segments: the segments' length using tradition method

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	double minCoords[3] = { split_surf_coords[LONDIM][0],split_surf_coords[LATDIM][0],split_surf_coords[RDIM][0] };
	double maxCoords[3] = { split_surf_coords[LONDIM][num_of_split_surfs[LONDIM] - 1],split_surf_coords[LATDIM][num_of_split_surfs[LATDIM] - 1],split_surf_coords[RDIM][num_of_split_surfs[RDIM] - 1] };
	if (!is_a_valid_ray(startPt, endPt, minCoords, maxCoords, NULL)) return 0;


	double normal = sqrt((startPt[XDIM] - endPt[XDIM])*(startPt[XDIM] - endPt[XDIM]) + (startPt[YDIM] - endPt[YDIM])
		*(startPt[YDIM] - endPt[YDIM]) + (startPt[ZDIM] - endPt[ZDIM])*(startPt[ZDIM] - endPt[ZDIM]));
	double len = 0;
	long nCount = 0;
	long dim10 = (num_of_split_surfs[LATDIM] - 1) * (num_of_split_surfs[LONDIM] - 1);
	double plCoeffs[4];
	if (!planeCoefficent(startPt, endPt, plCoeffs))
		//the ray is a vertical ray, i.e.,  it is consisit of three collinear points 
	{
		double startPtLLR[3];
		XYZ_to_LLR(startPt, startPtLLR);
		nCount = num_of_split_surfs[RDIM] - 1;
		for (int k = 0; k < nCount; k++)
		{
			INDEX_LENGTH* pnew = new INDEX_LENGTH;
			pnew->intercept = split_surf_coords[RDIM][k + 1] - split_surf_coords[RDIM][k];
			pnew->index[LONDIM] = lower_bound_index(split_surf_coords[LONDIM], num_of_split_surfs[LONDIM], startPtLLR[LONDIM]);
			pnew->index[LATDIM] = lower_bound_index(split_surf_coords[LATDIM], num_of_split_surfs[LATDIM], startPtLLR[LATDIM]);
			pnew->index[RDIM] = k;
			pnew->id = pnew->index[RDIM] * dim10 + pnew->index[LATDIM] * (num_of_split_surfs[LONDIM] - 1) + pnew->index[LONDIM];
			segments.push_back(pnew);
		}
		return nCount;
	}

	for (int k = 0; k < num_of_split_surfs[RDIM] - 1; k++)
	{
		minCoords[RDIM] = split_surf_coords[RDIM][k];
		maxCoords[RDIM] = split_surf_coords[RDIM][k + 1];
		for (int j = 0; j < num_of_split_surfs[LATDIM] - 1; j++)
		{
			minCoords[LATDIM] = split_surf_coords[LATDIM][j];
			maxCoords[LATDIM] = split_surf_coords[LATDIM][j + 1];
			for (int i = 0; i < num_of_split_surfs[LONDIM] - 1; i++)
			{
				minCoords[LONDIM] = split_surf_coords[LONDIM][i];
				maxCoords[LONDIM] = split_surf_coords[LONDIM][i + 1];
				if (!isIntersectedWithPlane(plCoeffs, minCoords, maxCoords)) continue;

				if ((len = segment_length_in_a_voxel(minCoords, maxCoords, startPt, endPt, normal)) > 0)
				{
					INDEX_LENGTH* pnew = new INDEX_LENGTH;
					pnew->intercept = len;
					pnew->index[LONDIM] = i;
					pnew->index[LATDIM] = j;
					pnew->index[RDIM] = k;
					pnew->id = pnew->index[RDIM] * dim10 + pnew->index[LATDIM] * (num_of_split_surfs[LONDIM] - 1) + pnew->index[LONDIM];
					segments.push_back(pnew);
					nCount++;
				}
			}
		}
	}
	return nCount;
}
int segments_of_a_ray_by_AP(double startPt[3], double endPt[3], double* split_surf_coords[3], int num_of_split_surfs[3], std::vector<INDEX_LENGTH*>& segments)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function: calucate the segments' length using Hong et al. 2017, Tomography Reconstruction of Ionospheric Electron Density with Empirical Orthonormal Functions Using Korea GNSS Network method, J. Astron. Space Sci. 34(1), 7-17 (2017)
//Input:
//		startPt: start point of a ray
//      endPt: end point of a ray
//      split_surf_coords: an two dimension array for the splitting surfaces' coordinates 
//      num_of_split_surfs: the number of splitting surfaces along each dimension
//Output:
//      segments: the segments' length using tradition method

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	double minCoords[3] = { split_surf_coords[LONDIM][0],split_surf_coords[LATDIM][0],split_surf_coords[RDIM][0] };
	double maxCoords[3] = { split_surf_coords[LONDIM][num_of_split_surfs[LONDIM] - 1],split_surf_coords[LATDIM][num_of_split_surfs[LATDIM] - 1],split_surf_coords[RDIM][num_of_split_surfs[RDIM] - 1] };
	if (!is_a_valid_ray(startPt, endPt, minCoords, maxCoords, NULL)) return 0;// test whether it is a valid ray

	double plCoeffs[4];
	long nCount = 0;
	long dim10 = (num_of_split_surfs[LATDIM] - 1) * (num_of_split_surfs[LONDIM] - 1);
	double geocenter[3] = { 0,0,0 };
	if(isCollinear(startPt, endPt, geocenter))
	//if (!planeCoefficent(startPt, endPt, plCoeffs))
		//the ray is a vertical ray, i.e.,  it is consisit of three collinear points 
	{
		double startPtLLR[3];
		XYZ_to_LLR(startPt, startPtLLR);
		nCount = num_of_split_surfs[RDIM] - 1;
		for (int k = 0; k < nCount; k++)
		{
			INDEX_LENGTH* pnew = new INDEX_LENGTH;
			pnew->intercept = split_surf_coords[RDIM][k+1]- split_surf_coords[RDIM][k];
			pnew->index[LONDIM] = lower_bound_index(split_surf_coords[LONDIM], num_of_split_surfs[LONDIM], startPtLLR[LONDIM]);
			pnew->index[LATDIM] = lower_bound_index(split_surf_coords[LATDIM], num_of_split_surfs[LATDIM], startPtLLR[LATDIM]);
			pnew->index[RDIM] = k;
			pnew->id = pnew->index[RDIM] * dim10 + pnew->index[LATDIM] * (num_of_split_surfs[LONDIM] - 1) + pnew->index[LONDIM];
			segments.push_back(pnew);
		}
		return nCount;
	}

	// the three points are non-collinear

	double planeboxMinLLR[3], planeboxMaxLLR[3];
	int minIndex[3], maxIndex[3];
	//get the bounding box of the plane in longitude and latitude coordinates
	if(!LLR_box_of_ray_plane(startPt, endPt, split_surf_coords, num_of_split_surfs, planeboxMinLLR, planeboxMaxLLR)) return 0;

	double normal = sqrt((startPt[XDIM] - endPt[XDIM]) * (startPt[XDIM] - endPt[XDIM]) + (startPt[YDIM] - endPt[YDIM])
		* (startPt[YDIM] - endPt[YDIM]) + (startPt[ZDIM] - endPt[ZDIM]) * (startPt[ZDIM] - endPt[ZDIM]));
	double len = 0;
	//get the index range of the bounding box 
	minIndex[LONDIM]=lower_bound_index(split_surf_coords[LONDIM], num_of_split_surfs[LONDIM], planeboxMinLLR[LONDIM]);
	maxIndex[LONDIM] = upper_bound_index(split_surf_coords[LONDIM], num_of_split_surfs[LONDIM], planeboxMaxLLR[LONDIM]);
	minIndex[LATDIM] = lower_bound_index(split_surf_coords[LATDIM], num_of_split_surfs[LATDIM], planeboxMinLLR[LATDIM]);
	maxIndex[LATDIM] = upper_bound_index(split_surf_coords[LATDIM], num_of_split_surfs[LATDIM], planeboxMaxLLR[LATDIM]);

	for (int k = 0; k < num_of_split_surfs[RDIM] - 1; k++)
	{
		minCoords[RDIM] = split_surf_coords[RDIM][k];
		maxCoords[RDIM] = split_surf_coords[RDIM][k + 1];
		for (int j = minIndex[LATDIM]; j <= maxIndex[LATDIM]; j++)
		{
			minCoords[LATDIM] = split_surf_coords[LATDIM][j];
			maxCoords[LATDIM] = split_surf_coords[LATDIM][j + 1];
			for (int i = minIndex[LONDIM]; i <= maxIndex[LONDIM]; i++)
			{
				minCoords[LONDIM] = split_surf_coords[LONDIM][i];
				maxCoords[LONDIM] = split_surf_coords[LONDIM][i + 1];

				//if (!isIntersectedWithPlane(plCoeffs, minCoords, maxCoords)) continue;
				if ((len = segment_length_in_a_voxel(minCoords, maxCoords, startPt, endPt, normal)) > 0)
				{
					INDEX_LENGTH* pnew = new INDEX_LENGTH;
					pnew->intercept = len;
					pnew->index[LONDIM] = i;
					pnew->index[LATDIM] = j;
					pnew->index[RDIM] = k;
					pnew->id = pnew->index[RDIM] * dim10 + pnew->index[LATDIM] * (num_of_split_surfs[LONDIM] - 1) + pnew->index[LONDIM];
					segments.push_back(pnew);
					nCount++;
				}
			}
		}
	}
	return nCount;
}
