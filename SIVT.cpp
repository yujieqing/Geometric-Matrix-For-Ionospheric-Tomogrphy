#include "SIVT.h"
#include "common.h"
#include <math.h>
#include <assert.h>
#include <algorithm>


bool cmpfunction0(NEXTVISIT& a, NEXTVISIT & b) {
	return (a.t_vector < b.t_vector);
}

void orientation_vector_of_a_ray(double start_pt_xyz[3], double end_pt_xyz[3], int orientVector[3], bool twoVal[2], double t_tangent[2])
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function: calucate the orientation vector for a given ray and the corresponding tangent point
//Input:
//		start_pt_xyz: the start point of a ray
//      end_pt_xyz: the end point of a ray
//Output:
//      orientVector: the orientation vector of the ray. Notice that the second value for latitudinal and radial dimension is always opposite to the first value, so we need to record only value for each dimension
//      twoVal: indicate that if there are two values for the orientation vector in the corresponding dimension. [0] is for latitudinal dimension, and [1] is for radial dimension
//      t_tangent: the tangent point along latitudinal and radial dimension//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	double end_pt_llr[3], start_pt_llr[3];
	XYZ_to_LLR(end_pt_xyz, end_pt_llr);
	XYZ_to_LLR(start_pt_xyz, start_pt_llr);

	orientVector[LONDIM] = 0;
	//calculate the vector from the start point to the end point of a ray
	double vec_start_to_end_pt[3] = { end_pt_xyz[XDIM] - start_pt_xyz[XDIM],end_pt_xyz[YDIM] - start_pt_xyz[YDIM], end_pt_xyz[ZDIM] - start_pt_xyz[ZDIM] };

	double vec_geocenter_2_start_pt[3] = { start_pt_xyz[XDIM], start_pt_xyz[YDIM],start_pt_xyz[ZDIM] };// vector from geocenter to the start point
	double cross_prod = vec_geocenter_2_start_pt[XDIM] * vec_start_to_end_pt[YDIM] - vec_geocenter_2_start_pt[YDIM] * vec_start_to_end_pt[XDIM];
	// if the cross product of the two vector greater than 0, the orientation of longitudinal dimension is positive, otherwise, negative 
	if (cross_prod == 0)	orientVector[LONDIM] = 0;
	else if (cross_prod > 0) orientVector[LONDIM] = 1;
	else orientVector[LONDIM] = -1;

	double ttmp = vec_start_to_end_pt[XDIM] * start_pt_xyz[XDIM] + vec_start_to_end_pt[YDIM] * start_pt_xyz[YDIM] + vec_start_to_end_pt[ZDIM] * start_pt_xyz[ZDIM];
	t_tangent[0] = start_pt_xyz[ZDIM] * ttmp;
	t_tangent[0] -= vec_start_to_end_pt[ZDIM] * (start_pt_xyz[XDIM] * start_pt_xyz[XDIM] + start_pt_xyz[YDIM] * start_pt_xyz[YDIM] + start_pt_xyz[ZDIM] * start_pt_xyz[ZDIM]);
	double denominator = vec_start_to_end_pt[ZDIM] * ttmp - start_pt_xyz[ZDIM] * (vec_start_to_end_pt[XDIM] * vec_start_to_end_pt[XDIM] + vec_start_to_end_pt[YDIM] * vec_start_to_end_pt[YDIM] + vec_start_to_end_pt[ZDIM] * vec_start_to_end_pt[ZDIM]);
	twoVal[0] = false;
	orientVector[LATDIM] = end_pt_llr[LATDIM] > start_pt_llr[LATDIM] ? 1 : -1;// normal case
	if (denominator != 0)//denominator should not be equal to 0
	{
		t_tangent[0] /= denominator; // the t-value of the tangent point in latitudinal dimension
		if (t_tangent[0] > 0 && t_tangent[0] < 1) // the exceptional case for latitudinal orientation. there are two values for the vector in this dimension
		{
			twoVal[0] = true;
			double tangent_pt_xyz[3], tangent_pt_llr[3];
			XYZ_from_t_vec(start_pt_xyz, end_pt_xyz, t_tangent[0], tangent_pt_xyz);
			XYZ_to_LLR(tangent_pt_xyz, tangent_pt_llr);
			orientVector[LATDIM] = tangent_pt_llr[LATDIM] > start_pt_llr[LATDIM] ? 1 : -1;
		}
	}


	// in this implementation, the ray is considered to be an invalid ray, if it has two valid intersections with the same sphere. Under this condition, there are only one value for the radial vector
	t_tangent[1] = 0;
	twoVal[1] = false;
	orientVector[RDIM] = 1;
	double vec_start_pt_2_geocenter[3] = { -start_pt_xyz[XDIM],-start_pt_xyz[YDIM],-start_pt_xyz[ZDIM] };
	double dot_prod = vec_start_pt_2_geocenter[XDIM] * vec_start_to_end_pt[XDIM] + vec_start_pt_2_geocenter[YDIM] * vec_start_to_end_pt[YDIM] + vec_start_pt_2_geocenter[ZDIM] * vec_start_to_end_pt[ZDIM];
	if (dot_prod == 0) { orientVector[RDIM] = 0; assert(dot_prod == 0); }
	else
		orientVector[RDIM] = dot_prod < 0 ? 1 : -1;

}


int t_vector_to_index_and_length(NEXTVISIT* t_vec_each_surf, int num_of_t_vector, double len_of_vector, const INDEX_LENGTH& first_voxel, int num_of_coords[3], INDEX_LENGTH*& segment_length)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function: calculate the voxel's index and the cooresponding segment length
//Input:
//		t_vec_each_surf: the t-values of the intersections on each splitting surface
//		num_of_t_vector: the size of the array t_vec_each_surf
//      len_of_vector: the length of the vector for a ray
//		first_voxel: the first intersecting voxel
//		num_of_coords: the number of splitting surface for each dimnesion
//Output:
//      segment_length: the voxel's index and the cooresponding segment length
//		return the number of intersecting voxels
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	long dim10 = (num_of_coords[LATDIM] - 1) * (num_of_coords[LONDIM] - 1);
	if (segment_length) delete[] segment_length;
	segment_length = new INDEX_LENGTH[num_of_t_vector + 1];
	segment_length[0].index[0] = first_voxel.index[0];
	segment_length[0].index[1] = first_voxel.index[1];
	segment_length[0].index[2] = first_voxel.index[2];
	segment_length[0].t_vec_first_pt = first_voxel.t_vec_first_pt;
	segment_length[0].id = segment_length[0].index[RDIM] * dim10 + segment_length[0].index[LATDIM] * (num_of_coords[LONDIM] - 1) + segment_length[0].index[LONDIM];


	int vox_index = 1;//starting from the first intersecting voxel
	int i = 0;//starting from the first t-value
	while (i < num_of_t_vector)
	{

		if ((i + 1) < num_of_t_vector && (t_vec_each_surf[i + 1].t_vector - t_vec_each_surf[i].t_vector) < THRESHOLD_T_OF_TWO_POINTS)
			//the ray is crossing two or three splitting surfaces at the same time
		{
			if ((i + 2) < num_of_t_vector && (t_vec_each_surf[i + 2].t_vector - t_vec_each_surf[i + 1].t_vector) < THRESHOLD_T_OF_TWO_POINTS)
				//crossing three splitting surfaces at the same time
			{
				segment_length[vox_index].t_vec_first_pt = (t_vec_each_surf[i].t_vector + t_vec_each_surf[i + 1].t_vector + t_vec_each_surf[i + 2].t_vector) / 3;
				//segment_length[vox_index-1].intercept =(t - segment_length[vox_index - 1].t_vec_first_pt) * len_of_vector;
				segment_length[vox_index].index[0] = segment_length[vox_index - 1].index[LONDIM] + t_vec_each_surf[i].next_visit[LONDIM] + t_vec_each_surf[i + 1].next_visit[LONDIM] + t_vec_each_surf[i + 2].next_visit[LONDIM];
				segment_length[vox_index].index[1] = segment_length[vox_index - 1].index[LATDIM] + t_vec_each_surf[i].next_visit[LATDIM] + t_vec_each_surf[i + 1].next_visit[LATDIM] + t_vec_each_surf[i + 2].next_visit[LATDIM];
				segment_length[vox_index].index[2] = segment_length[vox_index - 1].index[RDIM] + t_vec_each_surf[i].next_visit[RDIM] + t_vec_each_surf[i + 1].next_visit[RDIM] + t_vec_each_surf[i + 2].next_visit[RDIM];
				segment_length[vox_index].id = segment_length[vox_index].index[RDIM] * dim10 + segment_length[vox_index].index[LATDIM] * (num_of_coords[LONDIM] - 1) + segment_length[vox_index].index[LONDIM];
				vox_index++;
				i = i + 3;
			}
			else//crossing two splitting surfaces at the same time
			{
				segment_length[vox_index].t_vec_first_pt = (t_vec_each_surf[i].t_vector + t_vec_each_surf[i + 1].t_vector) / 2;
				//segment_length[vox_index].intercept = (t - segment_length[vox_index - 1].t_vec_first_pt)*len_of_vector;
				segment_length[vox_index].index[0] = segment_length[vox_index - 1].index[LONDIM] + t_vec_each_surf[i].next_visit[LONDIM] + t_vec_each_surf[i + 1].next_visit[LONDIM];
				segment_length[vox_index].index[1] = segment_length[vox_index - 1].index[LATDIM] + t_vec_each_surf[i].next_visit[LATDIM] + t_vec_each_surf[i + 1].next_visit[LATDIM];
				segment_length[vox_index].index[2] = segment_length[vox_index - 1].index[RDIM] + t_vec_each_surf[i].next_visit[RDIM] + t_vec_each_surf[i + 1].next_visit[RDIM];
				segment_length[vox_index].id = segment_length[vox_index].index[RDIM] * dim10 + segment_length[vox_index].index[LATDIM] * (num_of_coords[LONDIM] - 1) + segment_length[vox_index].index[LONDIM];

				vox_index++;
				i = i + 2;
			}
		}
		else//crossing only one splitting surfaces
		{
			segment_length[vox_index].t_vec_first_pt = t_vec_each_surf[i].t_vector;
			//segment_length[vox_index].intercept = (t - segment_length[vox_index - 1].t_vec_first_pt)*len_of_vector;
			segment_length[vox_index].index[0] = segment_length[vox_index - 1].index[LONDIM] + t_vec_each_surf[i].next_visit[LONDIM];
			segment_length[vox_index].index[1] = segment_length[vox_index - 1].index[LATDIM] + t_vec_each_surf[i].next_visit[LATDIM];
			segment_length[vox_index].index[2] = segment_length[vox_index - 1].index[RDIM] + t_vec_each_surf[i].next_visit[RDIM];
			segment_length[vox_index].id = segment_length[vox_index].index[RDIM] * dim10 + segment_length[vox_index].index[LATDIM] * (num_of_coords[LONDIM] - 1) + segment_length[vox_index].index[LONDIM];

			vox_index++;
			i = i + 1;
		}

	}
	int num_of_voxel = vox_index - 1;	//the last voxel is outside the range of the voxel model, and should be removed

	//calculte the absolute length for each segment according to two consecutive t-values
	for (int i = 0; i < num_of_voxel; i++)
	{
		segment_length[i].intercept = (segment_length[i + 1].t_vec_first_pt - segment_length[i].t_vec_first_pt)*len_of_vector;
	}

	return num_of_voxel;
}
int segments_of_a_ray_by_AP(double start_pt[3], double end_pt[3], double* split_surf_coords[3], int num_of_split_surfs[3], INDEX_LENGTH*& segment_length)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function: obtain the segments' lengths and owners of a ray and a voxel model
//Input:
//		start_pt: the start point of a ray
//		end_pt: the end point of a ray
//      split_surf_coords: an two dimension array for the splitting surfaces' coordinates 
//      num_of_split_surfs: the number of splitting surfaces along each dimension

//Output:
//      segment_length: the segments' lengths and owners of the ray and the voxel model, which is specified by the splitting surfaces
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{

	double minCoords[3] = { split_surf_coords[LONDIM][0],split_surf_coords[LATDIM][0],split_surf_coords[RDIM][0] };
	double maxCoords[3] = { split_surf_coords[LONDIM][num_of_split_surfs[LONDIM] - 1],split_surf_coords[LATDIM][num_of_split_surfs[LATDIM] - 1],split_surf_coords[RDIM][num_of_split_surfs[RDIM] - 1] };
	double t_firstLast[2];


	if (!is_a_valid_ray(start_pt, end_pt, minCoords, maxCoords, t_firstLast)) return 0;
	double  start_pt_xyz[3], end_pt_xyz[3], end_pt_llr[3], start_pt_llr[3];

	XYZ_from_t_vec(start_pt, end_pt, t_firstLast[0], start_pt_xyz);
	XYZ_from_t_vec(start_pt, end_pt, t_firstLast[1], end_pt_xyz);
	XYZ_to_LLR(end_pt_xyz, end_pt_llr);
	XYZ_to_LLR(start_pt_xyz, start_pt_llr);
	///////////////////////////////////////////////////////////////////////////////////////////////////
	//calculate the index range of the valid splitting surfaces
	///////////////////////////////////////////////////////////////////////////////////////////////////
	int index_of_low_bound_coord[3], index_of_high_bound_coord[3];
	INDEX_LENGTH low_pt_voxel;
	for (int i = 0; i < 3; i++)
	{

		if (start_pt_llr[i] < end_pt_llr[i])
		{
			index_of_low_bound_coord[i] = lower_bound_index(split_surf_coords[i], num_of_split_surfs[i], start_pt_llr[i]);
			if (index_of_low_bound_coord[i] == -1) index_of_low_bound_coord[i] = 0;//in case of any error in calculating an intersection
			low_pt_voxel.index[i] = index_of_low_bound_coord[i];

			index_of_high_bound_coord[i] = upper_bound_index(split_surf_coords[i], num_of_split_surfs[i], end_pt_llr[i]);
			if (index_of_high_bound_coord[i] == -1) index_of_high_bound_coord[i] = num_of_split_surfs[i] - 1;

		}
		else
		{
			index_of_low_bound_coord[i] = lower_bound_index(split_surf_coords[i], num_of_split_surfs[i], end_pt_llr[i]);
			if (index_of_low_bound_coord[i] == -1) index_of_low_bound_coord[i] = 0;

			index_of_high_bound_coord[i] = upper_bound_index(split_surf_coords[i], num_of_split_surfs[i], start_pt_llr[i]);
			if (index_of_high_bound_coord[i] == -1) index_of_high_bound_coord[i] = num_of_split_surfs[i] - 1;

			low_pt_voxel.index[i] = lower_bound_index(split_surf_coords[i], num_of_split_surfs[i], start_pt_llr[i]);
			if (low_pt_voxel.index[i] == -1) low_pt_voxel.index[i] = 0;

		}
	}
	low_pt_voxel.t_vec_first_pt = 0;
	low_pt_voxel.intercept = -1;
	///////////////////////////////////////////////////////////////////////////////////////////////////
	int orientVector[3];
	double t_tangent[2];
	bool twoVal[2];
	orientation_vector_of_a_ray(start_pt_xyz, end_pt_xyz, orientVector, twoVal, t_tangent);


	//////////////////////////////////////////////////////////////////////////////
	//calcuate the t-value for each valid intersection 
	/////////////////////////////////////////////////////////////////////////////
	int t_num_lon = index_of_high_bound_coord[LONDIM] - index_of_low_bound_coord[LONDIM] + 1;
	int t_num_lat = index_of_high_bound_coord[LATDIM] - index_of_low_bound_coord[LATDIM] + 1;
	int t_num_rad = index_of_high_bound_coord[RDIM] - index_of_low_bound_coord[RDIM] + 1;

	int t_num_all_possible = t_num_lon + t_num_lat * 2 + t_num_rad * 2;
	NEXTVISIT* t_vec_at_each_surf = new NEXTVISIT[t_num_all_possible];

	int t_num_total = 0;
	for (int i = index_of_low_bound_coord[LONDIM]; i <= index_of_high_bound_coord[LONDIM]; i++)
	{
		double t_vector;
		int num = t_vec_to_longitude_surface(start_pt_xyz, end_pt_xyz, split_surf_coords[LONDIM][i], t_vector);
		if (num < 1) continue;// no intersection
		t_vec_at_each_surf[t_num_total].t_vector = t_vector;
		t_vec_at_each_surf[t_num_total].next_visit[LONDIM] = orientVector[LONDIM];
		t_vec_at_each_surf[t_num_total].next_visit[LATDIM] = 0;
		t_vec_at_each_surf[t_num_total++].next_visit[RDIM] = 0;
	}
	for (int i = index_of_low_bound_coord[LATDIM]; i <= index_of_high_bound_coord[LATDIM]; i++)
	{

		double t_vector[2];
		int num = t_vec_to_latitude_surface(start_pt_xyz, end_pt_xyz, split_surf_coords[LATDIM][i], t_vector);
		if (num < 1) continue;//no intersection
		double tmpxyz[3], tmpllr[3];
		XYZ_from_t_vec(start_pt_xyz, end_pt_xyz, t_vector[0], tmpxyz);
		XYZ_to_LLR(tmpxyz, tmpllr);
		t_vec_at_each_surf[t_num_total].t_vector = t_vector[0];
		t_vec_at_each_surf[t_num_total].next_visit[LONDIM] = 0;
		t_vec_at_each_surf[t_num_total].next_visit[LATDIM] = orientVector[LATDIM];
		if (twoVal[0] && t_vector[0] > t_tangent[0])//current intersection locals on the second part of the ray
			t_vec_at_each_surf[t_num_total].next_visit[LATDIM] = -orientVector[LATDIM];
		t_vec_at_each_surf[t_num_total++].next_visit[RDIM] = 0;
		if (num == 2)
		{
			t_vec_at_each_surf[t_num_total].t_vector = t_vector[1];
			t_vec_at_each_surf[t_num_total].next_visit[LONDIM] = 0;
			t_vec_at_each_surf[t_num_total].next_visit[LATDIM] = orientVector[LATDIM];
			if (twoVal[0] && t_vector[0] > t_tangent[0])//current intersection locals on the second part of the ray
				t_vec_at_each_surf[t_num_total].next_visit[LATDIM] = -orientVector[LATDIM];
			t_vec_at_each_surf[t_num_total++].next_visit[RDIM] = 0;
		}
	}
	
	for (int i = index_of_low_bound_coord[RDIM]; i <= index_of_high_bound_coord[RDIM]; i++)
	{
		double t_vector[2];
		int num = t_vec_to_sphere(start_pt_xyz, end_pt_xyz, split_surf_coords[RDIM][i], t_vector);
		if (num < 1) continue;//no intersection
		t_vec_at_each_surf[t_num_total].t_vector = t_vector[0];
		t_vec_at_each_surf[t_num_total].next_visit[LONDIM] = 0;
		t_vec_at_each_surf[t_num_total].next_visit[LATDIM] = 0;
		t_vec_at_each_surf[t_num_total++].next_visit[RDIM] = orientVector[RDIM];
		if (num == 2)
		{
			t_vec_at_each_surf[t_num_total].t_vector = t_vector[1];
			t_vec_at_each_surf[t_num_total].next_visit[LONDIM] = 0;
			t_vec_at_each_surf[t_num_total].next_visit[LATDIM] = 0;
			t_vec_at_each_surf[t_num_total++].next_visit[RDIM] = -orientVector[RDIM];
		}

	}
	/////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	//sort the intersections according to their t-values
	/////////////////////////////////////////////////////////////////////////////
	std::sort(t_vec_at_each_surf, t_vec_at_each_surf + t_num_total, cmpfunction0);
	/////////////////////////////////////////////////////////////////////////////	


	//////////////////////////////////////////////////////////////////////////////
	////find the index of the first t-value , because that the t-value may smaller than that of the first intersection due to calucation error
	//////////////////////////////////////////////////////////////////////////////
	int first_index = 0;
	for (int i = 0; i < t_num_total; i++)
		if (t_vec_at_each_surf[i].t_vector - low_pt_voxel.t_vec_first_pt > THRESHOLD_T_OF_TWO_POINTS)
		{
			first_index = i;
			break;
		}
	//////////////////////////////////////////////////////////////////////////////


	double dist_vector = sqrt((start_pt_xyz[XDIM] - end_pt_xyz[XDIM])*(start_pt_xyz[XDIM] - end_pt_xyz[XDIM]) + (start_pt_xyz[YDIM] - end_pt_xyz[YDIM])
		*(start_pt_xyz[YDIM] - end_pt_xyz[YDIM]) + (start_pt_xyz[ZDIM] - end_pt_xyz[ZDIM])*(start_pt_xyz[ZDIM] - end_pt_xyz[ZDIM]));
	//get any two consecutive t-values, and calucate the length and the owner's index of the segment that are defined by these two t-values.
	int num = t_vector_to_index_and_length(t_vec_at_each_surf + first_index, t_num_total - first_index, dist_vector, low_pt_voxel, num_of_split_surfs, segment_length);
	delete[] t_vec_at_each_surf;
	return num;

}