#pragma once
#include "common.h"
#include <math.h>
#include <assert.h>
#include <algorithm>

void orientation_vector_of_a_ray(double low_pt_xyz[3], double high_pt_xyz[3], int orientVector[3], bool twoVal[2], double t_tangent[2]);

int t_vector_to_Index_length(NEXTVISIT* t_vec_each_surf, int num_of_t_vector, double lenOfpath, const INDEX_LENGTH& firstVoxel, int num_of_coords[3], INDEX_LENGTH*& voxels_intercepts);
int segments_of_a_ray_by_AP(double pt1[3], double pt2[3], double* coords[3], int num_of_coords[3], INDEX_LENGTH*& voxel_intercepts);