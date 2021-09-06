#pragma once
#include "common.h"
#include <vector>

double segment_length_in_a_voxel(double minCoords[3], double maxCoords[3], double startPt[3], double endPt[3], double vectNormal);
int segments_of_a_ray_by_tradition(double startPt[3], double endPt[3], double* coords[3], int num_of_coords[3], std::vector<INDEX_LENGTH*>& segments);
int segments_of_a_ray_by_AP(double startPt[3], double endPt[3], double* coords[3], int num_of_coords[3], std::vector<INDEX_LENGTH*>& segments);
