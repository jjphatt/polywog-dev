// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYWOG_FVPM_INTERPARTICLE_AREA_H
#define POLYWOG_FVPM_INTERPARTICLE_AREA_H

#include "core/point_cloud.h"
#include "model/neighbor_pairing.h"
#include "model/point_spacing_estimator.h"

// This class computes interparticle areas (denoted as beta_ij in the FVPM
// literature, e.g. Quinlan et al / CPC 185 (2014) 1554. It exactly integrates
// 
typedef struct fvpm_interparticle_area_t fvpm_interparticle_area_t;

// Creates an object to compute interparticle areas for a distribution of 
// points defined by a point cloud, using the given extent field to associate 
// spheres of radius h with points, the given particle overlap parameter h/dx,
// and the given point spacing estimator, which estimates dx in the vicinity 
// of a point.
fvpm_interparticle_area_t* sphere_fvpm_interparticle_area_new(point_cloud_t* cloud,
                                                              real_t* extents,
                                                              real_t overlap_parameter,
                                                              point_spacing_estimator_t* estimator);

// Creates an object to compute interparticle areas for a distribution of 
// points defined by a point cloud, using the given extent field to associate 
// cubes of side h with points, with the given particle overlap 
// parameter h/dx, and the given point spacing estimator, which estimates dx 
// in the vicinity of a point.
fvpm_interparticle_area_t* cube_fvpm_interparticle_area_new(point_cloud_t* cloud,
                                                            real_t* extents,
                                                            real_t overlap_parameter,
                                                            point_spacing_estimator_t* estimator);

// Frees the given interparticle area object.
void fvpm_interparticle_area_free(fvpm_interparticle_area_t* area);

// Computes the interparticle area vector beta_ij for the given nodes i and j, 
// placing the components into beta_ij.
void fvpm_interparticle_area_compute(fvpm_interparticle_area_t* area, 
                                     int i, int j,
                                     vector_t* beta_ij);

#endif
