// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polywog/fvpm_interparticle_area.h"

struct fvpm_interparticle_area_t 
{
  point_cloud_t* points;
  real_t* extents;
  real_t overlap;
  point_spacing_estimator_t* dx_estimator;
};

fvpm_interparticle_area_t* sphere_fvpm_interparticle_area_new(point_cloud_t* cloud,
                                                              real_t* extents,
                                                              real_t overlap_parameter,
                                                              point_spacing_estimator_t* dx_estimator)
{
  ASSERT(overlap_parameter > 0.0);

  fvpm_interparticle_area_t* area = polymec_malloc(sizeof(fvpm_interparticle_area_t));
  area->points = cloud;
  area->extents = extents;
  area->overlap = overlap_parameter;
  area->dx_estimator = dx_estimator;
  return area;
}

fvpm_interparticle_area_t* cube_fvpm_interparticle_area_new(point_cloud_t* cloud,
                                                            real_t* extents,
                                                            real_t overlap_parameter,
                                                            point_spacing_estimator_t* dx_estimator)
{
  ASSERT(overlap_parameter > 0.0);

  fvpm_interparticle_area_t* area = polymec_malloc(sizeof(fvpm_interparticle_area_t));
  area->points = cloud;
  area->extents = extents;
  area->overlap = overlap_parameter;
  area->dx_estimator = dx_estimator;
  return area;
}

void fvpm_interparticle_area_free(fvpm_interparticle_area_t* area)
{
  area->dx_estimator = NULL;
  polymec_free(area);
}

void fvpm_interparticle_area_compute(fvpm_interparticle_area_t* area, 
                                     int i, int j,
                                     vector_t* beta_ij)
{
}

