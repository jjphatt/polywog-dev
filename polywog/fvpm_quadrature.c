// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "integrators/gauss_rules.h"
#include "polywog/fvpm_quadrature.h"

surface_integral_t* fvpm_cube_surface_integral_new(point_cloud_t* cloud,
                                                   neighbor_pairing_t* pairing,
                                                   real_t* extents,
                                                   int num_points,
                                                   real_t side_to_extent_ratio)
{
}

volume_integral_t* fvpm_cube_volume_integral_new(point_cloud_t* cloud,
                                                 neighbor_pairing_t* pairing,
                                                 real_t* extents,
                                                 int num_points,
                                                 real_t side_to_extent_ratio)
{
}

surface_integral_t* fvpm_sphere_surface_integral_new(point_cloud_t* cloud,
                                                     neighbor_pairing_t* pairing,
                                                     real_t* extents,
                                                     int num_points,
                                                     real_t radius_to_extent_ratio)
{
}

volume_integral_t* fvpm_sphere_volume_integral_new(point_cloud_t* cloud,
                                                   neighbor_pairing_t* pairing,
                                                   real_t* extents,
                                                   int num_points,
                                                   real_t radius_to_extent_ratio)
{
}

