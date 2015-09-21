// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYWOG_FVPM_QUADRATURE_H
#define POLYWOG_FVPM_QUADRATURE_H

#include "core/point_cloud.h"
#include "integrators/volume_integral.h"
#include "integrators/surface_integral.h"

// This file contains several quadrature rules for evaluating volume and 
// surface integrals using the Finite Volume Point Method, which involves 
// integrals over intersections of subdomain pairs. These quadrature rules 
// are used in the construction of conservative meshless methods.

// All of these quadrature rules are associated with a point cloud and a 
// neighbor pairing, and the "ith" integral refers to the integral of the 
// region of intersection over the ith pair of points in the point cloud, 
// using the indexing scheme defined by the neighbor pairing.

// This quadrature rule computes surface integrals over the intersections of 
// pairs of cubes whose extents are defined by the ratio of its side to the 
// extent associated with a subdomain. The rule is associated with the given 
// point cloud, neighbor pairing, and (scalar) subdomain extents field, and 
// has the given number of points on a side.
surface_integral_t* fvpm_cube_surface_integral_new(point_cloud_t* cloud,
                                                   neighbor_pairing_t* pairing,
                                                   real_t* extents,
                                                   int num_points,
                                                   real_t side_to_extent_ratio);

// This quadrature rule computes volume integrals over the intersections of 
// pairs of cubes whose extents are defined by the ratio of its side to the 
// extent associated with a subdomain. The rule is associated with the given 
// point cloud, neighbor pairing, and (scalar) subdomain extents field, and 
// has the given number of points on a side.
volume_integral_t* fvpm_cube_volume_integral_new(point_cloud_t* cloud,
                                                 neighbor_pairing_t* pairing,
                                                 real_t* extents,
                                                 int num_points,
                                                 real_t side_to_extent_ratio);

// This quadrature rule computes surface integrals over the intersections of 
// pairs of spheres whose extents are defined by the ratio of its radius to 
// the extent associated with a subdomain. The rule is associated with the 
// given point cloud, neighbor pairing, and (scalar) subdomain extents field, 
// and has the given number of azimuthal points.
surface_integral_t* fvpm_sphere_surface_integral_new(point_cloud_t* cloud,
                                                     neighbor_pairing_t* pairing,
                                                     real_t* extents,
                                                     int num_points,
                                                     real_t radius_to_extent_ratio);

// This quadrature rule computes volume integrals over the intersections of 
// spheres whose extents are defined by the ratio of its radius to the extent 
// associated with a subdomain. The rule is associated with the given point 
// cloud, neighbor pairing, and (scalar) subdomain extents field, and has the 
// given number of radial points.
volume_integral_t* fvpm_sphere_volume_integral_new(point_cloud_t* cloud,
                                                   neighbor_pairing_t* pairing,
                                                   real_t* extents,
                                                   int num_points,
                                                   real_t radius_to_extent_ratio);

#endif
