// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "integrators/gauss_rules.h"
#include "polywog/fvpm_quadrature.h"

typedef struct
{
  point_cloud_t* cloud;
  neighbor_pairing_t* pairing;
  real_t* extents;
  int N;
  real_t ratio; 
} fvpm_simple_t;

static fvpm_simple_t* fvpm_simple_new(point_cloud_t* cloud,
                                      neighbor_pairing_t* pairing,
                                      real_t* extents,
                                      int num_points,
                                      real_t ratio)
{
  ASSERT(num_points > 0);
  ASSERT(ratio > 0.0);

  fvpm_simple_t* fvpm = polymec_malloc(sizeof(fvpm_simple_t));
  fvpm->cloud = cloud;
  fvpm->pairing = pairing;
  fvpm->extents = extents;
  fvpm->N = num_points;
  fvpm->ratio = ratio;
  return fvpm;
}

static void fvpm_simple_free(void* context)
{
  polymec_free(context);
}

static inline void get_cubes(fvpm_simple_t* fvpm, int i, int j, bbox_t* boxi, bbox_t* boxj)
{
  // Construct the two cubes and the intersection of their surfaces. This can 
  // be a point, a line segment, or a rectangular plane segment.
  point_t* xi = &fvpm->cloud->points[i];
  real_t hi = fvpm->extents[i];
  real_t Li = fvpm->ratio * hi;
  boxi->x1 = xi->x - 0.5*Li;
  boxi->x2 = xi->x + 0.5*Li;
  boxi->y1 = xi->y - 0.5*Li; 
  boxi->y2 = xi->y + 0.5*Li;
  boxi->z1 = xi->z - 0.5*Li;
  boxi->z2 = xi->z + 0.5*Li;

  point_t* xj = &fvpm->cloud->points[j];
  real_t hj = fvpm->extents[j];
  real_t Lj = fvpm->ratio * hj;
  boxj->x1 = xj->x - 0.5*Lj;
  boxj->x2 = xj->x + 0.5*Lj;
  boxj->y1 = xj->y - 0.5*Lj; 
  boxj->y2 = xj->y + 0.5*Lj;
  boxj->z1 = xj->z - 0.5*Lj;
  boxj->z2 = xj->z + 0.5*Lj;
}

static int cube_surf_num_quad_points(void* context, int k)
{
  fvpm_simple_t* fvpm = context;
  int i, j;
  neighbor_pairing_get(fvpm->pairing, k, &i, &j, NULL);
  bbox_t boxi, boxj, box_int;
  get_cubes(fvpm, i, j, &boxi, &boxj);
  get_cubes(fvpm, i, j, &boxi, &boxj);
  bbox_intersect_bbox(&boxi, &boxj, &box_int);
  if (bbox_is_empty_set(&box_int) || bbox_is_point(&box_int) || bbox_is_line(&box_int))
    return 0;
  else
    return fvpm->N * fvpm->N;
}

static void cube_surf_get_quad(void* context, int k, point_t* points, real_t* weights, vector_t* normals)
{
  fvpm_simple_t* fvpm = context;
  real_t gauss_pts[fvpm->N], gauss_wts[fvpm->N];
  get_gauss_legendre_points(fvpm->N, gauss_wts, gauss_wts);

  // Get the bounding boxes for the kth pair.
  int i, j;
  bbox_t boxi, boxj, box_int;
  get_cubes(fvpm, i, j, &boxi, &boxj);
  bbox_intersect_bbox(&boxi, &boxj, &box_int);
  ASSERT(!bbox_is_empty_set(&box_int) && !bbox_is_point(&box_int) && !bbox_is_line(&box_int));

  // FIXME: It remains to implement the quadrature for this rule.
}

surface_integral_t* fvpm_cube_surface_integral_new(point_cloud_t* cloud,
                                                   neighbor_pairing_t* pairing,
                                                   real_t* extents,
                                                   int num_points,
                                                   real_t side_to_extent_ratio)
{
  fvpm_simple_t* fvpm = fvpm_simple_new(cloud, pairing, extents, num_points, side_to_extent_ratio);
  surface_integral_vtable vtable = {.num_quad_points = cube_surf_num_quad_points,
                                    .get_quadrature = cube_surf_get_quad,
                                    .dtor = fvpm_simple_free};
  char name[1025];
  snprintf(name, 1024, "FVPM cube surface integral (N = %d, side/extent = %g)", 
           num_points, side_to_extent_ratio);
  return surface_integral_new(name, fvpm, vtable);
}

static int cube_vol_num_quad_points(void* context, int k)
{
  fvpm_simple_t* fvpm = context;
  int i, j;
  neighbor_pairing_get(fvpm->pairing, k, &i, &j, NULL);
  bbox_t boxi, boxj, box_int;
  get_cubes(fvpm, i, j, &boxi, &boxj);
  get_cubes(fvpm, i, j, &boxi, &boxj);
  bbox_intersect_bbox(&boxi, &boxj, &box_int);
  if (bbox_is_empty_set(&box_int) || bbox_is_point(&box_int) || bbox_is_line(&box_int) || bbox_is_plane(&box_int))
    return 0;
  else
    return fvpm->N * fvpm->N * fvpm->N;
}

static void cube_vol_get_quad(void* context, int k, point_t* points, real_t* weights)
{
  fvpm_simple_t* fvpm = context;
  real_t gauss_pts[fvpm->N], gauss_wts[fvpm->N];
  get_gauss_legendre_points(fvpm->N, gauss_wts, gauss_wts);

  // Get the bounding boxes for the kth pair.
  int i, j;
  bbox_t boxi, boxj, box_int;
  get_cubes(fvpm, i, j, &boxi, &boxj);
  bbox_intersect_bbox(&boxi, &boxj, &box_int);
  ASSERT(!bbox_is_empty_set(&box_int) && !bbox_is_point(&box_int) && !bbox_is_line(&box_int) && !bbox_is_plane(&box_int));

  // FIXME: It remains to implement the quadrature for this rule.
}

volume_integral_t* fvpm_cube_volume_integral_new(point_cloud_t* cloud,
                                                 neighbor_pairing_t* pairing,
                                                 real_t* extents,
                                                 int num_points,
                                                 real_t side_to_extent_ratio)
{
  fvpm_simple_t* fvpm = fvpm_simple_new(cloud, pairing, extents, num_points, side_to_extent_ratio);
  volume_integral_vtable vtable = {.num_quad_points = cube_vol_num_quad_points,
                                   .get_quadrature = cube_vol_get_quad,
                                   .dtor = fvpm_simple_free};
  char name[1025];
  snprintf(name, 1024, "FVPM cube volume integral (N = %d, side/extent = %g)", 
           num_points, side_to_extent_ratio);
  return volume_integral_new(name, fvpm, vtable);
}

static int sphere_surf_num_quad_points(void* context, int k)
{
  fvpm_simple_t* fvpm = context;
  int i, j;
  neighbor_pairing_get(fvpm->pairing, k, &i, &j, NULL);
  bbox_t boxi, boxj, box_int;
  get_cubes(fvpm, i, j, &boxi, &boxj);
  get_cubes(fvpm, i, j, &boxi, &boxj);
  bbox_intersect_bbox(&boxi, &boxj, &box_int);
  if (bbox_is_empty_set(&box_int) || bbox_is_point(&box_int) || bbox_is_line(&box_int))
    return 0;
  else
    return fvpm->N * fvpm->N;
}

static void sphere_surf_get_quad(void* context, int k, point_t* points, real_t* weights, vector_t* normals)
{
  fvpm_simple_t* fvpm = context;
  real_t gauss_pts[fvpm->N], gauss_wts[fvpm->N];
  get_gauss_legendre_points(fvpm->N, gauss_wts, gauss_wts);

  // Get the bounding boxes for the kth pair.
  int i, j;
  bbox_t boxi, boxj, box_int;
  get_cubes(fvpm, i, j, &boxi, &boxj);
  bbox_intersect_bbox(&boxi, &boxj, &box_int);
  ASSERT(!bbox_is_empty_set(&box_int) && !bbox_is_point(&box_int) && !bbox_is_line(&box_int));

  // FIXME: It remains to implement the quadrature for this rule.
}

surface_integral_t* fvpm_sphere_surface_integral_new(point_cloud_t* cloud,
                                                     neighbor_pairing_t* pairing,
                                                     real_t* extents,
                                                     int num_points,
                                                     real_t radius_to_extent_ratio)
{
  fvpm_simple_t* fvpm = fvpm_simple_new(cloud, pairing, extents, num_points, radius_to_extent_ratio);
  surface_integral_vtable vtable = {.num_quad_points = sphere_surf_num_quad_points,
                                    .get_quadrature = sphere_surf_get_quad,
                                    .dtor = fvpm_simple_free};
  char name[1025];
  snprintf(name, 1024, "FVPM sphere surface integral (N = %d, side/extent = %g)", 
           num_points, radius_to_extent_ratio);
  return surface_integral_new(name, fvpm, vtable);
}

static int sphere_vol_num_quad_points(void* context, int k)
{
  fvpm_simple_t* fvpm = context;
  int i, j;
  neighbor_pairing_get(fvpm->pairing, k, &i, &j, NULL);
  bbox_t boxi, boxj, box_int;
  get_cubes(fvpm, i, j, &boxi, &boxj);
  get_cubes(fvpm, i, j, &boxi, &boxj);
  bbox_intersect_bbox(&boxi, &boxj, &box_int);
  if (bbox_is_empty_set(&box_int) || bbox_is_point(&box_int) || bbox_is_line(&box_int) || bbox_is_plane(&box_int))
    return 0;
  else
    return fvpm->N * fvpm->N * fvpm->N;
}

static void sphere_vol_get_quad(void* context, int k, point_t* points, real_t* weights)
{
  fvpm_simple_t* fvpm = context;
  real_t gauss_pts[fvpm->N], gauss_wts[fvpm->N];
  get_gauss_legendre_points(fvpm->N, gauss_wts, gauss_wts);

  // Get the bounding boxes for the kth pair.
  int i, j;
  bbox_t boxi, boxj, box_int;
  get_cubes(fvpm, i, j, &boxi, &boxj);
  bbox_intersect_bbox(&boxi, &boxj, &box_int);
  ASSERT(!bbox_is_empty_set(&box_int) && !bbox_is_point(&box_int) && !bbox_is_line(&box_int) && !bbox_is_plane(&box_int));

  // FIXME: It remains to implement the quadrature for this rule.
}

volume_integral_t* fvpm_sphere_volume_integral_new(point_cloud_t* cloud,
                                                   neighbor_pairing_t* pairing,
                                                   real_t* extents,
                                                   int num_points,
                                                   real_t radius_to_extent_ratio)
{
  fvpm_simple_t* fvpm = fvpm_simple_new(cloud, pairing, extents, num_points, radius_to_extent_ratio);
  volume_integral_vtable vtable = {.num_quad_points = sphere_vol_num_quad_points,
                                    .get_quadrature = sphere_vol_get_quad,
                                    .dtor = fvpm_simple_free};
  char name[1025];
  snprintf(name, 1024, "FVPM sphere volume integral (N = %d, side/extent = %g)", 
           num_points, radius_to_extent_ratio);
  return volume_integral_new(name, fvpm, vtable);
}

