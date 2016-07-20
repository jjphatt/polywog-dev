// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYWOG_MULTICLOUD_H
#define POLYWOG_MULTICLOUD_H

#include "core/point_cloud.h"

// This class represents a multi-level solver that is able to create a 
// hierarchy of coarsened point clouds, given an original cloud, and 
// apply a Full Approximation Storage (FAS) method to solve a linear or 
// nonlinear problem on that hierarchy. Objects of this type are garbage-
// collected.
typedef struct multicloud_t multicloud_t;

// This class represents an iterative method that updates the solution to a 
// given linear or nonlinear problem on the points of a point cloud. Objects 
// of this type are garbage-collected.
typedef struct multicloud_iteration_t multicloud_iteration_t;

// This virtual table defines the behaviors required to implement a 
// multicloud iteration scheme.
typedef struct
{
  // Calculates the residual of the solution X for the given cloud at the 
  // given time t, storing it in the vector R.
  void (*residual)(void* context, point_cloud_t* cloud, real_t t, real_t* X, real_t* R);
  // Updates the solution X for the given cloud at the given time t.
  void (*update)(void* context, point_cloud_t* cloud, real_t t, real_t* X);
  // Destructor.
  void (*dtor)(void* context); 
} multicloud_iteration_vtable;

// This class represents a hierarchy of coarsened point clouds that are 
// created, manipulated, and used by multicloud solvers.
typedef struct multicloud_hierarchy_t multicloud_hierarchy_t;

// Creates a multicloud solver that uses the given iteration scheme to solve 
// a linear/nonlinear problem on the points of a point cloud.
multicloud_t* multicloud_new(multicloud_iteration_t* iteration_scheme);

// Creates and returns a coarsened multicloud hierarchy for the given point 
// cloud. This hierarchy can be used to perform calculations on the point 
// cloud as long as its relationship to its coarse counterparts remains the 
// same as when the hierarchy was created.
multicloud_hierarchy_t* multicloud_coarsen(multicloud_t* mc,
                                           point_cloud_t* cloud);

// Finds the solution corresponding to the points in the finest member of 
// the given multicloud hierarchy, and given the initial solution X. The 
// solution vector X is updated in place. The number of iterations is placed 
// in num_iters, and the residual norm is placed in res_norm.
void multicloud_solve(multicloud_t* mc,
                      multicloud_hierarchy_t* hierarchy,
                      real_t* X,
                      int* num_iters,
                      real_t* res_norm);

// Creates a new multicloud iteration scheme with the given context pointer 
// and virtual table.
multicloud_iteration_t* multicloud_iteration_new(const char* name,
                                                 void* context,
                                                 multicloud_iteration_vtable vtable);

// Creates a new coarsened point cloud hierarchy whose finest level is the 
// given point cloud. The hierarchy assumes control of all point clouds that 
// it contains EXCEPT this finest cloud.
multicloud_hierarchy_t* multicloud_hierarchy_new(point_cloud_t* cloud);

// Destroys the given point cloud hierarchy, destroying all point clouds 
// EXCEPT the finest one (with which it was constructed).
void multicloud_hierarchy_free(multicloud_hierarchy_t* hierarchy);

// Appends the given point cloud to the multicloud hierarchy as the next 
// coarsest level. The point cloud assumes control of this cloud.
void multicloud_hierarchy_append(multicloud_hierarchy_t* hierarchy,
                                 point_cloud_t* coarser_cloud);

#endif
