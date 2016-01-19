// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYWOG_SPH_H_UPDATER_H
#define POLYWOG_SPH_H_UPDATER_H

#include "polywog/sph_kernel.h"

// The SPH dynamics "H updater" is an object that computes the computes of 
// the smoothing tensor field H by examining the neighborhoods surrounding 
// each point, and adjusting until a desired number of neighbors per smoothing
// length is achieved along each principal axis. The approach for computing 
// this "ideal" value of H is described by Owen (Ref?).
typedef struct sph_H_updater_t sph_H_updater_t;

// Creates an SPH H updater object that uses the given number of neighbors 
// per smoothing length (n_per_h). This update should be performed at the 
// initialization of an SPH simulation and at the end of each time step.
// The scalar smoothing scale h is computed and placed in a diagonal H tensor.
sph_H_updater_t* isotropic_sph_H_updater_new(sph_kernel_t* W,
                                             real_t n_per_h);

// This version of the H updater computes the components of the anisotropic 
// H tensor used in ASPH. 
sph_H_updater_t* anisotropic_sph_H_updater_new(sph_kernel_t* W,
                                               real_t n_per_h);

// Destroys the given SPH H updater object.
void sph_H_updater_free(sph_H_updater_t* updater);

// Computes the new value of H at the point x given the neighbors {y}. The 
// old value of H may be placed in H to improve the performance of the update.
// The number of iterations and the maximum fractional change in H are 
// returned as diagnostics. This function returns true if the update succeeded 
// and false otherwise.
bool sph_H_updater_update(sph_H_updater_t* updater, 
                          point_t* x, point_t* ys, int num_neighbors,
                          sym_tensor2_t* H,
                          int* num_iterations,
                          real_t* max_fractional_change);

// Sets the maximum number of iterations used to find the new value of H.
void sph_H_updater_set_max_iterations(sph_H_updater_t* updater, 
                                      int max_iters);

// Sets the threshold for the fractional change in H that determines when an 
// update is completed.
void sph_H_updater_set_convergence_threshold(sph_H_updater_t* updater, 
                                             real_t fractional_change);

#endif
