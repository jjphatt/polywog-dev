// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYWOG_SPH_KERNEL_H
#define POLYWOG_SPH_KERNEL_H

#include "core/point.h"
#include "core/lookup1.h"

// An SPH kernel is a weighting function that determines the strength of 
// contributions from particles in a neighbor in a Smoothed Particle 
// Hydrodynamics calculation. Objects of this type are garbage-collected.
typedef struct sph_kernel_t sph_kernel_t;

// Creates an SPH kernel object with the given name, extent (in number of 
// smoothing lengths), context, compute function, and destructor.
sph_kernel_t* sph_kernel_new(const char* name,
                             real_t extent,
                             void* context,
                             void (*compute)(void* context, real_t eta_mag, real_t det_H, real_t* W, real_t* dWdeta),
                             void (*dtor)(void* context));

// Returns an internal pointer to a string containing the name of the kernel.
char* sph_kernel_name(sph_kernel_t* kernel);

// Returns the extent of the support of the SPH kernel, in units of 
// smoothing lengths.
real_t sph_kernel_extent(sph_kernel_t* kernel);

// Computes the value and gradient of the kernel at a displacement of x from 
// center, given the symmetric smoothing tensor H.
void sph_kernel_compute(sph_kernel_t* kernel, vector_t* x, sym_tensor2_t* H, real_t* W, vector_t* grad_W);

// Returns the sum of the contributions from this kernel in a neighborhood
// about a point with n_per_h neighbors per SPH smoothing scale.
real_t sph_kernel_sum(sph_kernel_t* kernel, real_t n_per_h);

// Creates and returns a cubic B-spline SPH kernel.
sph_kernel_t* b_spline_sph_kernel_new();

// Creates a kernel that uses a lookup table with linear or quadratic 
// interpolation with the given resolution to quickly evaluate values 
// precomputed by the given kernel.
sph_kernel_t* tabular_sph_kernel_new(sph_kernel_t* kernel, 
                                     lookup1_interpolation_t interpolation,
                                     int resolution);

#endif
