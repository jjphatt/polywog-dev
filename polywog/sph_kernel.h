// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYWOG_SPH_KERNEL_H
#define POLYWOG_SPH_KERNEL_H

#include "core/point.h"

// An SPH kernel is a weighting function that determines the strength of 
// contributions from particles in a neighbor in a Smoothed Particle 
// Hydrodynamics calculation.
typedef struct sph_kernel_t sph_kernel_t;

// Creates an SPH kernel object with the given name, context, and compute 
// function (and destructor).
sph_kernel_t* sph_kernel_new(const char* name,
                             void* context,
                             void (*compute)(void* context, vector_t* x, sym_tensor2_t* H, real_t* W, vector_t* grad_W),
                             void (*dtor)(void* context));

// Destroys the given SPH kernel object.
void sph_kernel_free(sph_kernel_t* kernel);

// Computes the value and gradient of the kernel at a displacement of x from 
// center, given the symmetric smoothing tensor H.
void sph_kernel_compute(sph_kernel_t* kernel, vector_t* x, sym_tensor2_t* H, real_t* W, vector_t* grad_W);

// Creates and returns a cubic B-spline SPH kernel.
sph_kernel_t* b_spline_sph_kernel_new();

// Creates and returns an SPH kernel that uses a look­up table to retrieve
// values from the given SPH kernel.
sph_kernel_t* tabular_sph_kernel_new(sph_kernel_t* kernel);

#endif
