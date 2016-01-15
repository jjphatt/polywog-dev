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

// Creates an SPH kernel object with the given name, that uses the given 
// function to compute time derivatives corresponding to terms in solution 
// vectors for points i and j, given the time t, solution vectors Ui and Uj, 
// values of the SPH kernel Wi and Wj, and their vector-valued gradients 
// grad_Wi and grad_Wj. The results of the computation for i and j are placed 
// in dUidt and dUjdt. 
sph_kernel_t* sph_kernel_new(const char* name,
                             void* context,
                             void (*compute)(void* context, point_t* x, real_t* H, real_t* W, vector_t* grad_W),
                             void (*dtor)(void* context));

// Destroys the given SPH kernel object.
void sph_kernel_free(sph_kernel_t* kernel);

// Computes the value and gradient of the kernel at the point x, given the 
// symmetric smoothing tensor H at that point.
void sph_kernel_compute(sph_kernel_t* kernel, point_t* x, real_t* H, real_t* W, vector_t* grad_W);

#endif
