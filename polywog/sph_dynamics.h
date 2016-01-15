// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYWOG_SPH_DYNAMICS_H
#define POLYWOG_SPH_DYNAMICS_H

#include "core/point.h"

// The SPH dynamics interface represents contributions to time evolution 
// equations resulting from a pairwise interaction between two points i and j.
typedef struct sph_dynamics_t sph_dynamics_t;

// Creates an SPH dynamics object with the given name, that uses the given 
// function to compute time derivatives corresponding to terms in solution 
// vectors for points i and j, given the time t, solution vectors Ui and Uj, 
// values of the SPH kernel Wi and Wj, and their vector-valued gradients 
// grad_Wi and grad_Wj. The results of the computation for i and j are placed 
// in dUidt and dUjdt. 
sph_dynamics_t* sph_dynamics_new(const char* name,
                                 void* context,
                                 void (*compute)(void* context, real_t t,
                                                 int i, int j,
                                                 real_t* Ui, real_t* Uj, 
                                                 real_t Wi, real_t Wj,
                                                 vector_t* grad_Wi, vector_t* grad_Wj,
                                                 real_t* dUidt, real_t* dUjdt),
                                 void (*dtor)(void* context));

// Destroys the given SPH dynamics object.
void sph_dynamics_free(sph_dynamics_t* dyn);

// Computes the time derivatives of the solution vectors for points i and j 
// at time t given all the vital SPH data.
void sph_dynamics_compute(sph_dynamics_t* dyn, real_t t,
                          int i, int j,
                          real_t* Ui, real_t* Uj,
                          real_t Wi, real_t Wj,
                          vector_t* grad_Wi, vector_t* grad_Wj,
                          real_t* dUidt, real_t* dUjdt);

#endif
