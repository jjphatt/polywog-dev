// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polywog/sph_dynamics.h"

struct sph_dynamics_t 
{
  char* name;
  void* context;
  void (*compute)(void* context, real_t t,
                  int i, int j,
                  real_t* Ui, real_t* Uj, 
                  real_t Wi, real_t Wj,
                  vector_t* grad_Wi, vector_t* grad_Wj,
                  real_t* dUidt, real_t* dUjdt);
  void (*dtor)(void* context);
};

sph_dynamics_t* sph_dynamics_new(const char* name,
                                 void* context,
                                 void (*compute)(void* context, real_t t,
                                                 int i, int j,
                                                 real_t* Ui, real_t* Uj, 
                                                 real_t Wi, real_t Wj,
                                                 vector_t* grad_Wi, vector_t* grad_Wj,
                                                 real_t* dUidt, real_t* dUjdt),
                                 void (*dtor)(void* context))
{
  ASSERT(compute != NULL);
  sph_dynamics_t* dyn = polymec_malloc(sizeof(sph_dynamics_t));
  dyn->name = string_dup(name);
  dyn->context = context;
  dyn->compute = compute;
  dyn->dtor = dtor;
  return dyn;
}

void sph_dynamics_free(sph_dynamics_t* dyn)
{
  if ((dyn->dtor != NULL) && (dyn->context != NULL))
    dyn->dtor(dyn->context);
  string_free(dyn->name);
  polymec_free(dyn);
}

void sph_dynamics_compute(sph_dynamics_t* dyn, real_t t,
                          int i, int j,
                          real_t* Ui, real_t* Uj,
                          real_t Wi, real_t Wj,
                          vector_t* grad_Wi, vector_t* grad_Wj,
                          real_t* dUidt, real_t* dUjdt)
{
  dyn->compute(dyn->context, t, i, j, Ui, Uj, Wi, Wj, grad_Wi, grad_Wj, dUidt, dUjdt);
}

