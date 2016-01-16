// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polywog/sph_kernel.h"

struct sph_kernel_t 
{
  char* name;
  real_t extent;
  void* context;
  void (*compute)(void* context, real_t eta, real_t det_H, real_t* W, real_t* dWdeta);
  void (*dtor)(void* context);
};

sph_kernel_t* sph_kernel_new(const char* name,
                             real_t extent,
                             void* context,
                             void (*compute)(void* context, real_t eta_mag, real_t det_H, real_t* W, real_t* dWdeta),
                             void (*dtor)(void* context))
{
  ASSERT(extent > 0.0);

  sph_kernel_t* kernel = polymec_malloc(sizeof(sph_kernel_t));
  kernel->name = string_dup(name);
  kernel->extent = extent;
  kernel->context = context;
  kernel->compute = compute;
  kernel->dtor = dtor;
  return kernel;
}

void sph_kernel_free(sph_kernel_t* kernel)
{
  if ((kernel->dtor != NULL) && (kernel->context != NULL))
    kernel->dtor(kernel->context);
  string_free(kernel->name);
  polymec_free(kernel);
}

real_t sph_kernel_extent(sph_kernel_t* kernel)
{
  return kernel->extent;
}

void sph_kernel_compute(sph_kernel_t* kernel, vector_t* x, sym_tensor2_t* H, real_t* W, vector_t* grad_W)
{
  // Compute eta = H o x.
  vector_t eta;
  sym_tensor2_dot_vector(H, x, &eta);
  real_t eta_mag = vector_mag(&eta);

  if (eta_mag <= kernel->extent)
  {
    // Compute W and dW/d|eta|.
    real_t det_H = sym_tensor2_det(H); 
    real_t dWdeta; 
    kernel->compute(kernel->context, eta_mag, det_H, W, &dWdeta);

    // Compute grad W.
    vector_t w = {.x = eta.x * dWdeta / eta_mag, 
                  .y = eta.y * dWdeta / eta_mag,
                  .z = eta.z * dWdeta / eta_mag};
    sym_tensor2_dot_vector(H, &w, grad_W);
  }
  else
  {
    *W = 0.0;
    grad_W->x = grad_W->y = grad_W->z = 0.0;
  }
}

static void b_spline_compute(void* context, real_t eta_mag, real_t det_H, real_t* W, real_t* dWdeta)
{
  if (eta_mag <= 1.0)
  {
    real_t eta2 = eta_mag * eta_mag;
    *W = (1.0 - 1.5*eta2 + 0.75*eta2*eta_mag) / M_PI;
    *dWdeta = (-1.5 + 2.25*eta2) / M_PI;
  }
  else if (eta_mag <= 2.0)
  {
    real_t term = 2.0 - eta_mag;
    *W = 0.25 * term * term * term / M_PI;
    *dWdeta = -0.75 * term * term / M_PI;
  }
  else
  {
    *W = *dWdeta = 0.0;
  }
}

sph_kernel_t* b_spline_sph_kernel_new()
{
  return sph_kernel_new("B-spline", 2.0, NULL, b_spline_compute, NULL);
}

