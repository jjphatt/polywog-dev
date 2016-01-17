// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gc/gc.h"
#include "polywog/sph_kernel.h"

struct sph_kernel_t 
{
  char* name;
  real_t extent;
  void* context;
  void (*compute)(void* context, real_t eta, real_t det_H, real_t* W, real_t* dWdeta);
  void (*dtor)(void* context);
};

static void sph_kernel_free(void* ctx, void* dummy)
{
  sph_kernel_t* kernel = ctx;
  if ((kernel->dtor != NULL) && (kernel->context != NULL))
    kernel->dtor(kernel->context);
  string_free(kernel->name);
}

sph_kernel_t* sph_kernel_new(const char* name,
                             real_t extent,
                             void* context,
                             void (*compute)(void* context, real_t eta_mag, real_t det_H, real_t* W, real_t* dWdeta),
                             void (*dtor)(void* context))
{
  ASSERT(extent > 0.0);

  sph_kernel_t* kernel = GC_MALLOC(sizeof(sph_kernel_t));
  kernel->name = string_dup(name);
  kernel->extent = extent;
  kernel->context = context;
  kernel->compute = compute;
  kernel->dtor = dtor;
  GC_register_finalizer(kernel, sph_kernel_free, kernel, NULL, NULL);
  return kernel;
}

char* sph_kernel_name(sph_kernel_t* kernel)
{
  return kernel->name;
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

real_t sph_kernel_sum(sph_kernel_t* kernel, real_t n_per_h)
{
  ASSERT(n_per_h >= 0.0);

  // Self contribution.
  real_t sum, dWdeta;
  kernel->compute(kernel->context, 0.0, 1.0, &sum, &dWdeta);

  // Neighbors.
  if (n_per_h > 0.0)
  {
    int n_steps = (int)round(kernel->extent * n_per_h);
    vector_t eta;
    real_t deta = 1.0 / n_steps;
    for (int i = 0; i <= n_steps; ++i)
    {
      eta.x = i * deta;
      for (int j = 0; j <= n_steps; ++j)
      {
        eta.y = j * deta;
        for (int k = 0; k <= n_steps; ++k)
        {
          if ((i == 0) && (j == 0) && (k == 0)) 
            continue;
          eta.z = k * deta;
          real_t eta_mag = vector_mag(&eta);
          real_t W, dW;
          kernel->compute(kernel->context, eta_mag, 1.0, &W, &dW);
          sum += 2.0 * W;
        }
      }
    }
  }

  return sum;
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

typedef struct
{
  lookup1_t* W_table;
  lookup1_t* dW_table;
} tabular_t;

static void tabular_compute(void* context, real_t eta_mag, real_t det_H, real_t* W, real_t* dWdeta)
{
  tabular_t* tabular = context;
  *W = det_H * lookup1_value(tabular->W_table, eta_mag);
  *dWdeta = det_H * lookup1_value(tabular->dW_table, eta_mag);
}

static void tabular_free(void* context)
{
  tabular_t* tabular = context;
  lookup1_free(tabular->W_table);
  lookup1_free(tabular->dW_table);
  polymec_free(tabular);
}

sph_kernel_t* tabular_sph_kernel_new(sph_kernel_t* kernel, 
                                     lookup1_interpolation_t interpolation,
                                     int resolution)
{
  ASSERT(resolution > 1);
  tabular_t* tabular = polymec_malloc(sizeof(tabular_t));

  // Tabulate the (normalized) values and gradients of the kernel.
  real_t W_values[resolution], dW_values[resolution];
  real_t extent = sph_kernel_extent(kernel);
  real_t deta = extent / resolution;
  for (int i = 0; i < resolution; ++i)
    kernel->compute(kernel->context, i*deta, 1.0, &W_values[i], &dW_values[i]);
  tabular->W_table = lookup1_new(0.0, extent, resolution, W_values, interpolation);
  tabular->dW_table = lookup1_new(0.0, extent, resolution, dW_values, interpolation);

  // Now create our proper SPH kernel.
  int name_len = strlen(kernel->name);
  char name[name_len+128];
  snprintf(name, name_len+127, "table(%s)", kernel->name);
  return sph_kernel_new(name, extent, tabular, tabular_compute, tabular_free);
}

