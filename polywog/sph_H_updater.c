// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "polywog/sph_H_updater.h"

struct sph_H_updater_t 
{
  sph_kernel_t* W;
  real_t n_per_h;
  int max_iters;
  real_t frac_change;
  bool anisotropic;

  // Lookup table mapping kernel sums to n_per_h values.
  lookup1_t* table;
};

static sph_H_updater_t* sph_H_updater_new(sph_kernel_t* W,
                                          real_t n_per_h,
                                          bool anisotropic)
{
  ASSERT(n_per_h > 0.0);

  sph_H_updater_t* updater = polymec_malloc(sizeof(sph_H_updater_t));
  updater->W = W;
  updater->n_per_h = n_per_h;
  updater->max_iters = 100;
  updater->frac_change = 0.05;
  updater->anisotropic = anisotropic;

  // Generate a lookup table.
  int table_size = 500;
  real_t min_nh = 0.0, max_nh = 10.0;
  real_t sumWijs[table_size];
  for (int i = 0; i < table_size; ++i)
  {
    real_t nh = i * (max_nh - min_nh) / table_size;
    // Our sums for kernels exclude the self contribution.
    sumWijs[i] = sph_kernel_sum(W, nh) - sph_kernel_sum(W, 0.0);
  }
  updater->table = lookup1_new(min_nh, max_nh, table_size, sumWijs, LOOKUP1_LINEAR);

  return updater;
}

sph_H_updater_t* isotropic_sph_H_updater_new(sph_kernel_t* W,
                                             real_t n_per_h)
{
  return sph_H_updater_new(W, n_per_h, false);
}

sph_H_updater_t* anisotropic_sph_H_updater_new(sph_kernel_t* W,
                                               real_t n_per_h)
{
  return sph_H_updater_new(W, n_per_h, true);
}

void sph_H_updater_free(sph_H_updater_t* updater)
{
  updater->W = NULL;
  lookup1_free(updater->table);
  polymec_free(updater);
}

static bool sph_H_updater_iterate(sph_H_updater_t* updater,
                                  point_t* x, point_t* ys, int num_neighbors,
                                  sym_tensor2_t* H, real_t* max_fractional_change)
{
  // Compute the normalized kernel sum at the point x.
  sym_tensor2_t Hi = {.xx = 1.0, .xy = 0.0, .xz = 0.0,
                                 .yy = 1.0, .yz = 0.0,
                                            .zz = 1.0};
  real_t sum;

  // Self contribution.
  vector_t zero = {0.0, 0.0, 0.0}, grad;
  sph_kernel_compute(updater->W, &zero, &Hi, &sum, &grad);

  // Neighbor contributions.
  for (int j = 0; j < num_neighbors; ++j)
  {
    real_t Wij;
    vector_t rij;
    point_displacement(x, &ys[j], &rij);
    sph_kernel_compute(updater->W, &rij, &Hi, &Wij, &grad);
    sum += Wij;
  }

  // Look up the n_per_h value corresponding to this sum and compare it 
  // to our target.

  // Have we fallen off the table?
  real_t min_nh, max_nh;
  lookup1_get_bounds(updater->table, &min_nh, &max_nh);
  if (sum > max_nh)
    return false;

  // Look up the value and compute our ratio of target to actual nh.
  real_t nh = lookup1_value(updater->table, sum);
  real_t target_nh = updater->n_per_h;
  real_t s = target_nh / nh;

  // Compute the new determinant for H, following Thakar et al (2000). 
  real_t a = (s <= 1.0) ? 0.4 * (1.0 + s*s) 
                        : 0.4 * (1.0 + 1.0/(s*s*s));
  real_t det_H = sym_tensor2_det(H);
  real_t new_det_H = det_H / pow(1.0 - a + a*s, 3.0);
  
  sym_tensor2_set_identity(H, 1.0);
  if (updater->anisotropic)
  {
    // We compute the "direction" of the tensor if we care about anisotropy.

    // First, compute the second moment.
    sym_tensor2_t X2;
    // FIXME

    // Find its minimum eigenvalue.
    real_t lambdas[3];
    sym_tensor2_get_eigenvalues(X2, lambdas);
    real_t lambda_min = MIN(lambdas[0], MIN(lambdas[1], lambdas[2]));

    // Compute a weighting that articulates a degree of confidence in 
    // the second moment calculation.
    real_t weight = MAX(0.0, MIN(1.0, 2.0/s - 1.0));
    if ((weight > 0.0) && 
        (sym_tensor2_det(X2) > 0.0) && 
        (lambda_min > 0.0))
    {
      // Compute the normalized "psi" tensor.

      // The new shape of H is just the inverse of the psi tensor.
      sym_tensor2_invert(psi, H);
    }
  }
  sym_tensor2_scale(H, (real_t)pow(new_det_H, 1.0/3.0));

  return true;
}

bool sph_H_updater_update(sph_H_updater_t* updater, 
                          point_t* x, point_t* ys, int num_neighbors,
                          sym_tensor2_t* H,
                          int* num_iterations,
                          real_t* max_fractional_change)
{
  bool failed = false;
  int num_iters = 0;
  while (!failed && *max_fractional_change > updater->frac_change)
  {
    if (!sph_H_updater_iterate(updater, x, ys, num_neighbors, H, max_fractional_change))
      failed = true;
    ++num_iters;
  }
  *num_iterations = num_iters;
  return !failed;
}

void sph_H_updater_set_max_iterations(sph_H_updater_t* updater, 
                                      int max_iters)
{
}

void sph_H_updater_set_convergence_threshold(sph_H_updater_t* updater, 
                                             real_t fractional_change)
{

}

