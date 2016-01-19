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

void sph_H_updater_update(sph_H_updater_t* updater, 
                          sym_tensor2_t* H, 
                          real_t zeroth_moment, 
                          sym_tensor2_t* second_moment,
                          sym_tensor2_t* new_H)
{
  // If the zeroth moment is "zeroish", we recommend doubling the number of 
  // neighbors per h.
  real_t s = 1.0;
  if (zeroth_moment < 1e-5)
  {
    s = 2.0;
  }
  else
  {
    // Look up the n_per_h value corresponding to the zeroth moment and compare 
    // it to our target.

    // Have we fallen off the table?
    real_t min_nh, max_nh;
    lookup1_get_bounds(updater->table, &min_nh, &max_nh);
    if (zeroth_moment > max_nh)
      s = 0.5; // cut in half!
    else
    {
      // Look up the value and compute our ratio of target to actual nh.
      // We allow values of this ratio between between 1/4 and 4.
      real_t nh = lookup1_value(updater->table, zeroth_moment);
      real_t target_nh = updater->n_per_h;
      s = MIN(4.0, MAX(0.25, target_nh / nh));
    }
  }

  // Compute the new determinant for H, following Thakar et al (2000). 
  real_t a = (s <= 1.0) ? 0.4 * (1.0 + s*s) 
                        : 0.4 * (1.0 + 1.0/(s*s*s));
  real_t det_H = sym_tensor2_det(H);
  real_t new_det_H = det_H / pow(1.0 - a + a*s, 3.0);
  
  sym_tensor2_set_identity(new_H, 1.0);
  if (updater->anisotropic)
  {
    // We compute the "direction" of the tensor if we care about anisotropy.

    // Find the minimum eigenvalue of the second moment.
    real_t lambdas[3];
    sym_tensor2_get_eigenvalues(second_moment, lambdas);
    real_t lambda_min = MIN(lambdas[0], MIN(lambdas[1], lambdas[2]));

    // Compute a weighting that articulates the importance of the second
    // moment as a function of s.
    real_t weight = MAX(0.0, MIN(1.0, 2.0/s - 1.0));
    if ((weight > 0.0) && 
        (sym_tensor2_det(second_moment) > 0.0) && 
        (lambda_min > 0.0))
    {
      // Compute the normalized "psi" tensor.
      real_t max_second_moment_elem = 
        MAX(second_moment->xx, 
            MAX(second_moment->xy, 
                MAX(second_moment->xz, 
                    MAX(second_moment->yy, 
                        MAX(second_moment->yz, second_moment->zz)))));
      sym_tensor2_t psi;
      sym_tensor2_copy(second_moment, &psi);
      sym_tensor2_scale(&psi, 1.0/max_second_moment_elem);
      if (sym_tensor2_det(&psi) < 1e-10)
        sym_tensor2_set_identity(&psi, 1.0);

      // The new shape of H is just the inverse of the psi tensor.
      sym_tensor2_invert(&psi, new_H);
    }
  }
  sym_tensor2_scale(new_H, (real_t)pow(new_det_H, 1.0/3.0));
}

