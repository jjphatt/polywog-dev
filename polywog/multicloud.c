// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gc/gc.h"
#include "core/string_utils.h"
#include "core/array.h"
#include "polywog/multicloud.h"

struct multicloud_t 
{
  multicloud_iteration_t* iter;
};

struct multicloud_iteration_t 
{
  char* name;
  void* context;
  multicloud_iteration_vtable vtable;
};

struct multicloud_hierarchy_t 
{
  ptr_array_t* clouds;
};

static void multicloud_free(void* ctx, void* dummy)
{
  multicloud_t* mc = ctx;
  mc->iter = NULL;
}

multicloud_t* multicloud_new(multicloud_iteration_t* iteration_scheme)
{
  multicloud_t* mc = GC_MALLOC(sizeof(multicloud_t));
  mc->iter = iteration_scheme;
  GC_register_finalizer(mc, multicloud_free, mc, NULL, NULL);
  return mc;
}

multicloud_hierarchy_t* multicloud_coarsen(multicloud_t* mc,
                                           point_cloud_t* cloud)
{
  multicloud_hierarchy_t* h = multicloud_hierarchy_new(cloud);
  return h;
}

void multicloud_solve(multicloud_t* mc,
                      multicloud_hierarchy_t* hierarchy,
                      real_t* X,
                      int* num_iters,
                      real_t* res_norm)
{
}

static void multicloud_iteration_free(void* ctx, void* dummy)
{
  multicloud_iteration_t* iter = ctx;
  string_free(iter->name);
  if ((iter->context != NULL) && (iter->vtable.dtor != NULL))
    iter->vtable.dtor(iter->context);
}

multicloud_iteration_t* multicloud_iteration_new(const char* name,
                                                 void* context,
                                                 multicloud_iteration_vtable vtable)
{
  multicloud_iteration_t* iter = GC_MALLOC(sizeof(multicloud_iteration_t));
  iter->name = string_dup(name);
  iter->context = context;
  iter->vtable = vtable;
  GC_register_finalizer(iter, multicloud_iteration_free, iter, NULL, NULL);
  return iter;
}

multicloud_hierarchy_t* multicloud_hierarchy_new(point_cloud_t* cloud)
{
  multicloud_hierarchy_t* h = polymec_malloc(sizeof(multicloud_hierarchy_t));
  h->clouds = ptr_array_new();
  ptr_array_append(h->clouds, cloud);
  return h;
}

void multicloud_hierarchy_free(multicloud_hierarchy_t* hierarchy)
{
  ptr_array_free(hierarchy->clouds);
  polymec_free(hierarchy);
}

void multicloud_hierarchy_append(multicloud_hierarchy_t* hierarchy,
                                 point_cloud_t* coarser_cloud)
{
  ptr_array_append_with_dtor(hierarchy->clouds, coarser_cloud, DTOR(point_cloud_free));
}

